# =============================================================================
#  loh_chain_analysis.R
#  Chain-based LOH event caller for ChimeraMapR
#
#  Integrates with the existing app.R pipeline. Call after compute_loh_map()
#  and compute_peak_pairs() have both completed.
#
#  Public entry point:
#    run_chain_analysis(loh_segments, fused_peaks, peak_pairs, rt_df,
#                       chr_span, params)
#      -> list(chains, events, unclaimed_peaks, unclaimed_loh, event_table)
#
#  Schema version: "1.0"
#  To add new per-token information (e.g. ploidy, CN state), add it to the
#  $meta list inside make_token() and bump CHAIN_SCHEMA_VERSION.
# =============================================================================

CHAIN_SCHEMA_VERSION <- "1.0"

# %||% may already be defined by app.R; define it here only if absent so the
# module works standalone (e.g. in tests or sourced independently).
if (!exists("%||%", mode = "function", inherits = TRUE)) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# -----------------------------------------------------------------------------
#  TUNABLE PARAMETERS
#  Edit here; passed as a named list through the pipeline.
# -----------------------------------------------------------------------------
default_chain_params <- function() {
  list(
    # Distance from pos=1 or chr_len within which an LOH is "telomeric"
    tel_tol_bp      = 5000L,

    # Max NA-only gap between two same-state fixed runs to merge them
    merge_gap_bp    = 5000L,

    # A peak is "associated" with an LOH token if peak_start..peak_end
    # overlaps it, OR snp_pos is within peak_pad_bp of the token boundary.
    # Chimeric-read peaks are sharp; 200 bp is a tight but appropriate default.
    peak_pad_bp     = 200L,

    # Minimum spanning reads to make a read-based call; below this ->
    # AMBIGUOUS(low_coverage)
    min_span        = 3L,

    # Fraction of spanning reads sharing the return pattern for NCO_GC
    homog_frac      = 0.80,

    # A peak is "roughly centered" if |peak_mid - token_mid| <=
    # center_tol * token_length
    center_tol      = 0.25,

    # Interstitial token is "small" relative to flanking fixed regions
    # when length < small_frac * sum(flank lengths)
    small_frac      = 0.10,

    # Depth ratio (real coverage / baseline, from compute_coverage_map();
    # falls back to an SNP-density proxy if no coverage map is supplied)
    # below which a terminal region is called TERMINAL_DELETION (R01) rather
    # than CO_TERM/TERMINAL_LOH. ~0.5 is expected for loss of one parental
    # homolog; 0.60 leaves margin above that for noise.
    depth_drop      = 0.60,

    # Fixed-LOH tokens with fewer SNPs than this threshold may be too short
    # for any individual chimeric read to span with enough consecutive
    # same-haplotype calls to generate a peak.  Such tokens are classified
    # as POSSIBLE_GC rather than UNCATEGORIZED_LOH.
    # Should be set to the same value as the Minimum Run Length parameter
    # used during chimeric-read detection (default 2).
    min_snps_for_peak = 2L
  )
}

# =============================================================================
#  TOKEN CONSTRUCTORS
# =============================================================================

# make_token() — the fundamental unit. All fields except type/state/start/end
# are optional. New information layers go into meta.
make_token <- function(type,           # "F"=fixed | "H"=HET | "G"=NA gap | "TEL"
                       state    = NA,  # "REF_fixed" | "ALT_fixed" | "HET" | NA
                       start    = NA_integer_,
                       end      = NA_integer_,
                       n_snps   = NA_integer_,
                       depth_ratio  = NA_real_,  # token median depth / chrom median
                       peak_over    = NULL,   # peak/pair whose window spans this token
                       peak_left    = NULL,   # junction peak entering this token
                       peak_right   = NULL,   # junction peak leaving this token
                       bridged_gap  = FALSE,  # TRUE if a sub-merge_gap NA gap was absorbed
                       meta         = list()  # extensible — add ploidy/CN/etc here
) {
  structure(
    list(
      type        = type,
      state       = state,
      start       = as.integer(start),
      end         = as.integer(end),
      length_bp   = if (!is.na(start) && !is.na(end)) as.integer(end - start) else NA_integer_,
      n_snps      = as.integer(n_snps),
      depth_ratio = as.numeric(depth_ratio),
      peak_over   = peak_over,
      peak_left   = peak_left,
      peak_right  = peak_right,
      bridged_gap = bridged_gap,
      meta        = meta
    ),
    class = "loh_token"
  )
}

make_chain <- function(tokens, chrom, chr_len,
                       snp_start      = NA_integer_,
                       snp_end        = NA_integer_,
                       schema_version = CHAIN_SCHEMA_VERSION) {
  structure(
    list(tokens         = tokens,
         chrom          = chrom,
         chr_len        = as.integer(chr_len),
         snp_start      = as.integer(snp_start),
         snp_end        = as.integer(snp_end),
         schema_version = schema_version),
    class = "loh_chain"
  )
}

print.loh_token <- function(x, ...) {
  cat(sprintf("[%s/%s %d-%d snps=%s dr=%.2f pk_l=%s pk_over=%s pk_r=%s]\n",
    x$type,
    ifelse(is.na(x$state), "NA", x$state),
    ifelse(is.na(x$start), -1L, x$start),
    ifelse(is.na(x$end),   -1L, x$end),
    ifelse(is.na(x$n_snps), "?", as.character(x$n_snps)),
    ifelse(is.na(x$depth_ratio), 0, x$depth_ratio),
    if (!is.null(x$peak_left))  x$peak_left$edge_type  else "-",
    if (!is.null(x$peak_over))  x$peak_over$edge_type  else "-",
    if (!is.null(x$peak_right)) x$peak_right$edge_type else "-"
  ))
  invisible(x)
}

# =============================================================================
#  STEP 1 — BUILD RAW CHAINS FROM loh_segments
#  One chain per chromosome; NA gaps and telomere sentinels are inserted.
# =============================================================================

# Length-weighted average depth_ratio from compute_coverage_map()'s
# coverage_segments overlapping [seg_start, seg_end] on chr_name. Returns NA
# when no coverage map was supplied, or it has no segment for this span.
.lookup_depth_ratio <- function(coverage_segments, chr_name, seg_start, seg_end) {
  if (is.null(coverage_segments) || nrow(coverage_segments) == 0) return(NA_real_)
  cs <- coverage_segments[as.character(chrom) == chr_name]
  cs <- cs[end >= seg_start & start <= seg_end]
  if (nrow(cs) == 0) return(NA_real_)
  ov_len <- pmax(pmin(cs$end, seg_end) - pmax(cs$start, seg_start) + 1L, 1L)
  sum(cs$depth_ratio * ov_len) / sum(ov_len)
}

# Mean real read depth (n_total) over [seg_start, seg_end] on chr_name, read
# directly from compute_coverage_map()$coverage_table. Returns NA if no
# coverage table was supplied or no position falls in range (e.g. a "G" gap
# token with no defined SNPs at all).
.lookup_mean_depth <- function(coverage_table, chr_name, seg_start, seg_end) {
  if (is.null(coverage_table) || nrow(coverage_table) == 0) return(NA_real_)
  vals <- coverage_table[as.character(chrom) == chr_name &
                          pos >= seg_start & pos <= seg_end, n_total]
  if (length(vals) == 0) return(NA_real_)
  mean(vals)
}

build_raw_chains <- function(loh_segments, chr_span, params,
                             fused_peaks = NULL, peak_pairs = NULL,
                             snp_peaks   = NULL, coverage_segments = NULL,
                             coverage_table = NULL) {

  if (is.null(loh_segments) || nrow(loh_segments) == 0)
    return(list())

  segs <- copy(loh_segments)
  segs[, chrom := as.character(chrom)]

  chr_dt <- copy(chr_span)
  chr_dt[, chrom := as.character(chrom)]

  # If fused_peaks is not available, promote snp_peaks to serve in its place.
  # snp_peaks has snp_pos, peak_start, peak_end, and haplotype_label — enough
  # for the chain to attach and classify peaks without requiring fusion to have run.
  # We synthesise the columns .get_chr_peaks expects.
  peaks_source <- if (!is.null(fused_peaks) && nrow(fused_peaks) > 0) {
    fused_peaks
  } else if (!is.null(snp_peaks) && nrow(snp_peaks) > 0) {
    sp <- copy(snp_peaks)
    sp[, chrom := as.character(chrom)]
    # Synthesise required columns if absent
    if (!"fusion_group"   %in% names(sp)) sp[, fusion_group   := .I]
    if (!"fused_pos_bp"   %in% names(sp)) sp[, fused_pos_bp   := as.integer(snp_pos)]
    if (!"fused_start_bp" %in% names(sp)) sp[, fused_start_bp := as.integer(peak_start)]
    if (!"fused_end_bp"   %in% names(sp)) sp[, fused_end_bp   := as.integer(peak_end)]
    if (!"n_sub_peaks"    %in% names(sp)) sp[, n_sub_peaks    := 1L]
    # haplotype_label IS the edge type for singleton peaks
    if (!"best_edge_type" %in% names(sp) && "haplotype_label" %in% names(sp))
      sp[, best_edge_type := haplotype_label]
    sp
  } else {
    NULL
  }

  segs <- copy(loh_segments)
  segs[, chrom := as.character(chrom)]

  chr_dt <- copy(chr_span)
  chr_dt[, chrom := as.character(chrom)]

  # Compute chromosome-level median depth from loh_segments balance_mean
  # (proxy for depth; a proper depth column is used if added via meta later)
  chrom_median_depth <- segs[
    loh_state %in% c("REF_fixed", "ALT_fixed", "HET"),
    .(med_depth = median(n_snps / pmax(length_bp, 1L) * 1000, na.rm = TRUE)),
    by = chrom
  ]

  chains <- lapply(unique(segs$chrom), function(chr_name) {

    chr_segs <- segs[chrom == chr_name][order(start)]
    chr_len  <- chr_dt[chrom == chr_name, length]
    if (length(chr_len) == 0) chr_len <- max(chr_segs$end) + 1L
    chr_len  <- as.integer(chr_len)

    med_d <- chrom_median_depth[chrom == chr_name, med_depth]
    if (length(med_d) == 0 || is.na(med_d)) med_d <- 1

    # Build peak lookup for this chromosome
    pk_list <- .get_chr_peaks(peaks_source, peak_pairs, chr_name)

    # -----------------------------------------------------------
    # Build token list: TEL_L, segments, gaps, TEL_R
    # -----------------------------------------------------------
    # Define telomere sentinels by the min/max SNP positions for this
    # chromosome rather than by physical chromosome ends.  This absorbs the
    # large unscored regions that always flank the called segments (regions
    # that have been trimmed from the data) into the TEL tokens so that the
    # first and last called segments are directly adjacent to their sentinel —
    # no spurious G gap between TEL and the terminal F.
    first_pos <- as.integer(min(chr_segs$start))
    last_pos  <- as.integer(max(chr_segs$end))

    left_tel_end <- max(1L, first_pos - 1L)

    tokens <- list()
    tokens[[1]] <- make_token("TEL", state = NA,
                               start = 1L, end = left_tel_end, n_snps = 0L)

    prev_end <- left_tel_end
    for (si in seq_len(nrow(chr_segs))) {
      seg <- chr_segs[si]

      # Insert explicit NA gap token if there's an unscored region
      if (seg$start > prev_end + 1L) {
        gap_tok <- make_token(
          type  = "G",
          state = NA,
          start = as.integer(prev_end + 1L),
          end   = as.integer(seg$start - 1L),
          n_snps = 0L
        )
        tokens[[length(tokens) + 1L]] <- gap_tok
      }

      # Depth ratio: prefer the real coverage-derived ratio from
      # compute_coverage_map() (passed in as `coverage_segments`) — it reads
      # actual per-position read depth, modelled with the same EM+HMM
      # approach as compute_loh_map() but on total depth instead of allele
      # balance.  Fall back to the SNP-density proxy below only when no
      # coverage map is supplied (e.g. tests that build chains by hand):
      # n_snps/length_bp tracks local SNP-definition density, not real
      # depth, and can be high or low in a deleted region purely by chance
      # of where the VCF happens to define SNPs.
      real_dr <- .lookup_depth_ratio(coverage_segments, chr_name, seg$start, seg$end)

      seg_density <- if (seg$length_bp > 0)
        (seg$n_snps / seg$length_bp * 1000) / med_d else NA_real_

      tok_depth_ratio <- if (!is.na(real_dr)) real_dr else seg_density

      tok_type <- if (seg$loh_state %in% c("REF_fixed", "ALT_fixed")) "F"
                  else if (seg$loh_state == "HET") "H"
                  else "G"

      tok <- make_token(
        type        = tok_type,
        state       = seg$loh_state,
        start       = as.integer(seg$start),
        end         = as.integer(seg$end),
        n_snps      = as.integer(seg$n_snps),
        depth_ratio = tok_depth_ratio
      )

      # Attach any peaks associated with this token
      tok <- .attach_peaks(tok, pk_list, params$peak_pad_bp)

      tokens[[length(tokens) + 1L]] <- tok
      prev_end <- seg$end
    }

    # Right TEL absorbs the trailing unscored region — no trailing G token.
    right_tel_start <- min(as.integer(chr_len), last_pos + 1L)
    tokens[[length(tokens) + 1L]] <- make_token(
      "TEL", state = NA, start = right_tel_start, end = chr_len, n_snps = 0L
    )

    # ── Local terminal depth-ratio: F segment vs. its own adjacent flank ──────
    # A single chromosome-wide depth baseline doesn't separate a true
    # terminal/arm deletion from an ordinary region that simply sequences a
    # bit shallower for reasons unrelated to copy number — this chromosome's
    # depth varies by ~2x at large scale even outside any deletion. Comparing
    # a terminal F token to the H/G flank immediately bordering it controls
    # for that: ratio near 1 means depth is consistent across the junction
    # (the rule R02b path); ratio well below 1 means a real drop (R01).
    # Stored in meta (not the general depth_ratio field) so it's only used
    # by the terminal rules and never disturbs the chromosome-wide proxy
    # used elsewhere. G-typed flanks (no defined SNPs) are skipped in favor
    # of the nearest H, via .nearest_nonfixed_right()/.nearest_nonfixed_left().
    real_idx <- which(vapply(tokens, `[[`, character(1), "type") != "TEL")
    if (length(real_idx) > 0) {
      first_i <- real_idx[1]
      last_i  <- real_idx[length(real_idx)]

      annotate_terminal <- function(idx, nb_idx) {
        if (is.null(nb_idx) || tokens[[idx]]$type != "F") return(invisible(NULL))
        d_f  <- .lookup_mean_depth(coverage_table, chr_name,
                                   tokens[[idx]]$start, tokens[[idx]]$end)
        d_nb <- .lookup_mean_depth(coverage_table, chr_name,
                                   tokens[[nb_idx]]$start, tokens[[nb_idx]]$end)
        if (!is.na(d_f) && !is.na(d_nb) && d_nb > 0)
          tokens[[idx]]$meta$terminal_depth_ratio <<- d_f / d_nb
        invisible(NULL)
      }

      annotate_terminal(first_i, .nearest_nonfixed_right(tokens, first_i))
      if (last_i != first_i)
        annotate_terminal(last_i, .nearest_nonfixed_left(tokens, last_i))
    }

    make_chain(tokens, chr_name, chr_len,
               snp_start = first_pos, snp_end = last_pos)
  })

  names(chains) <- sapply(chains, `[[`, "chrom")
  chains
}

# Attach peaks to a token based on positional overlap
.attach_peaks <- function(tok, pk_list, peak_pad_bp) {
  if (is.null(pk_list) || length(pk_list) == 0) return(tok)
  if (is.na(tok$start) || is.na(tok$end)) return(tok)

  for (pk in pk_list) {
    pk_start <- pk$fused_start_bp %||% pk$peak_start
    pk_end   <- pk$fused_end_bp   %||% pk$peak_end
    snp_pos  <- pk$fused_pos_bp   %||% pk$snp_pos

    if (is.na(pk_start) || is.na(pk_end)) next

    # "over" = peak window overlaps token span (with pad)
    overlaps <- pk_start <= (tok$end   + peak_pad_bp) &&
                pk_end   >= (tok$start - peak_pad_bp)

    # "left junction" = snp_pos is just to the left of token start
    at_left  <- !is.na(snp_pos) &&
                snp_pos >= (tok$start - peak_pad_bp) &&
                snp_pos <   tok$start

    # "right junction" = snp_pos is just to the right of token end
    at_right <- !is.na(snp_pos) &&
                snp_pos >  tok$end &&
                snp_pos <= (tok$end + peak_pad_bp)

    if (overlaps)  tok$peak_over  <- pk
    if (at_left)   tok$peak_left  <- pk
    if (at_right)  tok$peak_right <- pk
  }
  tok
}

# Flatten fused_peaks and peak_pairs into a simple list of peak descriptors
# for a single chromosome
.get_chr_peaks <- function(fused_peaks, peak_pairs, chr_name) {
  if (is.null(fused_peaks)) return(list())
  fp <- copy(fused_peaks)
  fp[, chrom := as.character(chrom)]
  chr_fp <- fp[chrom == chr_name]
  if (nrow(chr_fp) == 0) return(list())

  # One entry per fusion group (representative row)
  grp_rep <- chr_fp[, .(
    fused_pos_bp    = fused_pos_bp[1],
    fused_start_bp  = fused_start_bp[1],
    fused_end_bp    = fused_end_bp[1],
    n_sub_peaks     = if ("n_sub_peaks" %in% names(chr_fp)) n_sub_peaks[1] else 1L,
    # best_edge_type is "singleton" for un-fused peaks — fall back to
    # haplotype_label which holds the per-read classification (gene_conversion,
    # binary, etc) set during Run Analysis.  haplotype_label is the ground truth
    # for singleton peaks; best_edge_type only adds information for fused groups.
    best_edge_type  = {
      bet <- if ("best_edge_type"  %in% names(chr_fp)) best_edge_type[1]  else NA_character_
      hal <- if ("haplotype_label" %in% names(chr_fp)) haplotype_label[1] else NA_character_
      if (is.na(bet) || bet == "singleton") hal else bet
    },
    haplotype_label = if ("haplotype_label" %in% names(chr_fp)) haplotype_label[1] else NA_character_,
    n_spanning      = 0L,
    jaccard         = NA_real_,
    pair_edge_type  = NA_character_,
    # Per-read switch count from classify_peak_haplotype(), computed at
    # Run Analysis time for "binary"/"internal_crossover" labelled peaks —
    # the only direct read evidence available for peaks that never get a
    # peak_pairs row (e.g. a terminal binary peak with no eligible partner
    # within range, or an internal_crossover peak excluded from pairing
    # entirely). Used below as a fallback for n_spanning.
    n_read_support  = if ("n_read_support" %in% names(chr_fp)) n_read_support[1] else NA_integer_
  ), by = fusion_group]

  # Overlay spanning reads, jaccard, and pair_edge_type from peak_pairs.
  # A pair is "associated" with a peak if the peak's position falls anywhere
  # within the pair's span (snp_pos_a .. snp_pos_b), i.e. it is one of the
  # two endpoints or sits between them.
  if (!is.null(peak_pairs) && nrow(peak_pairs) > 0) {
    pp <- copy(peak_pairs)
    pp[, chrom := as.character(chrom)]
    chr_pp <- pp[chrom == chr_name]
    if (nrow(chr_pp) > 0) {
      for (ri in seq_len(nrow(grp_rep))) {
        pos <- grp_rep$fused_pos_bp[ri]
        # Match pairs where this peak is one endpoint or between the two endpoints
        matching <- chr_pp[snp_pos_a <= pos & snp_pos_b >= pos]
        if (nrow(matching) > 0) {
          best_pair <- matching[which.max(n_spanning)]
          grp_rep$n_spanning[ri]     <- best_pair$n_spanning
          grp_rep$jaccard[ri]        <- best_pair$jaccard
          grp_rep$pair_edge_type[ri] <- best_pair$edge_type   # always write
          if (is.na(grp_rep$best_edge_type[ri]))
            grp_rep$best_edge_type[ri] <- best_pair$edge_type
        }
      }
    }
  }

  # Peaks with no peak_pairs evidence at all (n_spanning still at the 0L
  # default — no eligible partner found, or excluded from pairing entirely)
  # fall back to their own per-read switch count so a real terminal binary
  # peak or standalone internal_crossover peak isn't reported as 0 support
  # just because it never formed a pair.
  no_pair_evidence <- grp_rep$n_spanning == 0L & !is.na(grp_rep$n_read_support)
  grp_rep$n_spanning[no_pair_evidence] <- grp_rep$n_read_support[no_pair_evidence]

  lapply(seq_len(nrow(grp_rep)), function(i) as.list(grp_rep[i]))
}

# =============================================================================
#  STEP 2 — CANONICALISE
#  Apply rewrite rules to fixpoint. Modifies the token list in-place.
#  Returns the updated token list.
# =============================================================================

canonicalise <- function(chain, params) {
  tokens <- chain$tokens

  # Zero-SNP G tokens represent coordinate space between called segments (i.e.
  # unscored HET context in real loh_segments output that only emits fixed rows).
  # They must be kept so motif rules R04/R07 see the non-fixed flanking context
  # around interstitial fixed tracts.  Small zero-SNP gaps between same-state F
  # tokens are absorbed by the F-G-F merge rule below; large ones are left in
  # place and act as G (non-fixed) flanks for het_bounded / double_gc matching.
  # Terminal F tokens that genuinely start at position 1 are still directly
  # adjacent to the left TEL sentinel (no gap token exists when seg$start == 1),
  # so R01/R02 continue to match correctly without this filter.
  changed <- TRUE

  while (changed) {
    changed <- FALSE
    i <- 1L

    while (i <= length(tokens)) {

      # Rule 1: F G F (same state, gap < merge_gap) -> merge into single F
      if (i <= length(tokens) - 2L) {
        a <- tokens[[i]]
        g <- tokens[[i + 1L]]
        b <- tokens[[i + 2L]]

        if (a$type == "F" && g$type == "G" && b$type == "F" &&
            !is.na(a$state) && !is.na(b$state) &&
            a$state == b$state &&
            !is.na(g$end) && !is.na(g$start) &&
            (g$end - g$start) < params$merge_gap_bp) {

          merged <- make_token(
            type        = "F",
            state       = a$state,
            start       = a$start,
            end         = b$end,
            n_snps      = as.integer(sum(c(a$n_snps, b$n_snps), na.rm = TRUE)),
            depth_ratio = mean(c(a$depth_ratio, b$depth_ratio), na.rm = TRUE),
            peak_over   = a$peak_over %||% b$peak_over,
            peak_left   = a$peak_left,
            peak_right  = b$peak_right,
            bridged_gap = TRUE,
            meta        = c(a$meta, b$meta)
          )
          tokens <- c(tokens[seq_len(i - 1L)],
                      list(merged),
                      tokens[seq.int(i + 3L, length(tokens))])
          changed <- TRUE
          next  # re-examine from same position
        }
      }

      # Rule 2: H G H (same state HET gap) -> merge
      if (i <= length(tokens) - 2L) {
        a <- tokens[[i]]
        g <- tokens[[i + 1L]]
        b <- tokens[[i + 2L]]

        if (a$type == "H" && g$type == "G" && b$type == "H" &&
            !is.na(g$end) && !is.na(g$start) &&
            (g$end - g$start) < params$merge_gap_bp) {

          merged <- make_token(
            type        = "H",
            state       = "HET",
            start       = a$start,
            end         = b$end,
            n_snps      = as.integer(sum(c(a$n_snps, b$n_snps), na.rm = TRUE)),
            depth_ratio = mean(c(a$depth_ratio, b$depth_ratio), na.rm = TRUE),
            peak_left   = a$peak_left,
            peak_right  = b$peak_right,
            bridged_gap = TRUE
          )
          tokens <- c(tokens[seq_len(i - 1L)],
                      list(merged),
                      tokens[seq.int(i + 3L, length(tokens))])
          changed <- TRUE
          next
        }
      }

      # Rule 3: non-fixed (H or G) on both sides of a small G gap -> merge.
      # In real data, loh_segments produces G tokens for unscored stretches
      # that sit between or adjacent to called HET segments.  Without this,
      # a chain like G [F] G never matches the H [F] H motif rules.
      # Merged token is H if either side is H, otherwise G.
      if (i <= length(tokens) - 2L) {
        a <- tokens[[i]]
        g <- tokens[[i + 1L]]
        b <- tokens[[i + 2L]]

        if (.is_non_fixed(a) && g$type == "G" && .is_non_fixed(b) &&
            !is.na(g$end) && !is.na(g$start) &&
            (g$end - g$start) < params$merge_gap_bp) {

          merged_type  <- if (a$type == "H" || b$type == "H") "H" else "G"
          merged_state <- if (merged_type == "H") "HET" else NA_character_
          merged <- make_token(
            type        = merged_type,
            state       = merged_state,
            start       = a$start,
            end         = b$end,
            n_snps      = as.integer(sum(c(a$n_snps, b$n_snps), na.rm = TRUE)),
            depth_ratio = mean(c(a$depth_ratio, b$depth_ratio), na.rm = TRUE),
            peak_left   = a$peak_left,
            peak_right  = b$peak_right,
            bridged_gap = TRUE
          )
          tokens <- c(tokens[seq_len(i - 1L)],
                      list(merged),
                      tokens[seq.int(i + 3L, length(tokens))])
          changed <- TRUE
          next
        }
      }

      i <- i + 1L
    }
  }

  chain$tokens <- tokens
  chain
}

# =============================================================================
#  SUB-PROCEDURE: CLASSIFY_TRACT
#  Determines NCO_GC, CO_GC, or AMBIGUOUS from a peak descriptor.
#  Called from motif rules — reads edge_type from the peak, falls back to
#  raw read patterns if edge_type is unavailable.
# =============================================================================

classify_tract <- function(peak, params) {

  if (is.null(peak))
    return(list(call = "AMBIGUOUS", reason = "no_peak", n_support = 0L))

  edge_type  <- peak$best_edge_type %||% peak$edge_type
  n_raw      <- peak$n_spanning %||% NA_integer_
  n_spanning <- if (is.null(n_raw) || is.na(n_raw)) NA_integer_ else as.integer(n_raw)

  if (is.na(edge_type) || is.null(edge_type))
    return(list(call = "AMBIGUOUS", reason = "no_edge_type",
                n_support = if (is.na(n_spanning)) 0L else n_spanning))

  # Only apply the min_span coverage gate when we actually have a real
  # spanning-read count (n_spanning > 0).  NA or 0 means no peak_pairs row
  # was found for this peak — the edge_type was set by the app's per-peak read
  # classifier and is itself sufficient evidence.  A genuinely low-count
  # situation has n_spanning > 0 but below the threshold.
  has_count <- !is.na(n_spanning) && n_spanning > 0L
  if (has_count && n_spanning < params$min_span)
    return(list(call = "AMBIGUOUS", reason = "low_coverage",
                n_support = n_spanning))

  result <- switch(edge_type,
    "gene_conversion"    = list(call = "NCO_GC", reason = edge_type, n_support = n_spanning),
    "crossover"          = list(call = "CO_GC",  reason = edge_type, n_support = n_spanning),
    "internal_crossover" = list(call = "CO_GC",  reason = edge_type, n_support = n_spanning),
    "binary"             = list(call = "AMBIGUOUS", reason = "binary_single_peak",
                                 n_support = n_spanning),
    "independent_events" = list(call = "AMBIGUOUS", reason = "independent_reads",
                                 n_support = n_spanning),
    "ambiguous"       = list(call = "AMBIGUOUS", reason = "mixed_pattern",
                              n_support = n_spanning),
    "unresolvable"    = list(call = "AMBIGUOUS", reason = "no_covering_reads",
                              n_support = n_spanning),
    list(call = "AMBIGUOUS", reason = paste0("unknown_edge:", edge_type),
         n_support = n_spanning)
  )
  result
}

# Junction-level complementarity check for two binary peaks flanking a
# fixed tract — the two-binary-peak case described in our discussion.
# Returns NCO_GC_LARGE | CO_GC | AMBIGUOUS(reason)
classify_two_binary_junction <- function(left_peak, right_peak,
                                          interstitial_state, params) {
  if (is.null(left_peak) || is.null(right_peak))
    return(list(call = "AMBIGUOUS", reason = "missing_flanking_peak", n_support = 0L))

  et_l  <- left_peak$best_edge_type  %||% left_peak$edge_type
  et_r  <- right_peak$best_edge_type %||% right_peak$edge_type
  ns_l_raw <- left_peak$n_spanning  %||% NA_integer_
  ns_r_raw <- right_peak$n_spanning %||% NA_integer_
  ns_l  <- if (is.null(ns_l_raw) || is.na(ns_l_raw)) NA_integer_ else as.integer(ns_l_raw)
  ns_r  <- if (is.null(ns_r_raw) || is.na(ns_r_raw)) NA_integer_ else as.integer(ns_r_raw)
  # ns is NA when neither peak has a real count; only used for the coverage gate
  ns    <- if (!is.na(ns_l) || !is.na(ns_r))
    as.integer(sum(c(ns_l, ns_r), na.rm = TRUE))
  else
    NA_integer_

  # Both must be binary for this resolver to be relevant
  both_binary <- isTRUE(et_l == "binary") && isTRUE(et_r == "binary")
  if (!both_binary) {
    # Fall back to the single-tract classifier on whichever peak spans better
    best_pk <- if (isTRUE((ns_l %||% 0L) >= (ns_r %||% 0L))) left_peak else right_peak
    return(classify_tract(best_pk, params))
  }

  # Only gate on coverage when we actually have a real count
  has_count <- !is.na(ns) && ns > 0L
  if (has_count && ns < params$min_span)
    return(list(call = "AMBIGUOUS", reason = "low_coverage",
                n_support = ns))

  # Complementarity: left peak should be AAABBB, right peak BBBAAAA
  # We infer orientation from peak position relative to the interstitial state.
  # "Complementary" means: left peak switches INTO interstitial state,
  # right peak switches OUT of interstitial state.
  # Both switching the same direction -> NCO_GC_large_tract (copy-neutral LOH)
  # Complementary -> CO_GC

  # Use the pair's edge_type from peak_pairs if the pair spans the tract
  # (this is populated by compute_peak_pairs when loh_in_gap == TRUE)
  pair_et <- left_peak$pair_edge_type %||% right_peak$pair_edge_type

  if (!is.null(pair_et) && !is.na(pair_et)) {
    return(switch(pair_et,
      "crossover"       = list(call = "CO_GC",            reason = "complementary_binary", n_support = ns),
      "gene_conversion" = list(call = "NCO_GC_LARGE",     reason = "homogeneous_binary",   n_support = ns),
      list(call = "AMBIGUOUS", reason = paste0("pair_edge:", pair_et), n_support = ns)
    ))
  }

  # No pair-level edge_type: can't distinguish CO from NCO at the junction level
  list(call = "AMBIGUOUS", reason = "binary_no_pair", n_support = ns)
}

# =============================================================================
#  HELPER PREDICATES
# =============================================================================

# TEL proximity is measured relative to the SNP-bounded ends of the chain
# (stored as chain$snp_start / chain$snp_end by build_raw_chains).  Falls
# back to positional 1 / chr_len for manually constructed test chains.
.is_telomeric <- function(token, chain, params) {
  if (is.na(token$start) || is.na(token$end)) return(FALSE)
  ref_start <- if (!is.null(chain$snp_start) && !is.na(chain$snp_start))
    chain$snp_start else 1L
  ref_end   <- if (!is.null(chain$snp_end)   && !is.na(chain$snp_end))
    chain$snp_end   else chain$chr_len
  token$start <= ref_start + params$tel_tol_bp ||
  token$end   >= ref_end   - params$tel_tol_bp
}

.opp_state <- function(state) {
  if (is.na(state)) return(NA_character_)
  if (state == "REF_fixed") "ALT_fixed" else "REF_fixed"
}

.has_peak <- function(token) {
  !is.null(token$peak_over) || !is.null(token$peak_left) || !is.null(token$peak_right)
}

# Depth ratio used by the terminal-event rules (R01/R02/R02b): the local
# ratio vs. the token's own adjacent flank (set by build_raw_chains() from a
# coverage_table, comparing the F token to its nearest H neighbor) when
# available, else the token's general depth_ratio (chromosome-wide real
# ratio, or the SNP-density proxy if no coverage map was supplied). A
# single chromosome-wide baseline doesn't separate a true terminal/arm
# deletion from a region that is merely low-depth for reasons unrelated to
# copy number when regional depth varies a lot across the chromosome —
# comparing a terminal segment to its own neighbor controls for that.
.terminal_depth_ratio <- function(tok) {
  tok$meta$terminal_depth_ratio %||% tok$depth_ratio
}

# TRUE for any token that is not fixed LOH — H (called HET) and G (unscored
# gap) are both treated as non-fixed context for motif matching.  In real data,
# loh_segments only contains SNP-position rows; the coordinate stretches between
# them become G tokens, not H tokens, even when the underlying biology is HET.
.is_non_fixed <- function(token) {
  token$type %in% c("H", "G")
}

# Safe n_spanning extractor — returns 0L for NULL or NA
.ns <- function(peak) {
  if (is.null(peak)) return(0L)
  v <- peak$n_spanning %||% 0L
  if (is.na(v)) 0L else as.integer(v)
}

.best_peak <- function(token) {
  # Return the peak with the most spanning reads, or peak_over if all equal
  candidates <- Filter(Negate(is.null),
                       list(token$peak_over, token$peak_left, token$peak_right))
  if (length(candidates) == 0) return(NULL)
  spans <- sapply(candidates, .ns)
  candidates[[which.max(spans)]]
}

.is_roughly_centered <- function(token, center_tol) {
  pk <- token$peak_over
  if (is.null(pk)) return(FALSE)
  pk_mid  <- ((pk$fused_start_bp %||% pk$peak_start) +
               (pk$fused_end_bp  %||% pk$peak_end)) / 2
  tok_mid <- (token$start + token$end) / 2
  abs(pk_mid - tok_mid) <= center_tol * token$length_bp
}

.make_event <- function(call, chrom, tokens_involved, evidence_peaks = NULL,
                        n_support = NA_integer_, notes = "") {
  # Compute genomic span from the tokens involved (excluding TEL sentinels)
  real_toks <- Filter(function(t) t$type != "TEL", tokens_involved)
  ev_start  <- if (length(real_toks) > 0) min(sapply(real_toks, `[[`, "start"), na.rm = TRUE) else NA_integer_
  ev_end    <- if (length(real_toks) > 0) max(sapply(real_toks, `[[`, "end"),   na.rm = TRUE) else NA_integer_

  pk_edge   <- if (!is.null(evidence_peaks))
    paste(unique(sapply(evidence_peaks, function(p)
      p$best_edge_type %||% p$edge_type %||% "?")), collapse = ";")
  else NA_character_

  list(
    event_class     = call,
    chrom           = chrom,
    start           = as.integer(ev_start),
    end             = as.integer(ev_end),
    length_bp       = as.integer(ev_end - ev_start),
    n_support       = as.integer(n_support),
    peak_edge_types = pk_edge,
    notes           = notes,
    tokens          = tokens_involved   # keep for downstream use / plotting
  )
}

# =============================================================================
#  STEP 3 — MOTIF SCANNER
#  Single left-to-right pass. Rules in priority order.
#  Each rule is a list(id, match_fn, fire_fn) where:
#    match_fn(tokens, i, chain, params) -> NULL | list(span_indices, ...)
#    fire_fn(match_result, chain, params) -> list(event, rewrite_tokens | NULL)
# =============================================================================

# ── Rule builders ─────────────────────────────────────────────────────────────

# Rule 1: Terminal deletion — T [F] ... (depth low)
rule_terminal_deletion <- list(
  id = "R01_terminal_deletion",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    # Forward: TEL [F] ...
    if (i + 1L <= n &&
        tokens[[i]]$type == "TEL" &&
        tokens[[i + 1L]]$type == "F" &&
        !is.na(.terminal_depth_ratio(tokens[[i + 1L]])) &&
        .terminal_depth_ratio(tokens[[i + 1L]]) < params$depth_drop)
      return(list(span = c(i, i + 1L), direction = "fwd",
                  f_tok = tokens[[i + 1L]]))

    # Reverse: ... [F] TEL
    if (i >= 2L &&
        tokens[[i]]$type == "F" &&
        tokens[[i + 1L]]$type == "TEL" &&
        !is.na(.terminal_depth_ratio(tokens[[i]])) &&
        .terminal_depth_ratio(tokens[[i]]) < params$depth_drop)
      return(list(span = c(i, i + 1L), direction = "rev",
                  f_tok = tokens[[i]]))

    NULL
  },
  fire_fn = function(m, chain, params) {
    ev <- .make_event("TERMINAL_DELETION", chain$chrom,
                      list(m$f_tok),
                      notes = sprintf("depth_ratio=%.2f < %.2f",
                                      .terminal_depth_ratio(m$f_tok), params$depth_drop))
    list(event = ev, rewrite = NULL, claims = list(peak = NULL, loh = m$f_tok))
  }
)

# Rule 2: Terminal LOH — T [F●] H  (depth ok, peak at internal end, orientation ok)
rule_terminal_loh <- list(
  id = "R02_terminal_loh",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)

    check_tl <- function(tel_i, f_i, h_i, direction) {
      tel <- tokens[[tel_i]]
      ftk <- tokens[[f_i]]
      htk <- tokens[[h_i]]
      if (tel$type != "TEL" || ftk$type != "F" || !.is_non_fixed(htk)) return(NULL)
      if (is.na(.terminal_depth_ratio(ftk)) ||
          .terminal_depth_ratio(ftk) < params$depth_drop) return(NULL)   # depth ok check
      pk <- ftk$peak_right %||% ftk$peak_left %||% ftk$peak_over
      if (is.null(pk)) return(NULL)
      list(span = c(tel_i, f_i, h_i), f_tok = ftk, peak = pk,
           direction = direction)
    }

    # Forward: TEL [F●] H
    if (i + 2L <= n) {
      m <- check_tl(i, i + 1L, i + 2L, "fwd")
      if (!is.null(m)) return(m)
    }
    # Reverse: H [●F] TEL
    if (i + 2L <= n) {
      m <- check_tl(i + 2L, i + 1L, i, "rev")
      if (!is.null(m)) return(m)
    }
    NULL
  },
  fire_fn = function(m, chain, params) {
    edge_type  <- m$peak$best_edge_type %||% m$peak$edge_type
    n_raw      <- m$peak$n_spanning %||% NA_integer_
    n_spanning <- if (is.null(n_raw) || is.na(n_raw)) NA_integer_ else as.integer(n_raw)
    has_count  <- !is.na(n_spanning) && n_spanning > 0L

    # Terminal LOH has a single haplotype boundary: the expected peak pattern
    # is "binary" (one clean switch, e.g. REF→ALT).  gene_conversion and
    # crossover patterns imply internal events, not a terminal boundary.
    call <- if (has_count && n_spanning < params$min_span) {
      "AMBIGUOUS(low_coverage)"
    } else if (!is.null(edge_type) && !is.na(edge_type) && edge_type == "binary") {
      "TERMINAL_LOH"
    } else {
      paste0("AMBIGUOUS(terminal_non_binary_peak:", edge_type %||% "NA", ")")
    }

    ev <- .make_event(call, chain$chrom, list(m$f_tok),
                      evidence_peaks = list(m$peak),
                      n_support = n_spanning,
                      notes = paste0("terminal; edge=", edge_type %||% "NA"))
    list(event = ev, rewrite = NULL,
         claims = list(peak = m$peak, loh = m$f_tok))
  }
)

# Rule 2b: Terminal, no peak — secondary depth-based assessment.
#
# A chimeric-read peak can only form where a read carries BOTH alleles in
# cis and switches between them. A true terminal/arm deletion removes the
# second haplotype's sequence entirely beyond the breakpoint, so no read can
# ever show that switch there — R02 (rule_terminal_loh) is structurally
# unable to fire for a genuine deletion, peak or no peak. When a terminal F
# token has no peak at all, the only remaining evidence is depth: a real
# drop (e.g. ~50% from losing one parental homolog) means a true
# TERMINAL_DELETION; depth that is consistent with the rest of the
# chromosome means the peak is simply missing — most likely a gap in SNP
# coverage swallowed the junction — and this is a terminal LOH/crossover
# (CO_TERM) that should go to manual review rather than being silently
# reported as a deletion or dropped as UNCATEGORIZED_LOH.
rule_terminal_no_peak <- list(
  id = "R02b_terminal_no_peak",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)

    check_tnp <- function(tel_i, f_i, h_i, direction) {
      tel <- tokens[[tel_i]]
      ftk <- tokens[[f_i]]
      htk <- tokens[[h_i]]
      if (tel$type != "TEL" || ftk$type != "F" || !.is_non_fixed(htk)) return(NULL)
      if (.has_peak(ftk)) return(NULL)                                       # R02 handles this case
      if (is.na(.terminal_depth_ratio(ftk))) return(NULL)                     # no depth evidence at all
      if (.terminal_depth_ratio(ftk) < params$depth_drop) return(NULL)        # R01 handles real deletions
      list(span = c(tel_i, f_i, h_i), f_tok = ftk, direction = direction)
    }

    # Forward: TEL [F] H
    if (i + 2L <= n) {
      m <- check_tnp(i, i + 1L, i + 2L, "fwd")
      if (!is.null(m)) return(m)
    }
    # Reverse: H [F] TEL
    if (i + 2L <= n) {
      m <- check_tnp(i + 2L, i + 1L, i, "rev")
      if (!is.null(m)) return(m)
    }
    NULL
  },
  fire_fn = function(m, chain, params) {
    ev <- .make_event("CO_TERM", chain$chrom, list(m$f_tok),
                      notes = sprintf(
                        paste0("no_chimera_peak; depth_ratio=%.2f >= %.2f; ",
                               "terminal LOH/crossover inferred from depth ",
                               "consistent with its adjacent flank ",
                               "(likely a SNP-coverage gap masking the junction peak)"),
                        .terminal_depth_ratio(m$f_tok), params$depth_drop))
    list(event = ev, rewrite = NULL, claims = list(peak = NULL, loh = m$f_tok))
  }
)

# Rule 3: TCO-captured-TCO — H ●[F]● F̄→TEL (opposite fixed reaches telomere,
#         peaks on both sides)
rule_tco_captured_tco <- list(
  id = "R03_tco_captured_tco",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    # Pattern: H [F] F̄{tel}  (with peaks on both junctions of F)
    # Positions:  i  i+1  i+2 (i+2 reaches telomere)
    if (i + 2L > n) return(NULL)
    htk  <- tokens[[i]]
    ftk  <- tokens[[i + 1L]]
    fbar <- tokens[[i + 2L]]
    if (!.is_non_fixed(htk) || ftk$type != "F" || fbar$type != "F") return(NULL)
    if (is.na(ftk$state) || is.na(fbar$state)) return(NULL)
    if (fbar$state != .opp_state(ftk$state)) return(NULL)
    if (!.is_telomeric(fbar, chain, params)) return(NULL)
    pk_l <- ftk$peak_left  %||% htk$peak_right
    pk_r <- ftk$peak_right %||% fbar$peak_left
    if (is.null(pk_l) || is.null(pk_r)) return(NULL)
    list(span = c(i, i + 2L), f_tok = ftk, fbar_tok = fbar,
         h_tok = htk, pk_l = pk_l, pk_r = pk_r)
  },
  fire_fn = function(m, chain, params) {
    ev <- .make_event("TCO_CAPTURED_TCO", chain$chrom,
                      list(m$h_tok, m$f_tok, m$fbar_tok),
                      evidence_peaks = list(m$pk_l, m$pk_r),
                      n_support = as.integer(.ns(m$pk_l) + .ns(m$pk_r)),
                      notes = "terminal CO over earlier terminal CO")
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk_l, m$pk_r),
                       loh  = list(m$f_tok)))
  }
)

# Rule 4 / 5: H ●[F]● H — interstitial HET-bounded conversion tract
rule_het_bounded <- list(
  id = "R04_het_bounded",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    if (i + 2L > n) return(NULL)
    l <- tokens[[i]]; ftk <- tokens[[i + 1L]]; r <- tokens[[i + 2L]]
    if (!.is_non_fixed(l) || ftk$type != "F" || !.is_non_fixed(r)) return(NULL)
    # Junction peaks can sit on the F token directly, or on the neighbouring
    # H tokens (the peak snp_pos is just outside the F boundaries)
    pk_l <- ftk$peak_left  %||% l$peak_right  %||% l$peak_over
    pk_r <- ftk$peak_right %||% r$peak_left   %||% r$peak_over
    pk   <- ftk$peak_over  %||% pk_l %||% pk_r
    if (is.null(pk)) return(NULL)
    list(span = c(i, i + 2L), f_tok = ftk, l_tok = l, r_tok = r,
         pk = pk, pk_l = pk_l, pk_r = pk_r)
  },
  fire_fn = function(m, chain, params) {
    # Single-peak case
    if (!is.null(m$pk) &&
        (is.null(m$pk_l) || identical(m$pk, m$pk_l)) &&
        (is.null(m$pk_r) || identical(m$pk, m$pk_r))) {
      tract <- classify_tract(m$pk, params)
      call  <- switch(tract$call,
        NCO_GC  = "NCO_GC",
        CO_GC   = "CO_GC",
        paste0("AMBIGUOUS(", tract$reason, ")")
      )
      ev <- .make_event(call, chain$chrom,
                        list(m$l_tok, m$f_tok, m$r_tok),
                        evidence_peaks = list(m$pk),
                        n_support = tract$n_support)
    } else {
      # Two-peak (binary) case: junction-level complementarity check
      tract <- classify_two_binary_junction(m$pk_l, m$pk_r,
                                             m$f_tok$state, params)
      call  <- switch(tract$call,
        NCO_GC       = "NCO_GC",
        NCO_GC_LARGE = "NCO_GC_LARGE",
        CO_GC        = "CO_GC",
        paste0("AMBIGUOUS(", tract$reason, ")")
      )
      ev <- .make_event(call, chain$chrom,
                        list(m$l_tok, m$f_tok, m$r_tok),
                        evidence_peaks = Filter(Negate(is.null),
                                                list(m$pk_l, m$pk_r)),
                        n_support = tract$n_support)
    }
    list(event = ev, rewrite = NULL,
         claims = list(peak = Filter(Negate(is.null),
                                     list(m$pk, m$pk_l, m$pk_r)),
                       loh  = list(m$f_tok)))
  }
)

# Rule 6: F ●[F̄]● F — small interstitial opposite-state fragment
#          (GC from a prior event captured by terminal CO)
rule_opp_sandwich <- list(
  id = "R06_opp_sandwich",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    if (i + 2L > n) return(NULL)
    l <- tokens[[i]]; ftk <- tokens[[i + 1L]]; r <- tokens[[i + 2L]]
    if (l$type != "F" || ftk$type != "F" || r$type != "F") return(NULL)
    if (is.na(l$state) || is.na(ftk$state) || is.na(r$state)) return(NULL)
    if (ftk$state == l$state || l$state != r$state) return(NULL)
    if (ftk$state != .opp_state(l$state)) return(NULL)
    # Check "small" criterion
    combined_flank <- (l$length_bp %||% 0L) + (r$length_bp %||% 0L)
    if ((ftk$length_bp %||% Inf) >= params$small_frac * combined_flank) return(NULL)
    pk_l <- ftk$peak_left  %||% l$peak_right  %||% l$peak_over
    pk_r <- ftk$peak_right %||% r$peak_left   %||% r$peak_over
    pk   <- ftk$peak_over  %||% pk_l %||% pk_r
    if (is.null(pk)) return(NULL)
    list(span = c(i, i + 2L), f_tok = ftk, l_tok = l, r_tok = r,
         pk = pk, pk_l = pk_l, pk_r = pk_r)
  },
  fire_fn = function(m, chain, params) {
    centered <- .is_roughly_centered(m$f_tok, params$center_tol) ||
                .is_roughly_centered(
                  make_token("F", start = m$f_tok$start, end = m$f_tok$end,
                             peak_over = m$pk_l), params$center_tol)

    tract <- if (!is.null(m$pk_l) && !is.null(m$pk_r))
      classify_two_binary_junction(m$pk_l, m$pk_r, m$f_tok$state, params)
    else
      classify_tract(m$pk, params)

    call <- switch(tract$call,
      NCO_GC       = "NCO_GC_in_terminal",
      NCO_GC_LARGE = "NCO_GC_in_terminal",
      CO_GC        = "CO_GC",
      paste0("AMBIGUOUS(", tract$reason, ")")
    )

    # Rewrite: combine flanking F tokens into one (they become adjacent after
    # consuming the tiny F̄ fragment)
    merged_flank <- make_token(
      type        = "F",
      state       = m$l_tok$state,
      start       = m$l_tok$start,
      end         = m$r_tok$end,
      n_snps      = as.integer(sum(c(m$l_tok$n_snps, m$r_tok$n_snps), na.rm = TRUE)),
      depth_ratio = mean(c(m$l_tok$depth_ratio, m$r_tok$depth_ratio), na.rm = TRUE),
      peak_left   = m$l_tok$peak_left,
      peak_right  = m$r_tok$peak_right,
      bridged_gap = FALSE,
      meta        = c(m$l_tok$meta, m$r_tok$meta)
    )

    # Record the merged token in the event so reconcile() can match it in
    # the post-rewrite chain (the originals no longer exist after the rewrite)
    ev <- .make_event(call, chain$chrom,
                      list(merged_flank, m$f_tok),
                      evidence_peaks = Filter(Negate(is.null),
                                              list(m$pk, m$pk_l, m$pk_r)),
                      n_support = tract$n_support,
                      notes = paste0("small_opp_fragment; centered=", centered))

    list(event = ev,
         rewrite = list(span = m$span, replacement = list(merged_flank)),
         claims  = list(peak = Filter(Negate(is.null), list(m$pk, m$pk_l, m$pk_r)),
                        loh  = list(m$f_tok)))
  }
)

# Rule 7: H [F] H [F̄] H — double GC
rule_double_gc <- list(
  id = "R07_double_gc",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    if (i + 4L > n) return(NULL)
    h1  <- tokens[[i]];     f1  <- tokens[[i + 1L]]
    h2  <- tokens[[i + 2L]]; f2  <- tokens[[i + 3L]]
    h3  <- tokens[[i + 4L]]
    if (!.is_non_fixed(h1) || f1$type != "F" || !.is_non_fixed(h2) ||
        f2$type != "F" || !.is_non_fixed(h3)) return(NULL)
    if (is.na(f1$state) || is.na(f2$state)) return(NULL)
    if (f1$state == f2$state) return(NULL)   # same state = not independent GCs
    pk1 <- f1$peak_over %||% f1$peak_left  %||% h1$peak_right %||%
                             f1$peak_right %||% h2$peak_left
    pk2 <- f2$peak_over %||% f2$peak_left  %||% h2$peak_right %||%
                             f2$peak_right %||% h3$peak_left
    if (is.null(pk1) || is.null(pk2)) return(NULL)
    list(span = c(i, i + 4L), f1_tok = f1, f2_tok = f2,
         h_toks = list(h1, h2, h3), pk1 = pk1, pk2 = pk2)
  },
  fire_fn = function(m, chain, params) {
    t1 <- classify_tract(m$pk1, params)
    t2 <- classify_tract(m$pk2, params)
    call <- if (t1$call == "NCO_GC" && t2$call == "NCO_GC") "DOUBLE_GC"
            else paste0("AMBIGUOUS(double_gc_unresolved:",
                        t1$reason, ";", t2$reason, ")")
    ev <- .make_event(call, chain$chrom,
                      c(m$h_toks, list(m$f1_tok, m$f2_tok)),
                      evidence_peaks = list(m$pk1, m$pk2),
                      n_support = as.integer(.ns(m$pk1) + .ns(m$pk2)),
                      notes = paste("f1:", t1$reason, "f2:", t2$reason))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk1, m$pk2),
                       loh  = list(m$f1_tok, m$f2_tok)))
  }
)

# Rule 8: F ●[ ]● F̄ (opposite fixed meeting at peak, no fixed token between)
#          -> CROSSOVER_NO_TRACT (clean crossover, tract below LOH resolution)
rule_crossover_no_tract <- list(
  id = "R08_crossover_no_tract",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    if (i + 1L > n) return(NULL)
    l <- tokens[[i]]; r <- tokens[[i + 1L]]
    # Directly adjacent F and F̄ (after canonicalisation, G between them is
    # below merge_gap and would have been merged; anything left is real)
    if (l$type != "F" || r$type != "F") return(NULL)
    if (is.na(l$state) || is.na(r$state)) return(NULL)
    if (l$state != .opp_state(r$state)) return(NULL)
    # Need a peak at the junction
    pk <- l$peak_right %||% r$peak_left
    if (is.null(pk)) return(NULL)
    et <- pk$best_edge_type %||% pk$edge_type
    if (is.null(et) || !et %in% c("crossover", "binary")) return(NULL)
    list(span = c(i, i + 1L), l_tok = l, r_tok = r, pk = pk)
  },
  fire_fn = function(m, chain, params) {
    ev <- .make_event("CROSSOVER_NO_TRACT", chain$chrom,
                      list(m$l_tok, m$r_tok),
                      evidence_peaks = list(m$pk),
                      n_support = .ns(m$pk),
                      notes = "adjacent_opposite_fixed; no_tract")
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk), loh = list()))
  }
)

# Rule 9: H ●[ ]● H (peak, no fixed token between) — sub-resolution tract
rule_subres_tract <- list(
  id = "R09_subres_tract",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)
    if (i + 1L > n) return(NULL)
    l <- tokens[[i]]; r <- tokens[[i + 1L]]
    if (!.is_non_fixed(l) || !.is_non_fixed(r)) return(NULL)
    pk <- l$peak_right %||% r$peak_left
    if (is.null(pk)) return(NULL)
    list(span = c(i, i + 1L), l_tok = l, r_tok = r, pk = pk)
  },
  fire_fn = function(m, chain, params) {
    tract <- classify_tract(m$pk, params)
    call  <- switch(tract$call,
      NCO_GC = "NCO_GC_subres",
      CO_GC  = "CO_GC_subres",
      paste0("AMBIGUOUS(", tract$reason, ")")
    )
    ev <- .make_event(call, chain$chrom, list(m$l_tok, m$r_tok),
                      evidence_peaks = list(m$pk),
                      n_support = tract$n_support,
                      notes = "no_fixed_tract; peak_only")
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk), loh = list()))
  }
)

# =============================================================================
#  INTERSTITIAL RULES (new peak model)
# =============================================================================

# ── Helpers: skip G tokens to reach nearest non-fixed context ─────────────────

# Walk left through G tokens to find the first non-fixed (H or G) token that
# has meaningful content.  G tokens skip over because build_raw_chains never
# calls .attach_peaks on them — peaks always land on H or F.  Stops (returns
# NULL) when TEL or another F is reached before any H is found.
.nearest_nonfixed_left <- function(tokens, i) {
  j <- i - 1L
  while (j >= 1L && tokens[[j]]$type == "G") j <- j - 1L
  if (j >= 1L && .is_non_fixed(tokens[[j]])) j else NULL
}

.nearest_nonfixed_right <- function(tokens, i) {
  n <- length(tokens)
  j <- i + 1L
  while (j <= n && tokens[[j]]$type == "G") j <- j + 1L
  if (j <= n && .is_non_fixed(tokens[[j]])) j else NULL
}

# Left junction peak for an F token: the peak that marks the H→F transition.
# Priority: F$peak_left > left_ctx$peak_right > left_ctx$peak_over (left of F)
#           > F$peak_over (if its snp_pos is left of F.start)
.left_junction_peak <- function(f_tok, left_ctx) {
  if (!is.null(f_tok$peak_left)) return(f_tok$peak_left)
  if (!is.null(left_ctx$peak_right)) return(left_ctx$peak_right)
  pk <- left_ctx$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    if (!is.na(pos) && pos < f_tok$start) return(pk)
  }
  pk <- f_tok$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    if (!is.na(pos) && pos < f_tok$start) return(pk)
  }
  NULL
}

# Right junction peak for an F token: the peak that marks the F→H transition.
# Priority: F$peak_right > F$peak_over (right of F.end) > right_ctx$peak_left
#           > right_ctx$peak_over (right of F.end)
.right_junction_peak <- function(f_tok, right_ctx) {
  if (!is.null(f_tok$peak_right)) return(f_tok$peak_right)
  pk <- f_tok$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    if (!is.na(pos) && pos > f_tok$end) return(pk)
  }
  if (!is.null(right_ctx$peak_left)) return(right_ctx$peak_left)
  pk <- right_ctx$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    if (!is.na(pos) && pos > f_tok$end) return(pk)
  }
  NULL
}

# ── Rule R10: Direct classification from peak edge type ───────────────────────
# Fires when an F token has a gene_conversion, crossover, or internal_crossover
# peak directly associated with it.  These peak types are self-classifying:
# gene_conversion → NCO_GC; crossover / internal_crossover → CO_GC.
# Uses classify_tract() to preserve the min_span coverage gate.
rule_peak_direct <- list(
  id = "R10_peak_direct",
  match_fn = function(tokens, i, chain, params) {
    tok <- tokens[[i]]
    if (tok$type != "F") return(NULL)
    pk <- .best_peak(tok)
    if (is.null(pk)) return(NULL)
    et <- pk$best_edge_type %||% pk$edge_type
    if (is.null(et) || is.na(et)) return(NULL)
    if (!et %in% c("gene_conversion", "crossover", "internal_crossover")) return(NULL)
    list(span = c(i, i), f_tok = tok, pk = pk)
  },
  fire_fn = function(m, chain, params) {
    tract <- classify_tract(m$pk, params)
    call  <- switch(tract$call,
      NCO_GC = "NCO_GC",
      CO_GC  = "CO_GC",
      paste0("AMBIGUOUS(", tract$reason, ")")
    )
    ev <- .make_event(call, chain$chrom, list(m$f_tok),
                      evidence_peaks = list(m$pk),
                      n_support = tract$n_support,
                      notes = paste0("peak_type=",
                                     m$pk$best_edge_type %||% m$pk$edge_type))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk), loh = m$f_tok))
  }
)

# ── Rule R11: Two binary flanking peaks (H-[F]-H, G-transparent) ──────────────
# Fires when an interstitial F is bracketed by a binary peak at each junction.
# G gaps between H and F are treated as transparent (the rule skips them).
# The left/right binary pair is classified via classify_two_binary_junction(),
# which reads pair_edge_type from compute_peak_pairs() to distinguish NCO/CO.
rule_two_binary_flanking <- list(
  id = "R11_two_binary_flanking",
  match_fn = function(tokens, i, chain, params) {
    tok <- tokens[[i]]
    if (tok$type != "F") return(NULL)

    li <- .nearest_nonfixed_left(tokens, i)
    ri <- .nearest_nonfixed_right(tokens, i)
    if (is.null(li) || is.null(ri)) return(NULL)

    left_ctx  <- tokens[[li]]
    right_ctx <- tokens[[ri]]

    pk_l <- .left_junction_peak(tok, left_ctx)
    pk_r <- .right_junction_peak(tok, right_ctx)
    if (is.null(pk_l) || is.null(pk_r)) return(NULL)

    # Both junction peaks must be binary
    et_l <- pk_l$best_edge_type %||% pk_l$edge_type
    et_r <- pk_r$best_edge_type %||% pk_r$edge_type
    if (!isTRUE(et_l == "binary") || !isTRUE(et_r == "binary")) return(NULL)

    # Guard: don't count the same peak twice
    pos_l <- pk_l$fused_pos_bp %||% pk_l$snp_pos
    pos_r <- pk_r$fused_pos_bp %||% pk_r$snp_pos
    if (!is.na(pos_l) && !is.na(pos_r) && pos_l == pos_r) return(NULL)

    list(span = c(i, ri), f_tok = tok, l_tok = left_ctx, r_tok = right_ctx,
         pk_l = pk_l, pk_r = pk_r)
  },
  fire_fn = function(m, chain, params) {
    tract <- classify_two_binary_junction(m$pk_l, m$pk_r, m$f_tok$state, params)
    call  <- switch(tract$call,
      NCO_GC       = "NCO_GC",
      NCO_GC_LARGE = "NCO_GC",
      CO_GC        = "CO_GC",
      paste0("AMBIGUOUS(", tract$reason, ")")
    )
    ev <- .make_event(call, chain$chrom,
                      list(m$l_tok, m$f_tok, m$r_tok),
                      evidence_peaks = list(m$pk_l, m$pk_r),
                      n_support = tract$n_support,
                      notes = paste0("two_binary_flanking; state=", m$f_tok$state))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk_l, m$pk_r), loh = list(m$f_tok)))
  }
)

# ── Rule R12: LOH-crossover (large interstitial LOH flanked by chimeric peaks) ─
# Fires when an F token with a large LOH is flanked by two junction peaks that
# share a pair record with edge_type = "crossover" from the LOH-crossover probing
# in compute_peak_pairs().  This is the fallback when R10 cannot fire because the
# flanking peaks' windows do not overlap the F token directly.
rule_loh_crossover <- list(
  id = "R12_loh_crossover",
  match_fn = function(tokens, i, chain, params) {
    tok <- tokens[[i]]
    if (tok$type != "F") return(NULL)

    li <- .nearest_nonfixed_left(tokens, i)
    ri <- .nearest_nonfixed_right(tokens, i)
    if (is.null(li) || is.null(ri)) return(NULL)

    left_ctx  <- tokens[[li]]
    right_ctx <- tokens[[ri]]

    pk_l <- .left_junction_peak(tok, left_ctx)
    pk_r <- .right_junction_peak(tok, right_ctx)
    if (is.null(pk_l) || is.null(pk_r)) return(NULL)

    # Both junction peaks must share a pair_edge_type = "crossover" record.
    pair_et_l <- pk_l$pair_edge_type
    pair_et_r <- pk_r$pair_edge_type
    if (!isTRUE(pair_et_l == "crossover") && !isTRUE(pair_et_r == "crossover"))
      return(NULL)

    # Guard: don't match the same peak on both sides.
    pos_l <- pk_l$fused_pos_bp %||% pk_l$snp_pos
    pos_r <- pk_r$fused_pos_bp %||% pk_r$snp_pos
    if (!is.na(pos_l) && !is.na(pos_r) && pos_l == pos_r) return(NULL)

    list(span = c(i, ri), f_tok = tok, l_tok = left_ctx, r_tok = right_ctx,
         pk_l = pk_l, pk_r = pk_r)
  },
  fire_fn = function(m, chain, params) {
    # n_spanning: both peaks share the same pair record, so avoid double-counting.
    ns_raw <- m$pk_l$n_spanning %||% m$pk_r$n_spanning %||% NA_integer_
    ns     <- if (is.null(ns_raw) || is.na(ns_raw)) NA_integer_ else as.integer(ns_raw)

    # Spanning a large LOH requires exceptional read length; accept n >= 2
    # (half of the default min_span floor of 3) because the fixed LOH allele
    # on the far side supplies the complementary haplotype evidence.
    loh_min_span <- max(2L, params$min_span - 1L)
    has_count <- !is.na(ns) && ns > 0L
    call <- if (has_count && ns < loh_min_span)
      paste0("AMBIGUOUS(low_coverage)")
    else
      "CO_GC"

    ev <- .make_event(call, chain$chrom,
                      list(m$l_tok, m$f_tok, m$r_tok),
                      evidence_peaks = list(m$pk_l, m$pk_r),
                      n_support = ns,
                      notes = paste0("loh_crossover; state=", m$f_tok$state))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk_l, m$pk_r), loh = list(m$f_tok)))
  }
)

# The ordered rule set — priority highest to lowest.
# R10, R11, R12 are the interstitial rules. Old R03-R09 remain disabled
# pending further redesign (R06 opp_sandwich, R07 double_gc, R08/R09).
MOTIF_RULES <- list(
  rule_terminal_deletion,       # R01
  rule_terminal_loh,            # R02
  rule_terminal_no_peak,        # R02b — terminal F with no peak: depth decides CO_TERM vs deletion
  rule_loh_crossover,           # R12 — crossover through large interstitial LOH (before R10)
  rule_peak_direct,             # R10 — gene_conversion / crossover / internal_crossover
  rule_two_binary_flanking      # R11 — two binary peaks flanking H-[F]-H
  # rule_tco_captured_tco,      # R03 — disabled
  # rule_het_bounded,           # R04 — disabled (replaced by R10/R11)
  # rule_opp_sandwich,          # R06 — disabled
  # rule_double_gc,             # R07 — disabled
  # rule_crossover_no_tract,    # R08 — disabled
  # rule_subres_tract           # R09 — disabled
)

# =============================================================================
#  STEP 4 — SCAN A SINGLE CHAIN
# =============================================================================

scan_chain <- function(chain, params, rules = MOTIF_RULES) {
  tokens  <- chain$tokens
  events  <- list()
  # Track claimed peaks by their position (proxy for identity)
  claimed_peak_pos <- numeric(0)

  i <- 1L
  while (i <= length(tokens)) {
    hit <- NULL

    for (rule in rules) {
      m <- rule$match_fn(tokens, i, chain, params)
      if (!is.null(m)) {
        # Check peak has not already been claimed
        pk_pos <- .peak_pos_from_match(m)
        if (any(pk_pos %in% claimed_peak_pos)) { m <- NULL; next }
        hit <- list(match = m, rule = rule)
        break
      }
    }

    if (!is.null(hit)) {
      fired  <- hit$rule$fire_fn(hit$match, chain, params)
      events <- c(events, list(fired$event))

      # Claim peaks
      new_pk_pos <- .peak_pos_from_match(hit$match)
      claimed_peak_pos <- unique(c(claimed_peak_pos, new_pk_pos))

      if (!is.null(fired$rewrite)) {
        span <- fired$rewrite$span
        repl <- fired$rewrite$replacement
        tokens <- c(tokens[seq_len(span[1] - 1L)],
                    repl,
                    if (span[2] < length(tokens))
                      tokens[seq.int(span[2] + 1L, length(tokens))]
                    else list())
        # Re-canonicalise after rewrite
        chain$tokens <- tokens
        chain        <- canonicalise(chain, params)
        tokens       <- chain$tokens
        i            <- max(1L, span[1] - 1L)   # back up to check new neighbours
      } else {
        i <- hit$match$span[2] + 1L
      }
    } else {
      i <- i + 1L
    }
  }

  list(chain  = chain,
       events = events,
       claimed_peak_pos = claimed_peak_pos)
}

.peak_pos_from_match <- function(m) {
  # Extract snp/fused positions from any peak fields in a match result
  pks <- Filter(Negate(is.null),
                list(m$pk %||% NULL, m$pk_l %||% NULL, m$pk_r %||% NULL,
                     m$pk1 %||% NULL, m$pk2 %||% NULL,
                     m$peak %||% NULL))
  if (length(pks) == 0) return(numeric(0))
  unlist(lapply(pks, function(p) {
    v <- p$fused_pos_bp %||% p$snp_pos %||% NA_real_
    if (is.na(v)) numeric(0) else v
  }))
}

# =============================================================================
#  STEP 5 — RECONCILE (Pass 3)
#  After scanning all chains, collect unclaimed LOH tokens and peaks.
# =============================================================================

reconcile <- function(scan_results, chains, fused_peaks, peak_pairs,
                      snp_peaks = NULL, params = default_chain_params()) {

  all_events <- unlist(lapply(scan_results, `[[`, "events"), recursive = FALSE)

  # Find unclaimed LOH segments
  claimed_ranges <- lapply(all_events, function(ev) {
    toks <- ev$tokens
    fixed_toks <- Filter(function(t) !is.null(t) && t$type == "F", toks)
    lapply(fixed_toks, function(t) list(chrom = ev$chrom, start = t$start, end = t$end))
  })
  claimed_ranges <- unlist(claimed_ranges, recursive = FALSE)

  unclaimed_loh <- list()
  for (cname in names(chains)) {
    chain <- chains[[cname]]
    for (tok in chain$tokens) {
      if (tok$type != "F") next
      tok_mid <- (tok$start + tok$end) / 2
      # A token is claimed if any event's claimed range contains its midpoint.
      # Using midpoint overlap (rather than exact start/end equality) handles
      # cases where a rewrite merged tokens and the event records the merged span.
      already_claimed <- any(sapply(claimed_ranges, function(cr) {
        cr$chrom == cname &&
          !is.na(cr$start) && !is.na(cr$end) &&
          cr$start <= tok_mid && cr$end >= tok_mid
      }))
      if (!already_claimed) {
        unclaimed_loh <- c(unclaimed_loh, list(list(
          chrom = cname, start = tok$start, end = tok$end,
          state = tok$state, n_snps = tok$n_snps
        )))
      }
    }
  }

  # Find unclaimed peaks — check fused_peaks first, fall back to snp_peaks
  all_claimed_pos <- unlist(lapply(scan_results, `[[`, "claimed_peak_pos"))

  peak_source_for_reconcile <- if (!is.null(fused_peaks) && nrow(fused_peaks) > 0) {
    fp <- copy(fused_peaks)
    fp[, chrom := as.character(chrom)]
    fp[, pos_col := fused_pos_bp]
    # best_edge_type is NA for peak classes excluded from pairing entirely
    # (e.g. "internal_crossover" — see FUSION_HEURISTICS$excluded_peak_classes
    # in compute_peak_pairs()), since get_peak_edge_info() never finds a pair
    # record for them. Fall back to haplotype_label, the per-peak read
    # classification, which IS the ground truth for these peaks — mirrors
    # the same fallback .get_chr_peaks() applies when building chains.
    fp[, et_col  := {
      bet <- if ("best_edge_type"  %in% names(fp)) best_edge_type else NA_character_
      hal <- if ("haplotype_label" %in% names(fp)) haplotype_label else NA_character_
      ifelse(is.na(bet) | bet == "singleton", hal, bet)
    }]
    fp[, ns_col := if ("n_read_support" %in% names(fp)) n_read_support else NA_integer_]
    fp[, .(chrom, pos_col, et_col, ns_col)]
  } else if (!is.null(snp_peaks) && nrow(snp_peaks) > 0) {
    sp <- copy(snp_peaks)
    sp[, chrom := as.character(chrom)]
    sp[, pos_col := as.integer(snp_pos)]
    sp[, et_col  := if ("haplotype_label" %in% names(sp)) haplotype_label else NA_character_]
    sp[, ns_col  := if ("n_read_support" %in% names(sp)) n_read_support else NA_integer_]
    sp[, .(chrom, pos_col, et_col, ns_col)]
  } else {
    NULL
  }

  unclaimed_peaks <- list()
  if (!is.null(peak_source_for_reconcile)) {
    for (ri in seq_len(nrow(peak_source_for_reconcile))) {
      pos <- peak_source_for_reconcile$pos_col[ri]
      if (!is.na(pos) && !pos %in% all_claimed_pos) {
        unclaimed_peaks <- c(unclaimed_peaks, list(list(
          chrom          = peak_source_for_reconcile$chrom[ri],
          snp_pos        = pos,
          edge_type      = peak_source_for_reconcile$et_col[ri],
          n_read_support = peak_source_for_reconcile$ns_col[ri]
        )))
      }
    }
  }

  # Emit UNCATEGORIZED / POSSIBLE_GC events for anything left over.
  # Fixed tokens with very few SNPs (< min_snps_for_peak) may be too short
  # for chimeric reads to generate a detectable peak even when a real gene
  # conversion occurred.  Flag these separately so the user knows to
  # consider lowering the Minimum Run Length parameter.
  uncat_events <- lapply(unclaimed_loh, function(u) {
    n <- u$n_snps %||% NA_integer_
    is_small <- !is.na(n) && n < params$min_snps_for_peak
    if (is_small) {
      list(event_class = "POSSIBLE_GC", chrom = u$chrom,
           start = u$start, end = u$end, length_bp = u$end - u$start,
           n_support = NA_integer_, peak_edge_types = NA_character_,
           notes = paste0("state=", u$state, "; n_snps=", n,
                          " < min_snps_for_peak (", params$min_snps_for_peak,
                          "); consider lowering the Minimum Run Length parameter"),
           tokens = list())
    } else {
      list(event_class = "UNCATEGORIZED_LOH", chrom = u$chrom,
           start = u$start, end = u$end, length_bp = u$end - u$start,
           n_support = NA_integer_, peak_edge_types = NA_character_,
           notes = paste0("state=", u$state), tokens = list())
    }
  })
  # Peaks whose own edge_type is self-classifying (gene_conversion /
  # crossover / internal_crossover) but that never attached to any token —
  # e.g. an internal_crossover sitting in a SNP-desert with no flanking LOH
  # to anchor a motif rule on — are promoted here using the same
  # classify_tract() logic the motif rules use, rather than left as an
  # unclassified UNCATEGORIZED_PEAK. Such peaks are excluded from
  # compute_peak_pairs()'s pairing step entirely (no peak_pairs row, so no
  # pairwise n_spanning) and have no flanking-LOH corroboration, so their
  # only direct read evidence is the per-read switch count computed at
  # classification time (n_read_support, from classify_peak_haplotype()).
  # They're named "_subres" (sub-resolution) and fall to "review" confidence
  # in build_event_table(), distinguishing them from the better-supported
  # CO_GC/NCO_GC calls fired by the motif rules.
  uncat_peak_events <- lapply(unclaimed_peaks, function(u) {
    tract <- classify_tract(list(best_edge_type = u$edge_type,
                                 n_spanning = u$n_read_support), params)
    if (tract$call %in% c("NCO_GC", "CO_GC")) {
      list(event_class = paste0(tract$call, "_subres"), chrom = u$chrom,
           start = as.integer(u$snp_pos), end = as.integer(u$snp_pos),
           length_bp = 0L, n_support = tract$n_support,
           peak_edge_types = u$edge_type %||% NA_character_,
           notes = "no_fixed_tract; peak_only", tokens = list())
    } else {
      list(event_class = "UNCATEGORIZED_PEAK", chrom = u$chrom,
           start = as.integer(u$snp_pos), end = as.integer(u$snp_pos),
           length_bp = 0L, n_support = NA_integer_,
           peak_edge_types = u$edge_type %||% NA_character_,
           notes = "", tokens = list())
    }
  })

  list(
    events          = c(all_events, uncat_events, uncat_peak_events),
    unclaimed_loh   = unclaimed_loh,
    unclaimed_peaks = unclaimed_peaks
  )
}

# =============================================================================
#  BUILD EVENT TABLE  (flat data.table for display/export)
# =============================================================================

build_event_table <- function(events) {
  if (length(events) == 0)
    return(data.table(
      event_class = character(), chrom = character(),
      start = integer(), end = integer(), length_kb = numeric(),
      n_support = integer(), peak_edge_types = character(),
      confidence = character(), notes = character()
    ))

  rows <- lapply(events, function(ev) {
    confidence <- if (ev$event_class %in% c("NCO_GC", "CO_GC", "TERMINAL_LOH",
                                             "CROSSOVER_NO_TRACT", "DOUBLE_GC",
                                             "TCO_CAPTURED_TCO"))
      "high"
    else
      "review"

    data.table(
      event_class     = ev$event_class,
      chrom           = ev$chrom,
      start           = as.integer(ev$start),
      end             = as.integer(ev$end),
      length_kb       = round((as.integer(ev$end) - as.integer(ev$start)) / 1000, 3),
      n_support       = as.integer(ev$n_support),
      peak_edge_types = ev$peak_edge_types %||% NA_character_,
      confidence      = confidence,
      notes           = ev$notes %||% ""
    )
  })

  out <- rbindlist(rows, fill = TRUE)
  setorder(out, chrom, start)
  out
}

# =============================================================================
#  PUBLIC ENTRY POINT
# =============================================================================

#' run_chain_analysis()
#'
#' @param loh_segments   data.table from compute_loh_map()$loh_segments
#' @param fused_peaks    data.table from compute_peak_pairs()$fused_peaks
#' @param peak_pairs     data.table from compute_peak_pairs()$peak_pairs
#' @param rt_df          chimeric read table from run_chimera_analysis()
#' @param chr_span       data.table from run_chimera_analysis()$chr_span
#' @param coverage_segments data.table from compute_coverage_map()$coverage_segments;
#'                       chromosome-wide depth_ratio (general-purpose fallback).
#' @param coverage_table data.table from compute_coverage_map()$coverage_table;
#'                       per-position real depth, used to compute each
#'                       terminal F token's depth ratio against its own
#'                       adjacent flank — what R01/R02b use to tell a true
#'                       terminal deletion apart from a terminal LOH/
#'                       crossover that is simply missing its junction peak.
#'                       NULL falls back to coverage_segments/SNP-density.
#' @param params         named list; defaults from default_chain_params()
#'
#' @return named list:
#'   $chains          per-chromosome loh_chain objects (canonicalised)
#'   $events          flat list of event objects
#'   $unclaimed_peaks list of unmatched peaks
#'   $unclaimed_loh   list of unmatched LOH segments
#'   $event_table     data.table ready for display / download

run_chain_analysis <- function(loh_segments,
                                fused_peaks  = NULL,
                                peak_pairs   = NULL,
                                snp_peaks    = NULL,
                                rt_df        = NULL,
                                chr_span,
                                coverage_segments = NULL,
                                coverage_table     = NULL,
                                params       = default_chain_params()) {

  message("  [chain] Building raw chains ...")
  chains <- build_raw_chains(loh_segments, chr_span, params,
                              fused_peaks, peak_pairs, snp_peaks,
                              coverage_segments, coverage_table)

  message("  [chain] Canonicalising chains ...")
  chains <- lapply(chains, canonicalise, params = params)

  message("  [chain] Scanning for motifs ...")
  scan_results <- lapply(names(chains), function(cname) {
    scan_chain(chains[[cname]], params)
  })
  names(scan_results) <- names(chains)

  # Update chains with any rewrites from the scan
  for (cname in names(chains)) {
    chains[[cname]] <- scan_results[[cname]]$chain
  }

  message("  [chain] Reconciling unclaimed tokens ...")
  rec <- reconcile(scan_results, chains, fused_peaks, peak_pairs, snp_peaks,
                   params = params)

  event_table <- build_event_table(rec$events)

  message(sprintf("  [chain] Done. %d events called (%d high confidence, %d review).",
    nrow(event_table),
    sum(event_table$confidence == "high"),
    sum(event_table$confidence == "review")
  ))

  list(
    chains          = chains,
    events          = rec$events,
    unclaimed_peaks = rec$unclaimed_peaks,
    unclaimed_loh   = rec$unclaimed_loh,
    event_table     = event_table
  )
}
