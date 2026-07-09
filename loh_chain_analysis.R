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
    # haplotype_label IS the edge type for singleton peaks (compound_binary
    # normalizes to "binary" here -- see normalize_edge_type_label())
    if (!"best_edge_type" %in% names(sp) && "haplotype_label" %in% names(sp))
      sp[, best_edge_type := normalize_edge_type_label(haplotype_label)]
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

    # ── Interstitial fixed tokens: flank-relative depth ratio ─────────────────
    # For every HET-bounded (non-terminal) fixed tract, store the ratio of its
    # depth to the higher of its two HET flanks. The interstitial-deletion rule
    # (Rd) uses this to separate a hemizygous deletion (depth ~halved) from a
    # copy-neutral LOH tract (depth preserved). Computed from the same
    # coverage_table the terminal ratios use; skipped (NA) when unavailable.
    for (fi in which(vapply(tokens, `[[`, character(1), "type") == "F")) {
      fr <- .interstitial_flank_depth_ratio(tokens, fi, chr_name, coverage_table)
      if (!is.na(fr)) tokens[[fi]]$meta$flank_depth_ratio <- fr
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
      hal <- normalize_edge_type_label(hal)
      if (is.na(bet) || bet == "singleton") hal else bet
    },
    haplotype_label = if ("haplotype_label" %in% names(chr_fp)) haplotype_label[1] else NA_character_,
    # Read-phasing evidence for peaks rescued from "undefined" by
    # classify_peak_haplotype()'s resolve_undefined_by_phase(): the finer
    # verdict ("crossover (reciprocal)" / "gene_conversion (NCO)" / ...) and
    # the fraction of spanning reads that switch haplotype across the island.
    # NA for peaks classified from consensus alone.
    phase_call        = if ("phase_call"        %in% names(chr_fp)) phase_call[1]        else NA_character_,
    phase_switch_frac = if ("phase_switch_frac" %in% names(chr_fp)) phase_switch_frac[1] else NA_real_,
    # Bounds of the excised FIX (LOH) island the phase rescue keyed on — the
    # actual conversion/crossover tract. Used to recompute the reported event
    # footprint (.make_event) so a rescued event covers this tract rather than
    # the broad peak window or a neighbouring LOH tract the scanner anchored to.
    phase_island_start = if ("phase_island_start" %in% names(chr_fp)) phase_island_start[1] else NA_integer_,
    phase_island_end   = if ("phase_island_end"   %in% names(chr_fp)) phase_island_end[1]   else NA_integer_,
    n_spanning        = 0L,
    jaccard           = NA_real_,
    pair_edge_type    = NA_character_,
    pair_fusion_mode  = NA_character_,
    # Position of the OTHER peak that pair_edge_type actually describes (see
    # below). Lets consumers of pair_edge_type (classify_two_binary_junction)
    # verify the borrowed verdict is really about the junction they're
    # evaluating, not a different neighboring pair.
    pair_partner_pos  = NA_real_,
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
          grp_rep$n_spanning[ri]       <- best_pair$n_spanning
          grp_rep$jaccard[ri]          <- best_pair$jaccard
          grp_rep$pair_edge_type[ri]   <- best_pair$edge_type
          grp_rep$pair_fusion_mode[ri] <- if ("fusion_mode" %in% names(best_pair))
                                            best_pair$fusion_mode else NA_character_
          # Only a defined "partner" when pos is literally one of the pair's
          # two endpoints (the normal case) — leave NA if pos merely falls
          # inside the pair's span, since there's no single peak to name.
          grp_rep$pair_partner_pos[ri] <- if (!is.na(best_pair$snp_pos_a) && best_pair$snp_pos_a == pos)
                                             best_pair$snp_pos_b
                                           else if (!is.na(best_pair$snp_pos_b) && best_pair$snp_pos_b == pos)
                                             best_pair$snp_pos_a
                                           else NA_real_
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

      # Rule 0: G H_tiny G -> G
      # A single-observation HET island (n_snps <= min_snps_for_peak) sandwiched
      # between two gap tokens is a segmentation artifact — the HMM assigned one
      # ambiguous SNP as HET in what is otherwise a fixed LOH block.  Absorbing it
      # into a plain gap lets the F-G-F rule below consolidate same-state fixed
      # segments that would otherwise be permanently split, which in turn lets
      # terminal rules (R03) see the correct telomere-reaching extent.
      if (i <= length(tokens) - 2L) {
        a <- tokens[[i]]
        h <- tokens[[i + 1L]]
        b <- tokens[[i + 2L]]

        if (a$type == "G" && h$type == "H" && b$type == "G" &&
            !.has_peak(h) &&
            !is.na(h$n_snps) && h$n_snps <= params$min_snps_for_peak &&
            (b$end - a$start) < params$merge_gap_bp) {

          merged <- make_token(
            type        = "G",
            state       = NA_character_,
            start       = a$start,
            end         = b$end,
            n_snps      = 0L,
            depth_ratio = NA_real_,
            bridged_gap = TRUE
          )
          tokens <- c(tokens[seq_len(i - 1L)],
                      list(merged),
                      tokens[seq.int(i + 3L, length(tokens))])
          changed <- TRUE
          next
        }
      }

      # Rule 0b: F [G H G] F (same state, no peak on H, span < merge_gap) -> F
      # When two same-state F tokens are separated by a G-H-G triplet whose H
      # carries no chimeric-read peak and whose total span is within merge_gap_bp,
      # the H is an artifact of the segmentation model (a handful of minority-allele
      # reads in an otherwise fixed block), not a real biological HET event.  The
      # same-state flanking context is the key prior: a genuine interruption would
      # either be larger or would carry a peak.  Absorbing all five tokens into one
      # F lets R01/R02 see the correct telomere-reaching extent without raising the
      # general min_snps_for_peak threshold that Rule 0 uses.
      if (i <= length(tokens) - 4L) {
        a  <- tokens[[i]]
        g1 <- tokens[[i + 1L]]
        h  <- tokens[[i + 2L]]
        g2 <- tokens[[i + 3L]]
        b  <- tokens[[i + 4L]]

        if (a$type == "F" && g1$type == "G" && h$type == "H" && g2$type == "G" && b$type == "F" &&
            !is.na(a$state) && !is.na(b$state) && a$state == b$state &&
            !.has_peak(h) &&
            (g2$end - g1$start) < params$merge_gap_bp) {

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
                      tokens[seq.int(i + 5L, length(tokens))])
          changed <- TRUE
          next
        }
      }

      # Rule 0c: F G F_opp G F (same outer states, opposite-state middle ≤ min_snps,
      # no peak on middle, total inner span < merge_gap) -> merged F
      # A single- or two-SNP fixed segment of opposite state sandwiched between two
      # same-state fixed regions is a HMM segmentation artifact (one ambiguous read
      # misclassifying a position), not a true LOH boundary switch.  Absorbing it
      # lets downstream rules (R02c, R03) treat the composite LOH as a single token.
      if (i <= length(tokens) - 4L) {
        a  <- tokens[[i]]
        g1 <- tokens[[i + 1L]]
        b  <- tokens[[i + 2L]]
        g2 <- tokens[[i + 3L]]
        cc <- tokens[[i + 4L]]

        if (a$type == "F" && g1$type == "G" && b$type == "F" && g2$type == "G" && cc$type == "F" &&
            !is.na(a$state) && !is.na(b$state) && !is.na(cc$state) &&
            a$state == cc$state && a$state != b$state &&
            !is.na(b$n_snps) && b$n_snps <= params$min_snps_for_peak &&
            !.has_peak(b) &&
            (g2$end - g1$start) < params$merge_gap_bp) {

          merged <- make_token(
            type        = "F",
            state       = a$state,
            start       = a$start,
            end         = cc$end,
            n_snps      = as.integer(sum(c(a$n_snps, cc$n_snps), na.rm = TRUE)),
            depth_ratio = mean(c(a$depth_ratio, cc$depth_ratio), na.rm = TRUE),
            peak_over   = a$peak_over %||% b$peak_over %||% cc$peak_over,
            peak_left   = a$peak_left,
            peak_right  = cc$peak_right,
            bridged_gap = TRUE,
            meta        = c(a$meta, cc$meta)
          )
          tokens <- c(tokens[seq_len(i - 1L)],
                      list(merged),
                      tokens[seq.int(i + 5L, length(tokens))])
          changed <- TRUE
          next
        }
      }

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
  #
  # Exception 1: automatic pair classifications (LOH-based crossovers/GCs where
  # fusion_mode = "automatic") use LOH pattern evidence, not chimeric-read
  # spanning evidence.  n_spanning for these reflects incidental spanning reads
  # over the LOH tract, not the primary evidence, so don't gate on it.
  #
  # Exception 2: gene_conversion and internal_crossover peaks are excluded from
  # compute_peak_pairs() entirely (FUSION_HEURISTICS$excluded_peak_classes).
  # Any n_spanning they carry was copied from an unrelated flanking pair whose
  # span happens to cover this peak's position — not evidence for or against
  # this peak.  The per-read haplotype pattern from classify_peak_haplotype()
  # is the authoritative call for these types, so bypass the coverage gate.
  auto_pair      <- isTRUE(peak$pair_fusion_mode == "automatic") &&
                    isTRUE(peak$pair_edge_type %in% c("crossover", "gene_conversion"))
  self_classifying <- isTRUE(edge_type %in% c("gene_conversion", "internal_crossover"))
  has_count <- !is.na(n_spanning) && n_spanning > 0L
  if (has_count && n_spanning < params$min_span && !auto_pair && !self_classifying)
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

  # Only gate on coverage when we actually have a real count.
  # Bypass if either flanking peak was auto-classified from LOH evidence.
  auto_pair <- isTRUE(left_peak$pair_fusion_mode  == "automatic") ||
               isTRUE(right_peak$pair_fusion_mode == "automatic")
  has_count <- !is.na(ns) && ns > 0L
  if (has_count && ns < params$min_span && !auto_pair)
    return(list(call = "AMBIGUOUS", reason = "low_coverage",
                n_support = ns))

  # Complementarity: left peak should be AAABBB, right peak BBBAAAA
  # We infer orientation from peak position relative to the interstitial state.
  # "Complementary" means: left peak switches INTO interstitial state,
  # right peak switches OUT of interstitial state.
  # Both switching the same direction -> NCO_GC_large_tract (copy-neutral LOH)
  # Complementary -> CO_GC

  # Use the pair's edge_type from peak_pairs if the pair spans the tract
  # (this is populated by compute_peak_pairs when loh_in_gap == TRUE) — but
  # only trust it when it actually describes the relationship between THESE
  # two peaks. pair_edge_type is attached to a peak from whichever pair best
  # matches that peak's own position; for a peak with a neighbor on both
  # sides, that pair is not necessarily the one formed with the peak on the
  # far side of THIS interstitial tract (e.g. it may describe a small
  # fragment on the peak's other flank instead). pair_partner_pos records
  # which peak that stored verdict is actually about, so we can confirm it
  # lines up with the other junction peak before trusting it.
  # (Plain %||% doesn't work here: an explicitly-NA field is not NULL, so
  # %||% never falls through to the other side.)
  left_pos  <- left_peak$fused_pos_bp  %||% left_peak$snp_pos
  right_pos <- right_peak$fused_pos_bp %||% right_peak$snp_pos

  left_valid  <- !is.null(left_peak$pair_edge_type) && !is.na(left_peak$pair_edge_type) &&
                 !is.na(left_peak$pair_partner_pos %||% NA_real_) && !is.na(right_pos) &&
                 left_peak$pair_partner_pos == right_pos
  right_valid <- !is.null(right_peak$pair_edge_type) && !is.na(right_peak$pair_edge_type) &&
                 !is.na(right_peak$pair_partner_pos %||% NA_real_) && !is.na(left_pos) &&
                 right_peak$pair_partner_pos == left_pos

  pair_et <- if (left_valid) left_peak$pair_edge_type
             else if (right_valid) right_peak$pair_edge_type
             else NA_character_

  if (!is.na(pair_et)) {
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

# Flank-relative depth ratio for an INTERSTITIAL fixed token: the token's mean
# real read depth divided by the *higher* of its two nearest non-fixed (HET)
# flank depths (skipping unscored G gaps). A hemizygous deletion of one homolog
# yields a fixed-allele tract whose total depth is ~half its flanking diploid
# regions; a copy-neutral LOH tract (gene conversion / crossover) keeps full
# depth, giving a ratio near 1. The MAX flank — not the average — is the
# reference so a genuine ~half-depth deletion is still detected when one flank
# is itself depth-depressed (e.g. a short sub-telomeric HET stretch that
# sequences shallow for reasons unrelated to copy number). Returns NA unless a
# coverage_table was supplied and both flanks are HET-resolvable.
.interstitial_flank_depth_ratio <- function(tokens, idx, chr_name, coverage_table) {
  li <- .nearest_nonfixed_left(tokens, idx)
  ri <- .nearest_nonfixed_right(tokens, idx)
  if (is.null(li) || is.null(ri)) return(NA_real_)   # not HET-bounded both sides
  d_f <- .lookup_mean_depth(coverage_table, chr_name,
                            tokens[[idx]]$start, tokens[[idx]]$end)
  d_l <- .lookup_mean_depth(coverage_table, chr_name,
                            tokens[[li]]$start, tokens[[li]]$end)
  d_r <- .lookup_mean_depth(coverage_table, chr_name,
                            tokens[[ri]]$start, tokens[[ri]]$end)
  ref <- suppressWarnings(max(d_l, d_r, na.rm = TRUE))
  if (is.na(d_f) || !is.finite(ref) || ref <= 0) return(NA_real_)
  d_f / ref
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

.edge_peak <- function(tok, side) {
  direct <- if (side == "left") tok$peak_left else tok$peak_right
  if (!is.null(direct)) return(direct)
  pk <- tok$peak_over
  if (is.null(pk)) return(NULL)
  pos <- pk$fused_pos_bp %||% pk$snp_pos
  if (is.na(pos)) return(NULL)
  if (side == "left"  && pos <= tok$start) return(pk)
  if (side == "right" && pos >= tok$end)   return(pk)
  NULL
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

  # Read-phasing evidence carried from the associated peak(s): the largest
  # per-read switch fraction across the event's peaks (NA when none was
  # phase-resolved), i.e. the fraction of island-spanning reads that switched
  # haplotype -- direct support for a CO/GC call rescued from "undefined".
  # Gather peaks from evidence_peaks AND from the involved tokens' attached
  # peaks (peak_over/left/right) so the fraction surfaces regardless of which
  # rule fired the event or whether it passed evidence_peaks explicitly.
  tok_peaks  <- unlist(lapply(tokens_involved, function(t)
    list(t$peak_over, t$peak_left, t$peak_right)), recursive = FALSE)
  all_peaks  <- c(evidence_peaks, tok_peaks)
  phase_frac <- if (length(all_peaks) > 0) {
    fr <- suppressWarnings(as.numeric(unlist(lapply(all_peaks, function(p)
      if (is.null(p)) NA_real_ else p$phase_switch_frac %||% NA_real_))))
    if (all(is.na(fr))) NA_real_ else max(fr, na.rm = TRUE)
  } else NA_real_

  # ── Re-anchor the footprint of a phase-rescued event to its excised island ─
  # A peak recovered from "undefined" carries the FIX (LOH) island bounds the
  # phasing keyed on. The motif scanner may have matched an adjacent tract
  # (e.g. the REF_fixed tract flanking a het-bounded ALT island), giving a span
  # that doesn't even overlap the real event. When any involved peak supplies
  # island bounds, report the event over that tract instead.
  isl_starts <- suppressWarnings(as.integer(unlist(lapply(all_peaks, function(p)
    if (is.null(p)) NA_integer_ else p$phase_island_start %||% NA_integer_))))
  isl_ends   <- suppressWarnings(as.integer(unlist(lapply(all_peaks, function(p)
    if (is.null(p)) NA_integer_ else p$phase_island_end   %||% NA_integer_))))
  if (any(!is.na(isl_starts)) && any(!is.na(isl_ends))) {
    ev_start <- min(isl_starts, na.rm = TRUE)
    ev_end   <- max(isl_ends,   na.rm = TRUE)
  }

  list(
    event_class       = call,
    chrom             = chrom,
    start             = as.integer(ev_start),
    end               = as.integer(ev_end),
    length_bp         = as.integer(ev_end - ev_start),
    n_support         = as.integer(n_support),
    peak_edge_types   = pk_edge,
    phase_switch_frac = phase_frac,
    notes             = notes,
    tokens            = tokens_involved   # keep for downstream use / plotting
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

# Rule 1i: Interstitial deletion — H [F] H with a real depth drop.
# A hemizygous deletion of one homolog leaves a fixed-allele tract (the retained
# homolog's SNPs, so the region reads as REF_fixed or ALT_fixed) whose total
# read depth is ~half that of the flanking diploid HET regions. That depth drop
# is what separates it from a copy-neutral LOH tract (gene conversion /
# crossover), which keeps full depth. Uses SNP-site coverage only — no CNV or
# breakpoint modelling — so it fires only on a clear drop below params$depth_drop
# vs. the higher HET flank, and is reported at "review" confidence. Placed above
# the recombination rules so a genuinely deleted tract is not mis-called as a
# gene conversion when a coincident peak sits inside it.
rule_interstitial_deletion <- list(
  id = "Rd_interstitial_deletion",
  match_fn = function(tokens, i, chain, params) {
    tok <- tokens[[i]]
    if (tok$type != "F") return(NULL)
    # HET-bounded on both sides (skipping unscored G gaps). A TEL-adjacent
    # fixed tract is a terminal deletion (R01), handled above — not here.
    if (is.null(.nearest_nonfixed_left(tokens, i)) ||
        is.null(.nearest_nonfixed_right(tokens, i))) return(NULL)
    dr <- tok$meta$flank_depth_ratio
    if (is.null(dr) || is.na(dr) || dr >= params$depth_drop) return(NULL)
    list(span = i, f_tok = tok, ratio = dr)
  },
  fire_fn = function(m, chain, params) {
    ev <- .make_event("DELETION", chain$chrom, list(m$f_tok),
                      notes = sprintf(paste0("interstitial hemizygous deletion; ",
                                             "flank_depth_ratio=%.2f < %.2f (SNP-site depth)"),
                                      m$ratio, params$depth_drop))
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
      # Use junction-peak helpers so we catch peaks that land inside the
      # adjacent non-fixed token rather than directly on F.
      pk <- if (direction == "rev")
        .left_junction_peak(ftk, htk, params)    # H [●F] TEL: junction at H→F boundary
      else
        .right_junction_peak(ftk, htk, params)   # TEL [F●] H: junction at F→H boundary
      # If htk is an unscored G gap (no peaks), look one position further to
      # the H token that carries the junction peak.  This handles the common
      # subtelomeric case: H ends at the last SNP, an unscored gap follows,
      # then the terminal F starts — the peak lives on H, not on G.
      if (is.null(pk) && htk$type == "G") {
        neighbor_i <- if (direction == "rev") h_i - 1L else h_i + 1L
        if (neighbor_i >= 1L && neighbor_i <= n && tokens[[neighbor_i]]$type == "H") {
          pk <- if (direction == "rev")
            .left_junction_peak(ftk, tokens[[neighbor_i]], params)
          else
            .right_junction_peak(ftk, tokens[[neighbor_i]], params)
        }
      }
      if (is.null(pk)) return(NULL)

      # If the boundary peak's edge_type describes an interior excursion
      # (gene_conversion/crossover) rather than a single clean terminal
      # switch, and R06 (rule_opp_sandwich) can resolve ftk as a same-state
      # run with a small embedded opposite-state fragment further down the
      # chain, defer to R06 instead of reporting AMBIGUOUS here. Without this,
      # R02 — evaluated from the TEL token, one position before ftk — claims
      # ftk before the scanner ever reaches ftk's own index, so R06 (which
      # only matches when anchored AT the F token) never gets the chance its
      # priority-list position implies it should have. This only fires for
      # the forward direction; R06 looks rightward from f_i, which does not
      # mirror the reverse (H [●F] TEL) orientation, so the reverse case is
      # unaffected and keeps its prior behavior.
      edge_type_here <- pk$best_edge_type %||% pk$edge_type
      if (direction == "fwd" &&
          !is.null(edge_type_here) && !is.na(edge_type_here) &&
          edge_type_here %in% c("gene_conversion", "crossover") &&
          !is.null(rule_opp_sandwich$match_fn(tokens, f_i, chain, params))) {
        return(NULL)
      }

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

    # A terminal LOH tract reaching the telomere IS a terminal LOH/crossover —
    # that is the defining evidence, independent of the junction peak's fine
    # structure. A clean "binary" switch is the textbook case, but a
    # gene_conversion / crossover edge type at the boundary does not negate the
    # terminal event: it almost always reflects a tiny artifactual HET (or
    # opposite-state) island sitting right at the junction, which biases the
    # per-peak read classifier away from "binary" even though the tract still
    # runs to the telomere (see RAD5_v2_04 chrXIV @363616, chrXVI @156908).
    # So any resolved switch-type junction peak → CO_TERM; only a genuinely
    # unresolvable edge type is left ambiguous. "undefined" (too little SNP
    # coverage to type the switch) stays CO_TERM_PROBABLE (review confidence).
    known_switch <- !is.null(edge_type) && !is.na(edge_type) &&
      edge_type %in% c("binary", "gene_conversion", "crossover", "internal_crossover")
    call <- if (has_count && n_spanning < params$min_span) {
      "AMBIGUOUS(low_coverage)"
    } else if (is.null(edge_type) || is.na(edge_type) || edge_type == "undefined") {
      "CO_TERM_PROBABLE"
    } else if (known_switch) {
      "CO_TERM"
    } else {
      paste0("AMBIGUOUS(terminal_non_binary_peak:", edge_type, ")")
    }

    ev <- .make_event(call, chain$chrom, list(m$f_tok),
                      evidence_peaks = list(m$peak),
                      n_support = n_spanning,
                      notes = paste0("terminal; edge=", edge_type %||% "NA"))
    list(event = ev, rewrite = NULL,
         claims = list(peak = m$peak, loh = m$f_tok))
  }
)

# Rule 2g: Terminal LOH whose internal haplotype switch falls in a SNP gap.
#
# Pattern (forward):  H [F:X ●binary] {G} [F:Y] ... {G} [F:Z]{tel}
#         (reverse): {tel}[F:Z] ... [F:Y] {G} [●binary F:X] H
#
# A terminal LOH tract can carry an internal allele switch (X→Y, e.g. ALT→REF)
# whose junction lands in a wide SNP-desert gap, so no chimeric-read peak can
# ever mark it — while the tract's PROXIMAL (H-facing) boundary does carry a
# clean binary switch peak. R02 can't fire (the telomere-reaching F is peakless
# and non-adjacent to the H); R03 can't fire (it needs a peak on BOTH of the
# small LOH's junctions). The tract still runs to the telomere, so it is a
# terminal LOH/crossover. We walk from the binary-peaked proximal F across G
# gaps and further F tokens (of ANY state — the whole point is that the switch
# is real but unobservable) to a telomere-reaching F, and call CO_TERM.
# Requires >= 2 F tokens so it never competes with R02's single-F terminal.
rule_terminal_loh_gapped <- list(
  id = "R02g_terminal_loh_gapped",
  match_fn = function(tokens, i, chain, params) {
    tok <- tokens[[i]]
    if (tok$type != "F") return(NULL)

    # Collect a run of F tokens starting at i, walking one direction through
    # G gaps only (no intervening H — that would be a real, peak-able boundary
    # this rule must not swallow). Returns the {toks, last_idx} reached.
    walk_run <- function(dir) {
      run       <- list(tok)
      last_idx  <- i
      repeat {
        nf <- if (dir == "fwd") .nearest_fixed_right(tokens, last_idx)
              else              .nearest_fixed_left(tokens, last_idx)
        if (is.null(nf)) break
        # only G tokens may separate the F run (an intervening H is a real,
        # peak-able boundary this rule must not swallow)
        lo <- min(last_idx, nf) + 1L; hi <- max(last_idx, nf) - 1L
        if (hi >= lo && any(vapply(tokens[lo:hi], function(t) t$type == "H", logical(1))))
          break
        run <- c(run, list(tokens[[nf]])); last_idx <- nf
      }
      list(run = run, last_idx = last_idx)
    }

    try_dir <- function(dir) {
      # Proximal binary junction peak on the H-facing side of tok.
      ctx_i <- if (dir == "fwd") .nearest_nonfixed_left(tokens, i)
               else              .nearest_nonfixed_right(tokens, i)
      if (is.null(ctx_i)) return(NULL)
      pk <- if (dir == "fwd") .left_junction_peak(tok, tokens[[ctx_i]], params)
            else               .right_junction_peak(tok, tokens[[ctx_i]], params)
      if (is.null(pk)) return(NULL)
      et <- pk$best_edge_type %||% pk$edge_type
      if (!isTRUE(et == "binary")) return(NULL)
      if (is.na(.terminal_depth_ratio(tok)) ||
          .terminal_depth_ratio(tok) < params$depth_drop) return(NULL)

      w <- walk_run(dir)
      if (length(w$run) < 2L) return(NULL)                       # need a gap-split continuation
      distal <- tokens[[w$last_idx]]
      if (!.is_telomeric(distal, chain, params)) return(NULL)    # must reach the telomere
      if (.has_peak(distal)) return(NULL)                        # distal switch is genuinely peakless

      span <- range(c(i, w$last_idx, ctx_i))
      list(span = span, f_toks = w$run, peak = pk, direction = dir)
    }

    m <- try_dir("fwd"); if (!is.null(m)) return(m)
    try_dir("rev")
  },
  fire_fn = function(m, chain, params) {
    n_raw      <- m$peak$n_spanning %||% NA_integer_
    n_spanning <- if (is.null(n_raw) || is.na(n_raw)) NA_integer_ else as.integer(n_raw)
    ev <- .make_event("CO_TERM", chain$chrom, m$f_toks,
                      evidence_peaks = list(m$peak),
                      n_support = n_spanning,
                      notes = paste0("terminal LOH; internal allele switch in a ",
                                     "SNP-gap (no junction peak); proximal boundary binary"))
    list(event = ev, rewrite = NULL,
         claims = list(peak = m$peak, loh = m$f_toks))
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

# Rule 2c: Terminal LOH, TEL-adjacent HET interposed —
#   TEL [H] [F●], TEL [H] [G] [F●], [●F] [H] TEL, or [●F] [G] [H] TEL.
#
# Some chromosomes carry a HET region directly adjacent to a telomere that is
# not a true biological heterozygosity but an alignment artifact: reads from
# other chromosomes (or the other parental homolog's sub-telomeric repeats)
# map to the reference here, inflating apparent heterozygosity.  When the true
# LOH (F) starts just inside this artifactual HET zone and has a binary peak on
# its distal boundary, the standard TEL [F●] H (R02) pattern cannot fire because
# the HET layer sits between TEL and F.  This rule recognises the extended
# TEL [H] [F●] pattern and calls it CO_TERM_PROBABLE (review confidence) — but
# only when F itself reaches within tel_tol_bp of the telomere, i.e. the
# interposed HET is a thin subtelomeric artifact. A wide interposed HET leaves
# F interstitial, and the event is left to R11b (GC_ONE_SIDED) instead.
# G-transparent variants (G gap between H and F) are also checked to handle
# cases where an unscored coordinate gap separates the HET from the LOH token.
rule_tel_adjacent_het_loh <- list(
  id = "R02c_tel_adjacent_het_loh",
  match_fn = function(tokens, i, chain, params) {
    n <- length(tokens)

    check_tahl <- function(tel_i, h_i, f_i, direction) {
      tel <- tokens[[tel_i]]
      htk <- tokens[[h_i]]
      ftk <- tokens[[f_i]]
      if (tel$type != "TEL" || !.is_non_fixed(htk) || ftk$type != "F") return(NULL)
      # The interposed HET must be a THIN subtelomeric misalignment artifact, so
      # the LOH (F) itself has to sit within tel_tol_bp of the chromosome end on
      # the TEL-facing side. A wide interposed HET (real diploid arm) leaves F
      # far from the telomere — that is an interstitial LOH flanked on one side
      # by a peak, not a terminal event; let a downstream interstitial rule
      # (R11b → GC_ONE_SIDED) handle it rather than over-calling CO_TERM_PROBABLE.
      f_reaches_tel <- if (direction == "fwd") {
        ref_start <- if (!is.null(chain$snp_start) && !is.na(chain$snp_start))
          chain$snp_start else 1L
        ftk$start <= ref_start + params$tel_tol_bp
      } else {
        ref_end <- if (!is.null(chain$snp_end) && !is.na(chain$snp_end))
          chain$snp_end else chain$chr_len
        ftk$end >= ref_end - params$tel_tol_bp
      }
      if (!f_reaches_tel) return(NULL)
      # Depth must be consistent with LOH, not a deletion.
      if (is.na(.terminal_depth_ratio(ftk)) ||
          .terminal_depth_ratio(ftk) < params$depth_drop) return(NULL)
      # No peak at the H-F (proximal/TEL-side) boundary — if there is one, a
      # different rule should handle it (e.g. R02 if TEL were adjacent to F).
      pk_prox <- if (direction == "fwd")
        .left_junction_peak(ftk, htk, params)    # TEL [H] [●F]: H→F boundary
      else
        .right_junction_peak(ftk, htk, params)   # [F●] [H] TEL: F→H boundary
      if (!is.null(pk_prox)) return(NULL)
      # F must have a binary peak on its distal (away-from-TEL) boundary.
      pk_distal <- if (direction == "fwd") {
        next_i <- f_i + 1L
        if (next_i > n) return(NULL)
        .right_junction_peak(ftk, tokens[[next_i]], params)
      } else {
        prev_i <- f_i - 1L
        if (prev_i < 1L) return(NULL)
        .left_junction_peak(ftk, tokens[[prev_i]], params)
      }
      if (is.null(pk_distal)) return(NULL)
      et <- pk_distal$best_edge_type %||% pk_distal$edge_type
      if (!isTRUE(et == "binary")) return(NULL)
      list(span      = c(min(tel_i, h_i, f_i), max(tel_i, h_i, f_i)),
           f_tok     = ftk,
           h_tok     = htk,
           peak      = pk_distal,
           direction = direction)
    }

    # Forward: TEL [H] [F●]
    if (i + 2L <= n) {
      m <- check_tahl(i, i + 1L, i + 2L, "fwd")
      if (!is.null(m)) return(m)
    }
    # Forward G-transparent: TEL [H] [G] [F●]
    # Handles the case where an unscored gap sits between the subtelomeric H and F.
    if (i + 3L <= n && tokens[[i + 2L]]$type == "G") {
      m <- check_tahl(i, i + 1L, i + 3L, "fwd")
      if (!is.null(m)) return(m)
    }
    # Reverse: [●F] [H] TEL
    if (i + 2L <= n) {
      m <- check_tahl(i + 2L, i + 1L, i, "rev")
      if (!is.null(m)) return(m)
    }
    # Reverse G-transparent: [●F] [G] [H] TEL
    # Handles the case where an unscored gap sits between F and the subtelomeric H.
    if (i + 3L <= n && tokens[[i + 1L]]$type == "G") {
      m <- check_tahl(i + 3L, i + 2L, i, "rev")
      if (!is.null(m)) return(m)
    }
    NULL
  },
  fire_fn = function(m, chain, params) {
    edge_type  <- m$peak$best_edge_type %||% m$peak$edge_type
    n_raw      <- m$peak$n_spanning %||% NA_integer_
    n_spanning <- if (is.null(n_raw) || is.na(n_raw)) NA_integer_ else as.integer(n_raw)
    has_count  <- !is.na(n_spanning) && n_spanning > 0L

    call <- if (has_count && n_spanning < params$min_span)
      "AMBIGUOUS(low_coverage)"
    else
      "CO_TERM_PROBABLE"

    ev <- .make_event(call, chain$chrom, list(m$f_tok),
                      evidence_peaks = list(m$peak),
                      n_support = n_spanning,
                      notes = paste0("tel_adjacent_het; edge=", edge_type %||% "NA",
                                     "; LOH->HET->TEL without junction peak; ",
                                     "probable subtelomeric misalignment"))
    list(event = ev, rewrite = NULL,
         claims = list(peak = m$peak, loh = m$f_tok))
  }
)

# Rule 3: TCO-captured-TCO — H ●[F]● F̄→TEL (opposite fixed reaches telomere,
#         peaks on both sides)
rule_tco_captured_tco <- list(
  id = "R03_tco_captured_tco",
  match_fn = function(tokens, i, chain, params) {
    # Pattern: H [F] F̄{tel}  (with peaks on both junctions of F).
    # G-transparent like R06/R10/R11: real loh_segments data rarely tiles a
    # chromosome with zero gaps, so ftk/fbar are usually separated from H and
    # each other by small unscored G tokens rather than sitting at literal
    # i, i+1, i+2 offsets.
    check <- function(h_idx, f_idx, fb_idx) {
      htk <- tokens[[h_idx]]; ftk <- tokens[[f_idx]]; fbar <- tokens[[fb_idx]]
      if (!.is_non_fixed(htk)) return(NULL)
      if (is.na(ftk$state) || is.na(fbar$state)) return(NULL)
      if (fbar$state != .opp_state(ftk$state)) return(NULL)

      # fbar_pk_ctx: the token whose left boundary carries the junction peak.
      # This is always the FIRST fbar token (the one immediately adjacent to
      # the small-LOH ftk).  When the terminal LOH is segmented into multiple
      # same-state F tokens by a tiny HET island (e.g. a single ambiguous SNP
      # in the full-genome HMM), the first segment still holds peak_left while
      # the telomere-reaching continuation has no peak.
      fbar_pk_ctx  <- fbar
      fbar_extra   <- list()   # intermediate same-state F tokens (if any)
      fb_tel_idx   <- fb_idx

      if (!.is_telomeric(fbar, chain, params)) {
        cont <- .find_telomeric_continuation(tokens, fb_idx, fbar$state, chain, params)
        if (is.null(cont)) return(NULL)
        # intermediate tokens between fbar and the telomeric end (usually empty)
        fbar_extra <- lapply(cont[-length(cont)], `[[`, "tok")
        fb_tel_idx <- cont[[length(cont)]]$idx
        fbar       <- cont[[length(cont)]]$tok   # telomere-reaching token
      }

      # .left_junction_peak/.right_junction_peak fall back to a neighbour's
      # peak_over (position-checked), not just peak_left/peak_right — real
      # junction peaks land there far more often than directly on the token.
      pk_l <- .left_junction_peak(ftk, htk, params)
      pk_r <- .right_junction_peak(ftk, fbar_pk_ctx, params)  # peak on original boundary
      if (is.null(pk_l) || is.null(pk_r)) return(NULL)
      list(f_tok = ftk, fbar_tok = fbar, fbar_pk_ctx = fbar_pk_ctx,
           fbar_extra = fbar_extra, fb_tel_idx = fb_tel_idx,
           h_tok = htk, pk_l = pk_l, pk_r = pk_r)
    }

    if (!.is_non_fixed(tokens[[i]])) return(NULL)

    # Early: anchor on htk itself, walk outward toward the telomere. Fires
    # on the scanner's first pass through htk's position — only works if
    # fbar is already at its final extent by then.
    fi  <- .nearest_fixed_right(tokens, i)
    fbi <- if (!is.null(fi)) .nearest_fixed_right(tokens, fi) else NULL
    if (!is.null(fi) && !is.null(fbi)) {
      m <- check(i, fi, fbi)
      if (!is.null(m)) return(c(m, list(span = c(i, m$fb_tel_idx))))
    }

    # Late: anchor on the non-fixed gap between ftk and fbar. fbar is very
    # often the product of one or more R06 (opp_sandwich) merges that
    # haven't happened yet on the scanner's first pass through htk's
    # position — but scan_chain() backs up to just before each merge's
    # span after every rewrite, so this gap position gets re-checked every
    # time, with fbar's fully-merged, telomere-reaching extent by the last
    # one. Tried second so the early form (cheaper, no merge required)
    # still wins outright when it's already sufficient.
    fi <- .nearest_fixed_left(tokens, i)
    if (!is.null(fi)) {
      hi  <- .nearest_nonfixed_left(tokens, fi)
      fbi <- .nearest_fixed_right(tokens, i)
      if (!is.null(hi) && !is.null(fbi)) {
        m <- check(hi, fi, fbi)
        if (!is.null(m)) return(c(m, list(span = c(hi, m$fb_tel_idx))))
      }
    }

    NULL
  },
  fire_fn = function(m, chain, params) {
    # Build the full list of LOH tokens covering f_tok through the telomeric end,
    # including any intermediate same-state F tokens produced by HET-island splits.
    has_cont  <- !identical(m$fbar_pk_ctx$start, m$fbar_tok$start)
    fbar_toks <- if (has_cont) {
      c(list(m$fbar_pk_ctx), m$fbar_extra, list(m$fbar_tok))
    } else {
      list(m$fbar_tok)
    }
    all_loh <- c(list(m$f_tok), fbar_toks)
    ev <- .make_event("TCO_CAPTURED_TCO", chain$chrom,
                      all_loh,
                      evidence_peaks = list(m$pk_l, m$pk_r),
                      n_support = as.integer(.ns(m$pk_l) + .ns(m$pk_r)),
                      notes = "terminal CO over earlier terminal CO")
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk_l, m$pk_r),
                       loh  = all_loh))
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

# Rule 6: F ●[F̄]● F [●F̄]● F ... — a chain of small interstitial fragments
#          alternating against a common background state
#          (GC from a prior event captured by terminal CO)
#
#   A single embedded excursion (L [Fbar] R, L/R sharing one state, Fbar the
#   other) is the base case. But consecutive excursions chain: L [F1] R1
#   [F2] R2 [F3] R3 ... as long as each Fk is small relative to what
#   immediately follows it (Rk), each Rk in turn becomes the new "L" for
#   testing whether the NEXT gap also holds a small fragment, rather than
#   being silently absorbed as part of one giant merge. Each Fk gets its own
#   reported event — e.g. R06 walking L=ALT [F1=REF] R1=ALT [F2=REF] R2=ALT
#   reports F1 and F2 as two separate excursions, using R1 purely as the
#   shared background connecting them (R1 itself is not reported: it plays
#   the same role L does for F1, and the same role L plays for F2 relative to
#   R1 — every interior token alternates between "fragment" and "confirming
#   background" roles as the walk proceeds, and only fragment-role tokens are
#   reported). Only the walk's outermost L and final R get merged into one
#   background token; every fragment token in between is consumed and
#   reported without being folded into that merge.
rule_opp_sandwich <- list(
  id = "R06_opp_sandwich",
  match_fn = function(tokens, i, chain, params) {
    l <- tokens[[i]]
    if (l$type != "F") return(NULL)

    fragments  <- list()
    cur_bg     <- l
    cur_bg_idx <- i

    repeat {
      # G-transparent: real loh_segments data almost never tiles a chromosome
      # with zero gaps, so the fragment and its flanks are usually separated
      # by small unscored G tokens rather than being strictly adjacent list
      # entries. Skip over them to find the actual fixed neighbours.
      fi <- .nearest_fixed_right(tokens, cur_bg_idx)
      if (is.null(fi)) break
      ri <- .nearest_fixed_right(tokens, fi)
      if (is.null(ri)) break

      ftk <- tokens[[fi]]; r <- tokens[[ri]]
      if (is.na(cur_bg$state) || is.na(ftk$state) || is.na(r$state)) break
      if (ftk$state == cur_bg$state || cur_bg$state != r$state) break
      if (ftk$state != .opp_state(cur_bg$state)) break
      # Check "small" criterion
      combined_flank <- (cur_bg$length_bp %||% 0L) + (r$length_bp %||% 0L)
      if ((ftk$length_bp %||% Inf) >= params$small_frac * combined_flank) break

      # cur_bg can be a large, already-resolved accumulation from earlier
      # steps of this same walk — its size reflects history, not whether ftk
      # is truly a small excursion here and now. r, by contrast, is always
      # the next not-yet-merged raw segment, so it reflects the actual local
      # scale at this point in the scan. When ftk is not smaller than r, ftk
      # is not a small anomaly embedded in cur_bg's background — it is at
      # least as large as r itself, meaning ftk is really the leading edge of
      # a *different*, potentially much larger region that happens to look
      # "small" only because it is being measured against cur_bg's inflated
      # combined_flank. Stopping the walk here leaves that boundary for
      # R02/R03 to close out instead of silently misreporting the true
      # transition point. Comparing against r alone (not cur_bg, not
      # cur_bg+r) is what lets a later scan position reach ftk directly on
      # its own turn, where this same rule can correctly grow ftk's own
      # state against the (now fresh) r beyond it.
      if ((ftk$length_bp %||% Inf) >= (r$length_bp %||% 0L)) break

      # Peak evidence must be ftk's OWN attachment, not borrowed from its
      # neighbours. Those can be large (or, after an earlier step of this
      # same walk, already a cumulative merge of several prior segments) —
      # falling back to their peak_over/peak_left/peak_right grabs a peak
      # that really belongs to a far boundary (an unrelated, earlier event),
      # not the junction adjacent to ftk. A wide fused peak window can
      # legitimately overlap (and so get attached as peak_over to) both ftk
      # and its immediate neighbour; trusting only ftk's own fields means
      # that peak is claimed by ftk's event, not stolen by the flank's merge
      # — which would otherwise starve ftk of the evidence it needs to ever
      # fire.
      pk_l <- ftk$peak_left
      pk_r <- ftk$peak_right
      pk   <- .best_peak(ftk)

      # Fragments with too few SNPs to ever generate a chimeric-read peak
      # (the same floor reconcile() already uses for POSSIBLE_GC) are still
      # eligible to be absorbed — with no peak, there's nothing to mis-claim,
      # and leaving a single-SNP noise blip unmerged would otherwise block
      # the cascade from ever reaching the far telomere.
      too_small_for_peak <- !is.na(ftk$n_snps) && ftk$n_snps < params$min_snps_for_peak

      # A peak attached as peak_over may have its SNP position (fused_pos_bp)
      # outside the fragment's own [start, end] span — this happens when a
      # wide chimeric-peak window overlaps the fragment edge but the
      # underlying SNP sits in the adjacent flank. For binary peaks this is a
      # flank-boundary marker that belongs to R03, not interior evidence for
      # this fragment. Claiming it here would block R03 from later using it
      # as the merged-flank's right-junction peak. Treat the fragment as
      # having no internal peak evidence so it is absorbed without consuming
      # the boundary peak.
      # Exception: self-classifying peaks (gene_conversion, crossover,
      # internal_crossover) whose SNP sits just past ftk$end ARE the right-
      # junction marker for this gene-conversion motif. If left unclaimed
      # they propagate to the merged-flank's peak_right and get consumed by
      # R10, producing a spuriously large NCO_GC that spans the entire merged
      # region. Claim them here so R06 fires the correct NCO_GC_in_terminal
      # instead.
      if (!is.null(pk) && !is.na(pk$fused_pos_bp) &&
          (pk$fused_pos_bp < ftk$start || pk$fused_pos_bp > ftk$end)) {
        et_pk <- pk$best_edge_type %||% pk$edge_type
        sc_pk <- isTRUE(et_pk %in% c("gene_conversion", "crossover", "internal_crossover"))
        if (!sc_pk) {
          pk             <- NULL
          too_small_for_peak <- TRUE
        }
      }

      if (is.null(pk) && !too_small_for_peak) break

      fragments[[length(fragments) + 1L]] <- list(
        f_tok = ftk, pk = pk, pk_l = pk_l, pk_r = pk_r,
        too_small_for_peak = too_small_for_peak
      )

      # r becomes the new confirming background: test whether the NEXT gap
      # (past r) also holds a small fragment instead of assuming everything
      # from here on is uniform background to be swallowed whole.
      cur_bg     <- r
      cur_bg_idx <- ri
    }

    if (length(fragments) == 0) return(NULL)

    list(span = c(i, cur_bg_idx), l_tok = l, r_tok = cur_bg, fragments = fragments)
  },
  fire_fn = function(m, chain, params) {

    events <- lapply(m$fragments, function(frag) {
      centered <- .is_roughly_centered(frag$f_tok, params$center_tol) ||
                  .is_roughly_centered(
                    make_token("F", start = frag$f_tok$start, end = frag$f_tok$end,
                               peak_over = frag$pk_l), params$center_tol)

      tract <- if (frag$too_small_for_peak && is.null(frag$pk))
        list(call = "POSSIBLE_GC", reason = "too_small_for_peak", n_support = 0L)
      else if (!is.null(frag$pk_l) && !is.null(frag$pk_r))
        classify_two_binary_junction(frag$pk_l, frag$pk_r, frag$f_tok$state, params)
      else
        classify_tract(frag$pk, params)

      call <- switch(tract$call,
        NCO_GC       = "NCO_GC_in_terminal",
        NCO_GC_LARGE = "NCO_GC_in_terminal",
        CO_GC        = "CO_GC",
        POSSIBLE_GC  = "POSSIBLE_GC",
        paste0("AMBIGUOUS(", tract$reason, ")")
      )

      # Report only this fragment's own span, not the merged flank's — the
      # merged flank is a structural rewrite for later rules (R02/R02b etc)
      # to evaluate as one unit, not part of this event. Including the
      # cumulative l_tok/r_tok span here would make each fragment's reported
      # span balloon to cover all the unrelated territory peeled before it,
      # and would duplicate whatever event later claims the fully-merged
      # flank.
      .make_event(call, chain$chrom,
                  list(frag$f_tok),
                  evidence_peaks = Filter(Negate(is.null),
                                          list(frag$pk, frag$pk_l, frag$pk_r)),
                  n_support = tract$n_support,
                  notes = paste0("small_opp_fragment; centered=", centered))
    })

    # Rewrite: combine the walk's outermost flanking F tokens into one (every
    # fragment token in between has now been consumed and reported above).
    # l_tok/r_tok's own outer-edge peak can be attached as peak_over rather
    # than peak_left/peak_right (the G gap to the next segment is often wider
    # than peak_pad_bp) — .edge_peak() recovers it by position so a terminal
    # rule downstream (R02/R02b/R03) doesn't see a peakless token where a
    # real junction peak exists.
    merged_flank <- make_token(
      type        = "F",
      state       = m$l_tok$state,
      start       = m$l_tok$start,
      end         = m$r_tok$end,
      n_snps      = as.integer(sum(c(m$l_tok$n_snps, m$r_tok$n_snps), na.rm = TRUE)),
      depth_ratio = mean(c(m$l_tok$depth_ratio, m$r_tok$depth_ratio), na.rm = TRUE),
      peak_left   = .edge_peak(m$l_tok, "left"),
      # Fall back to l_tok's own peak_right when r_tok has no outer-edge peak.
      # Without this, each R06 cascade step that merges two same-state flanks
      # loses the accumulated right-boundary peak from all prior merges,
      # leaving the final merged token peakless on its right side and starving
      # R03 of the pk_r it needs to fire.
      peak_right  = .edge_peak(m$r_tok, "right") %||% m$l_tok$peak_right,
      bridged_gap = FALSE,
      meta        = c(m$l_tok$meta, m$r_tok$meta)
    )

    all_peaks <- unlist(lapply(m$fragments, function(frag)
      Filter(Negate(is.null), list(frag$pk, frag$pk_l, frag$pk_r))), recursive = FALSE)
    all_loh   <- lapply(m$fragments, `[[`, "f_tok")

    list(events  = events,
         rewrite = list(span = m$span, replacement = list(merged_flank)),
         claims  = list(peak = all_peaks, loh = all_loh))
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

# Same skip-over-G traversal as .nearest_nonfixed_left/right, but looking for
# the next FIXED (F) token instead — used by R06 (opp_sandwich) so a small
# embedded fragment separated from its flanks by an unscored gap (the normal
# case in real data; loh_segments rarely tiles a chromosome with zero gaps)
# is still recognised as adjacent. Stops (returns NULL) if H or TEL is
# reached before any F is found.
.nearest_fixed_left <- function(tokens, i) {
  j <- i - 1L
  while (j >= 1L && tokens[[j]]$type == "G") j <- j - 1L
  if (j >= 1L && tokens[[j]]$type == "F") j else NULL
}

.nearest_fixed_right <- function(tokens, i) {
  n <- length(tokens)
  j <- i + 1L
  while (j <= n && tokens[[j]]$type == "G") j <- j + 1L
  if (j <= n && tokens[[j]]$type == "F") j else NULL
}

# Walk right from tokens[[start_i+1]], skipping G tokens and tiny H tokens
# (n_snps <= min_snps_for_peak), collecting consecutive same-state F tokens.
# Returns a list of {tok, idx} pairs ending with the last F of `state` that
# is telomeric, or NULL if no telomeric same-state F is reachable.
# Used by R03 to bridge a single-SNP HET island that splits what is logically
# one terminal LOH block into two same-state F tokens.
.find_telomeric_continuation <- function(tokens, start_i, state, chain, params) {
  n         <- length(tokens)
  j         <- start_i + 1L
  collected <- list()
  while (j <= n) {
    tk <- tokens[[j]]
    if (tk$type == "G") { j <- j + 1L; next }
    if (tk$type == "H" && !is.na(tk$n_snps) &&
        tk$n_snps <= params$min_snps_for_peak) { j <- j + 1L; next }
    if (tk$type == "F" && identical(tk$state, state)) {
      collected[[length(collected) + 1L]] <- list(tok = tk, idx = j)
      j <- j + 1L; next
    }
    break
  }
  if (length(collected) == 0) return(NULL)
  if (!.is_telomeric(collected[[length(collected)]]$tok, chain, params)) return(NULL)
  collected
}

# Left junction peak for an F token: the peak that marks the H→F transition.
# Priority: F$peak_left > left_ctx$peak_right > left_ctx$peak_over (left of F)
#           > F$peak_over (if its snp_pos is left of F.start)
#
# peak_over is attached to a token whenever the peak's window overlaps ANY
# part of that token's span (see .attach_peaks) — for a short token that's
# effectively "at the boundary", but for a long H/G run it can be a peak
# sitting deep inside, near the token's OTHER (far) end, describing a
# completely different junction. Bound the *_over fallback to merge_gap_bp:
# real junction peaks sit just past a short unscored gap (a few hundred to a
# few thousand bp); anything farther is evidence for a different boundary,
# not this one.
.left_junction_peak <- function(f_tok, left_ctx, params) {
  tol <- params$merge_gap_bp %||% 5000L
  if (!is.null(f_tok$peak_left)) return(f_tok$peak_left)
  if (!is.null(left_ctx$peak_right)) return(left_ctx$peak_right)
  pk <- left_ctx$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    # <=, not <: a junction peak commonly lands exactly at the boundary
    # (that's where the allele switch is), not strictly inside the gap.
    if (!is.na(pos) && pos <= f_tok$start && f_tok$start - pos <= tol) return(pk)
  }
  pk <- f_tok$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    if (!is.na(pos) && pos <= f_tok$start && f_tok$start - pos <= tol) return(pk)
  }
  NULL
}

# Right junction peak for an F token: the peak that marks the F→H transition.
# Priority: F$peak_right > F$peak_over (right of F.end) > right_ctx$peak_left
#           > right_ctx$peak_over (right of F.end)
# See .left_junction_peak for why the *_over fallback is distance-bounded.
.right_junction_peak <- function(f_tok, right_ctx, params) {
  tol <- params$merge_gap_bp %||% 5000L
  if (!is.null(f_tok$peak_right)) return(f_tok$peak_right)
  pk <- f_tok$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    # >=, not >: see .left_junction_peak — boundary-exact peaks are normal.
    if (!is.na(pos) && pos >= f_tok$end && pos - f_tok$end <= tol) return(pk)
  }
  if (!is.null(right_ctx$peak_left)) return(right_ctx$peak_left)
  pk <- right_ctx$peak_over
  if (!is.null(pk)) {
    pos <- pk$fused_pos_bp %||% pk$snp_pos
    if (!is.na(pos) && pos >= f_tok$end && pos - f_tok$end <= tol) return(pk)
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

    pk_l <- .left_junction_peak(tok, left_ctx, params)
    pk_r <- .right_junction_peak(tok, right_ctx, params)
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
    call  <- if (identical(tract$call, "AMBIGUOUS") && identical(tract$reason, "binary_no_pair")) {
      # Both flanking peaks are independently self-classifying ("binary"),
      # confirming a real LOH tract bridges this gap — it's just wider than
      # any read can span, so no read evidence exists to call NCO vs CO.
      # Report the confirmed gene-conversion-type event instead of leaving it
      # fully ambiguous (or letting the peaks fall through separately to
      # reconcile()'s peak-only "_subres" path).
      "GC_UNRESOLVED"
    } else {
      switch(tract$call,
        NCO_GC       = "NCO_GC",
        NCO_GC_LARGE = "NCO_GC",
        CO_GC        = "CO_GC",
        paste0("AMBIGUOUS(", tract$reason, ")")
      )
    }
    ev <- .make_event(call, chain$chrom,
                      list(m$f_tok),
                      evidence_peaks = list(m$pk_l, m$pk_r),
                      n_support = tract$n_support,
                      notes = paste0("two_binary_flanking; state=", m$f_tok$state))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk_l, m$pk_r), loh = list(m$f_tok)))
  }
)

# ── Rule R11b: One-sided binary peak (H-[F]-H, only one junction confirmed) ───
# Fires when an interstitial F is bracketed by a confirmed binary peak on
# exactly ONE side, and the other side has NO peak at all (not merely a
# non-binary one — that's a different, unresolved situation left to review).
# This covers a real biological case: one breakpoint produced a chimeric read
# peak, but the other breakpoint falls where the reference makes SNP calling
# unreliable (repeat/low-complexity sequence, coverage dropout, etc.), so no
# read can ever be scored there. R11 requires both flanks and won't fire;
# without this rule the token would fall through to reconcile() as a plain
# UNCATEGORIZED_LOH, silently discarding the one real peak's evidence.
#
# Deliberately reads pk$haplotype_label (the raw, pre-normalization type)
# rather than best_edge_type here: a compound_binary peak marks a patchy,
# repair-flipped switch whose immediate host fragment is often just that
# noise, not a genuine standalone island. Since scan_chain() claims peaks by
# position as soon as any rule uses them, letting compound_binary qualify
# here would let this rule grab the peak for the noisy fragment before the
# scan reaches the real, farther boundary it actually marks — starving
# .left_junction_peak()'s bounded peak_over fallback (and thus R03) of the
# peak it needs to chain the real tract. best_edge_type (which normalizes
# compound_binary to "binary") is still what downstream fusion/pairing and
# event reporting see — only this rule's own eligibility gate is stricter.
rule_one_sided_binary <- list(
  id = "R11b_one_sided_binary",
  match_fn = function(tokens, i, chain, params) {
    tok <- tokens[[i]]
    if (tok$type != "F") return(NULL)

    li <- .nearest_nonfixed_left(tokens, i)
    ri <- .nearest_nonfixed_right(tokens, i)
    if (is.null(li) || is.null(ri)) return(NULL)

    left_ctx  <- tokens[[li]]
    right_ctx <- tokens[[ri]]

    pk_l <- .left_junction_peak(tok, left_ctx, params)
    pk_r <- .right_junction_peak(tok, right_ctx, params)

    et_l <- pk_l$haplotype_label %||% pk_l$best_edge_type %||% pk_l$edge_type
    et_r <- pk_r$haplotype_label %||% pk_r$best_edge_type %||% pk_r$edge_type
    l_binary <- !is.null(pk_l) && isTRUE(et_l == "binary")
    r_binary <- !is.null(pk_r) && isTRUE(et_r == "binary")

    if (l_binary && is.null(pk_r)) {
      list(span = c(i, i), f_tok = tok, pk = pk_l, side = "left")
    } else if (r_binary && is.null(pk_l)) {
      list(span = c(i, i), f_tok = tok, pk = pk_r, side = "right")
    } else {
      NULL
    }
  },
  fire_fn = function(m, chain, params) {
    ev <- .make_event("GC_ONE_SIDED", chain$chrom, list(m$f_tok),
                      evidence_peaks = list(m$pk),
                      n_support = .ns(m$pk),
                      notes = paste0("one_sided_binary; confirmed_side=", m$side,
                                     "; no peak found on the other side"))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk), loh = list(m$f_tok)))
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

    pk_l <- .left_junction_peak(tok, left_ctx, params)
    pk_r <- .right_junction_peak(tok, right_ctx, params)
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
                      list(m$f_tok),
                      evidence_peaks = list(m$pk_l, m$pk_r),
                      n_support = ns,
                      notes = paste0("loh_crossover; state=", m$f_tok$state))
    list(event = ev, rewrite = NULL,
         claims = list(peak = list(m$pk_l, m$pk_r), loh = list(m$f_tok)))
  }
)

# The ordered rule set — priority highest to lowest.
# R06 goes first: it matches three consecutive F tokens (small opposite-state
# fragment sandwiched between two same-state flanks) and rewrites the chain by
# merging the flanks into one token. Its match pattern requires tokens[[i]] to
# be type F, which never overlaps the TEL-anchored patterns R01/R02/R02b need
# at that position — but it DOES overlap R10/R12, which can also match a long
# flank directly via a junction peak shared with the small fragment it
# borders. Without R06 running first, R10 claims that junction peak for the
# long flank itself (mislabeling it as a small NCO_GC), which both produces a
# wrong call AND starves R06 of the peak it needs to ever fire. Putting R06
# first lets it claim and peel the small fragment, merging the flanks; the
# scanner backs up and re-canonicalises after every rewrite (see
# scan_chain()), so a chain of several such fragments (e.g. a long terminal
# CO interrupted by multiple short embedded GCs) gets peeled one at a time
# until what remains touching the telomere is a single unified token, letting
# R02/R02b correctly fire on the full span.
# R10, R11, R12 are the interstitial rules. R03 is re-enabled below; old
# R04/R07-R09 remain disabled pending further redesign.
# R03 goes before R02/R02b: its match pattern (H [F] F-bar->TEL) shares its
# right-hand junction peak with R02/R02b's single-F terminal pattern — both
# can fire on the same boundary once an earlier R06 merge has extended the
# outer F all the way to the telomere. R03 is the more complete call (it
# also accounts for the inner, now-captured F instead of leaving it an
# orphaned UNCATEGORIZED_LOH), so it must claim that peak first.
MOTIF_RULES <- list(
  rule_opp_sandwich,             # R06 — small opposite-state fragment nested in a same-state run
  rule_terminal_deletion,        # R01
  rule_interstitial_deletion,    # Rd — HET-bounded fixed tract at ~half depth (hemizygous deletion)
  rule_tco_captured_tco,         # R03 — terminal CO over an earlier terminal CO
  rule_terminal_loh,             # R02
  rule_terminal_loh_gapped,      # R02g — terminal LOH whose internal switch falls in a SNP gap
  rule_tel_adjacent_het_loh,     # R02c — CO_TERM_PROBABLE: TEL [H] [F●], HET interposed by misalignment
  # rule_terminal_no_peak disabled: without a chimeric peak at the TEL boundary
  # there is no positive evidence for a crossover mechanism.  Leave unclaimed →
  # UNCATEGORIZED_LOH → no symbol on the overview map.
  rule_loh_crossover,            # R12 — crossover through large interstitial LOH (before R10)
  rule_peak_direct,              # R10 — gene_conversion / crossover / internal_crossover
  rule_two_binary_flanking,      # R11 — two binary peaks flanking H-[F]-H
  rule_one_sided_binary          # R11b — one confirmed binary peak, other side genuinely peak-less
  # rule_het_bounded,            # R04 — disabled (replaced by R10/R11)
  # rule_double_gc,              # R07 — disabled
  # rule_crossover_no_tract,     # R08 — disabled
  # rule_subres_tract            # R09 — disabled
)

# =============================================================================
#  STEP 4 — SCAN A SINGLE CHAIN  (candidate generation + global resolver)
# =============================================================================
#
#  Design (see the module header / MOTIF_RULES priority notes): instead of a
#  greedy left-to-right pass where the FIRST rule to reach a shared peak (by
#  scan position) wins it, the scan is split into three phases:
#
#    1. PROPOSE  — every rule is run at every anchor position, read-only. Each
#                  match becomes a *candidate* carrying the event it would fire,
#                  the evidence it would consume (peak positions + F-token
#                  identities), whether it rewrites the chain, and a score.
#    2. RESOLVE  — a single global selector picks the highest-scoring set of
#                  candidates that do not contend for the same evidence. Rule
#                  priority (MOTIF_RULES order) becomes the score, so conflicts
#                  are settled by an explicit, inspectable objective rather than
#                  by which anchor the cursor happened to reach first.
#    3. APPLY    — chosen events are committed. Structural rewrites (R06 flank
#                  merges) are applied one at a time, re-canonicalising and
#                  regenerating candidates, so "events within events" resolve by
#                  peeling the chain to a fixpoint before the terminal /
#                  interstitial rules are committed on the simplified chain.
#
#  Two candidates CONFLICT when they would claim the same chimeric-read peak
#  (peaks are shared observations that back a single event) OR the same fixed
#  (F) LOH token (each LOH tract gets one call). This reproduces the exclusivity
#  the old advancing cursor gave for free — e.g. R01 and R10 can no longer both
#  claim one terminal F — without making it depend on scan order.

# Evidence keys a fired candidate would consume: peak positions + F-token ids.
# Used by the resolver to detect conflicts and by the caller to claim peaks.
.fired_claim_keys <- function(fired, peak_pos) {
  keys <- character(0)
  if (length(peak_pos))
    keys <- c(keys, paste0("pk:", format(peak_pos, scientific = FALSE, trim = TRUE)))
  loh <- fired$claims$loh
  if (!is.null(loh)) {
    # claims$loh is either a single token (has $type) or a list of tokens
    toks <- if (!is.null(loh$type)) list(loh) else loh
    for (t in toks) {
      if (!is.null(t) && !is.null(t$start) && !is.na(t$start))
        keys <- c(keys, paste0("loh:", t$start, "-", t$end))
    }
  }
  keys
}

# Priority score: earlier rules in MOTIF_RULES score higher, so the resolver
# prefers the same rule the old priority-ordered pass preferred when two
# candidates contend for the same evidence.
.candidate_score <- function(rule_idx, n_rules) {
  (n_rules - rule_idx + 1L) * 1000L
}

# PROPOSE: run every rule at every position over the current token list,
# read-only. Skips matches whose peaks are already claimed. Fires each match
# once (fire_fn is pure) to learn the event, rewrite spec, and claim keys.
.generate_candidates <- function(tokens, chain, params, rules, claimed_peak_pos) {
  cands <- list()
  for (ri in seq_along(rules)) {
    rule <- rules[[ri]]
    for (i in seq_along(tokens)) {
      m <- rule$match_fn(tokens, i, chain, params)
      if (is.null(m)) next
      pk_pos <- .peak_pos_from_match(m)
      if (length(pk_pos) && any(pk_pos %in% claimed_peak_pos)) next
      fired <- rule$fire_fn(m, chain, params)
      cands[[length(cands) + 1L]] <- list(
        rule_idx = ri, match = m, fired = fired,
        peak_pos = pk_pos,
        keys     = .fired_claim_keys(fired, pk_pos),
        is_rewrite = !is.null(fired$rewrite),
        score    = .candidate_score(ri, length(rules))
      )
    }
  }
  cands
}

# RESOLVE: greedy weighted selection — accept candidates by descending score,
# skipping any whose evidence keys overlap an already-accepted candidate.
# Ties (equal score, i.e. same rule) break by generation order, which is
# position order — matching the old "lower position first" behaviour.
.resolve_candidates <- function(cands) {
  if (!length(cands)) return(list())
  ord <- order(vapply(cands, `[[`, numeric(1), "score"), decreasing = TRUE)
  accepted <- list()
  used     <- character(0)
  for (k in ord) {
    ck <- cands[[k]]$keys
    if (length(ck) && any(ck %in% used)) next
    accepted[[length(accepted) + 1L]] <- cands[[k]]
    used <- c(used, ck)
  }
  accepted
}

scan_chain <- function(chain, params, rules = MOTIF_RULES) {
  tokens <- chain$tokens
  chain$tokens <- tokens
  events <- list()
  claimed_peak_pos <- numeric(0)

  commit <- function(cand) {
    fired <- cand$fired
    new_events <- if (!is.null(fired$events)) fired$events else list(fired$event)
    events <<- c(events, new_events)
    claimed_peak_pos <<- unique(c(claimed_peak_pos, cand$peak_pos))
  }

  repeat {
    cands  <- .generate_candidates(tokens, chain, params, rules, claimed_peak_pos)
    if (!length(cands)) break
    chosen <- .resolve_candidates(cands)
    if (!length(chosen)) break

    # Phase 1 — structural simplification. Apply the highest-scoring rewrite
    # (chosen is already score-ordered), then re-canonicalise and regenerate so
    # the terminal/interstitial rules see the peeled, unified chain. Only the
    # rewrite's own events are committed here; other chosen candidates are
    # deferred to the next round on the simplified chain.
    rw <- Find(function(c) c$is_rewrite, chosen)
    if (!is.null(rw)) {
      commit(rw)
      span <- rw$fired$rewrite$span
      repl <- rw$fired$rewrite$replacement
      tokens <- c(tokens[seq_len(span[1] - 1L)],
                  repl,
                  if (span[2] < length(tokens))
                    tokens[seq.int(span[2] + 1L, length(tokens))]
                  else list())
      chain$tokens <- tokens
      chain        <- canonicalise(chain, params)
      tokens       <- chain$tokens
      next
    }

    # Phase 2 — no rewrites remain: commit every chosen (non-conflicting) event.
    for (c in chosen) commit(c)
    break
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
  # R06 (rule_opp_sandwich) can carry several fragments' peaks in m$fragments
  # instead of the single pk/pk_l/pk_r fields other rules use.
  if (!is.null(m$fragments)) {
    frag_pks <- unlist(lapply(m$fragments, function(frag)
      Filter(Negate(is.null), list(frag$pk, frag$pk_l, frag$pk_r))),
      recursive = FALSE)
    pks <- c(pks, frag_pks)
  }
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
      hal <- normalize_edge_type_label(hal)
      ifelse(is.na(bet) | bet == "singleton", hal, bet)
    }]
    fp[, ns_col := if ("n_read_support" %in% names(fp)) n_read_support else NA_integer_]
    fp[, .(chrom, pos_col, et_col, ns_col)]
  } else if (!is.null(snp_peaks) && nrow(snp_peaks) > 0) {
    sp <- copy(snp_peaks)
    sp[, chrom := as.character(chrom)]
    sp[, pos_col := as.integer(snp_pos)]
    sp[, et_col  := normalize_edge_type_label(
      if ("haplotype_label" %in% names(sp)) haplotype_label else NA_character_
    )]
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

  # Peaks whose own edge_type is self-classifying (gene_conversion /
  # crossover / internal_crossover) but that never attached to any token —
  # e.g. an internal_crossover sitting in a SNP-desert with no flanking LOH
  # to anchor a motif rule on — are promoted here using the same
  # classify_tract() logic the motif rules use. A peak that resolves to a
  # real call (NCO_GC/CO_GC) IS an event, with real read support behind it,
  # so it's added to `events` and removed from `unclaimed_peaks` — showing
  # the same peak simultaneously as a called event *and* as "unclaimed" is
  # exactly the confusing double-reporting this reconciliation used to do.
  # They're named "_subres" (sub-resolution) and fall to "review" confidence
  # in build_event_table(), distinguishing them from the better-supported
  # CO_GC/NCO_GC calls fired by the motif rules. Peaks that classify_tract()
  # genuinely cannot resolve (binary singleton, independent_events, no
  # edge_type, low coverage, ...) are not events — they stay in
  # `unclaimed_peaks` for manual review, tagged with the reason.
  uncat_peak_events   <- list()
  still_unclaimed_pks <- list()
  for (u in unclaimed_peaks) {
    tract <- classify_tract(list(best_edge_type = u$edge_type,
                                 n_spanning = u$n_read_support), params)
    if (tract$call %in% c("NCO_GC", "CO_GC")) {
      # If phase-rescued, report the event over the excised island tract;
      # otherwise it is a point event at the peak's SNP position.
      u_isl_s <- u$phase_island_start
      u_isl_e <- u$phase_island_end
      ev_st <- if (!is.null(u_isl_s) && !is.na(u_isl_s)) as.integer(u_isl_s) else as.integer(u$snp_pos)
      ev_en <- if (!is.null(u_isl_e) && !is.na(u_isl_e)) as.integer(u_isl_e) else as.integer(u$snp_pos)
      uncat_peak_events <- c(uncat_peak_events, list(list(
        event_class = paste0(tract$call, "_subres"), chrom = u$chrom,
        start = ev_st, end = ev_en,
        length_bp = as.integer(ev_en - ev_st), n_support = tract$n_support,
        peak_edge_types = u$edge_type %||% NA_character_,
        phase_switch_frac = u$phase_switch_frac %||% NA_real_,
        notes = "no_fixed_tract; peak_only", tokens = list())))
    } else {
      u$reason <- tract$reason
      still_unclaimed_pks <- c(still_unclaimed_pks, list(u))
    }
  }
  unclaimed_peaks <- still_unclaimed_pks

  # LOH tokens that never attached to any motif-scanned event and have no
  # corroborating peak are genuinely unresolved — too short to raise a
  # detectable peak, sitting in a SNP desert, or just lacking chimeric-read
  # evidence either way. Unlike self-classifying peaks above, an orphan LOH
  # tract has no positive evidence for any specific mechanism, so it is NOT
  # promoted into `events` (that used to add UNCATEGORIZED_LOH/POSSIBLE_GC
  # rows to the events table while the *same* interval also sat in
  # `unclaimed_loh` — duplicate reporting of an un-scored interval as if it
  # were a called event). It's annotated with a reason and left only in
  # `unclaimed_loh` for manual review.
  unclaimed_loh <- lapply(unclaimed_loh, function(u) {
    n <- u$n_snps %||% NA_integer_
    is_small <- !is.na(n) && n < params$min_snps_for_peak
    u$reason <- if (is_small)
      paste0("n_snps=", n, " < min_snps_for_peak (", params$min_snps_for_peak,
             "); consider lowering the Minimum Run Length parameter")
    else
      "no attached peak"
    u
  })

  list(
    events          = c(all_events, uncat_peak_events),
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
      phase_switch_frac = numeric(),
      confidence = character(), notes = character()
    ))

  rows <- lapply(events, function(ev) {
    confidence <- if (ev$event_class %in% c("NCO_GC", "CO_GC", "CO_TERM",
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
      phase_switch_frac = ev$phase_switch_frac %||% NA_real_,
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
#'   $events          flat list of event objects. Self-classifying peaks
#'                    (gene_conversion/crossover/internal_crossover) that never
#'                    attached to a motif-scanned token are included here too
#'                    (as *_subres, "review" confidence) -- they ARE events,
#'                    just not shown twice by also appearing below.
#'   $unclaimed_peaks "Other events": peaks classify_tract() could not resolve
#'                    (binary singleton, independent_events, low coverage, ...),
#'                    for manual review. Not included in $events/$event_table.
#'   $unclaimed_loh   "Other events": LOH segments with no attached peak and no
#'                    motif-scanned event, for manual review. Not included in
#'                    $events/$event_table.
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
