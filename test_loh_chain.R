# =============================================================================
#  test_loh_chain.R
#
#  Self-contained test script for loh_chain_analysis.R
#
#  Usage:
#    Rscript test_loh_chain.R
#    — or —
#    source("test_loh_chain.R")   # from an interactive R session
#
#  Each test section prints PASS / FAIL. Final line reports total score.
#  No external test framework required; only data.table.
# =============================================================================

suppressMessages(library(data.table))

# %||% needed before loh_chain_analysis.R is sourced
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Locate loh_chain_analysis.R relative to this script, whether run via
# Rscript test_loh_chain.R  or  source("test_loh_chain.R")
.script_dir <- local({
  # Rscript: commandArgs gives --file=<path>
  argv <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", argv, value = TRUE)
  if (length(file_flag) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_flag[1])))
  } else {
    # Interactive source(): walk the call stack for an 'ofile' binding
    ofile <- NULL
    for (i in seq_len(sys.nframe())) {
      e <- sys.frame(i)
      if (exists("ofile", envir = e, inherits = FALSE)) {
        ofile <- get("ofile", envir = e, inherits = FALSE)
        break
      }
    }
    if (!is.null(ofile)) dirname(normalizePath(ofile)) else "."
  }
})

source(file.path(.script_dir, "loh_chain_analysis.R"), local = FALSE)

# -----------------------------------------------------------------------------
#  Test harness
# -----------------------------------------------------------------------------
.passed <- 0L
.failed <- 0L

expect <- function(label, expr) {
  result <- tryCatch(expr, error = function(e) {
    cat(sprintf("  FAIL  %s\n        ERROR: %s\n", label, conditionMessage(e)))
    .failed <<- .failed + 1L
    return(invisible(NULL))
  })
  if (isTRUE(result)) {
    cat(sprintf("  PASS  %s\n", label))
    .passed <<- .passed + 1L
  } else if (!is.null(result)) {
    cat(sprintf("  FAIL  %s\n        got: %s\n", label,
                paste(capture.output(print(result)), collapse = " ")))
    .failed <<- .failed + 1L
  }
}

section <- function(title) {
  cat(sprintf("\n── %s ──\n", title))
}

# Shorthand: pull the event_class vector from a run_chain_analysis result
event_classes <- function(res) res$event_table$event_class

# =============================================================================
#  SHARED FIXTURES
#  These are the minimal column sets each function requires.
# =============================================================================

# chr_span: one row per chromosome
# Columns: chrom (character), length (integer)
make_chr_span <- function(chroms, lengths) {
  data.table(chrom = chroms, length = as.integer(lengths))
}

# loh_segments: one row per contiguous same-state run
# Columns required by build_raw_chains:
#   chrom, start, end, length_bp, n_snps, loh_state,
#   balance_mean, balance_sd
make_loh_seg <- function(chrom, start, end, state,
                          n_snps = 20L, balance_mean = NULL, balance_sd = 0.01) {
  if (is.null(balance_mean))
    balance_mean <- if (state == "REF_fixed") 0.95 else
                   if (state == "ALT_fixed") 0.05 else 0.50
  data.table(
    chrom        = as.character(chrom),
    start        = as.integer(start),
    end          = as.integer(end),
    length_bp    = as.integer(end - start),
    n_snps       = as.integer(n_snps),
    loh_state    = state,
    balance_mean = balance_mean,
    balance_sd   = balance_sd
  )
}

loh_rbind <- function(...) rbindlist(list(...))

# fused_peaks: one row per peak (singleton fusion group = one sub-peak)
# Columns required by .get_chr_peaks:
#   chrom, fusion_group, fused_pos_bp, fused_start_bp, fused_end_bp,
#   n_sub_peaks, haplotype_label, best_edge_type
make_peak <- function(chrom, pos, start = NULL, end = NULL,
                      edge_type = "gene_conversion", n_sub = 1L,
                      hap_label = NA_character_, fusion_group = NULL,
                      window_half = 3000L) {
  if (is.null(start)) start <- pos - window_half
  if (is.null(end))   end   <- pos + window_half
  fg <- fusion_group %||% pos   # use position as a unique group id if not given
  data.table(
    chrom           = as.character(chrom),
    fusion_group    = as.integer(fg),
    fused_pos_bp    = as.integer(pos),
    fused_start_bp  = as.integer(start),
    fused_end_bp    = as.integer(end),
    n_sub_peaks     = as.integer(n_sub),
    haplotype_label = hap_label,
    best_edge_type  = edge_type,
    snp_pos         = as.integer(pos),   # needed by reconcile()
    peak_start      = as.integer(start),
    peak_end        = as.integer(end)
  )
}

peaks_rbind <- function(...) rbindlist(list(...))

# peak_pairs: one row per evaluated adjacent pair
# Columns required by .get_chr_peaks:
#   chrom, snp_pos_a, snp_pos_b, edge_type, n_spanning, jaccard
make_pair <- function(chrom, pos_a, pos_b,
                      edge_type  = "gene_conversion",
                      n_spanning = 5L,
                      jaccard    = 0.40) {
  data.table(
    chrom      = as.character(chrom),
    snp_pos_a  = as.integer(pos_a),
    snp_pos_b  = as.integer(pos_b),
    edge_type  = edge_type,
    n_spanning = as.integer(n_spanning),
    jaccard    = jaccard
  )
}

pairs_rbind <- function(...) rbindlist(list(...))

params <- default_chain_params()

# =============================================================================
#  SECTION 1 — TOKEN CONSTRUCTION
# =============================================================================
section("Token construction")

expect("make_token sets length_bp", {
  t <- make_token("F", state = "REF_fixed", start = 1000L, end = 50000L)
  t$length_bp == 49000L
})

expect("make_token meta field is an empty list by default", {
  t <- make_token("F", state = "REF_fixed", start = 1L, end = 1000L)
  is.list(t$meta) && length(t$meta) == 0L
})

expect("meta field accepts arbitrary content without error", {
  t <- make_token("F", state = "REF_fixed", start = 1L, end = 1000L,
                  meta = list(cn_state = 2L, ploidy_call = "diploid"))
  t$meta$cn_state == 2L && t$meta$ploidy_call == "diploid"
})

expect("make_chain stores schema_version", {
  ch <- make_chain(list(), "chrI", 500000L)
  ch$schema_version == CHAIN_SCHEMA_VERSION
})

# =============================================================================
#  SECTION 2 — CLASSIFY_TRACT
# =============================================================================
section("classify_tract")

expect("NULL peak -> AMBIGUOUS(no_peak)", {
  r <- classify_tract(NULL, params)
  r$call == "AMBIGUOUS" && r$reason == "no_peak"
})

expect("gene_conversion edge -> NCO_GC", {
  pk <- list(best_edge_type = "gene_conversion", n_spanning = 5L)
  r  <- classify_tract(pk, params)
  r$call == "NCO_GC"
})

expect("crossover edge -> CO_GC", {
  pk <- list(best_edge_type = "crossover", n_spanning = 5L)
  r  <- classify_tract(pk, params)
  r$call == "CO_GC"
})

expect("low spanning reads -> AMBIGUOUS(low_coverage)", {
  pk <- list(best_edge_type = "gene_conversion", n_spanning = 1L)
  r  <- classify_tract(pk, params)
  r$call == "AMBIGUOUS" && r$reason == "low_coverage"
})

expect("independent_events -> AMBIGUOUS(independent_reads)", {
  pk <- list(best_edge_type = "independent_events", n_spanning = 8L)
  r  <- classify_tract(pk, params)
  r$call == "AMBIGUOUS" && r$reason == "independent_reads"
})

expect("binary single peak -> AMBIGUOUS(binary_single_peak)", {
  pk <- list(best_edge_type = "binary", n_spanning = 6L)
  r  <- classify_tract(pk, params)
  r$call == "AMBIGUOUS" && r$reason == "binary_single_peak"
})

# =============================================================================
#  SECTION 3 — CLASSIFY_TWO_BINARY_JUNCTION
# =============================================================================
section("classify_two_binary_junction")

expect("two binary peaks with crossover pair -> CO_GC", {
  # pair_partner_pos must point at each other's own position — that's what
  # tells classify_two_binary_junction() this pair_edge_type genuinely
  # describes THIS junction rather than one peak's relationship to some
  # other, unrelated neighbor.
  pk_l <- list(best_edge_type = "binary", n_spanning = 4L, fused_pos_bp = 100000,
               pair_edge_type = "crossover", pair_partner_pos = 200000)
  pk_r <- list(best_edge_type = "binary", n_spanning = 4L, fused_pos_bp = 200000,
               pair_edge_type = "crossover", pair_partner_pos = 100000)
  r <- classify_two_binary_junction(pk_l, pk_r, "ALT_fixed", params)
  r$call == "CO_GC"
})

expect("two binary peaks with gene_conversion pair -> NCO_GC_LARGE", {
  pk_l <- list(best_edge_type = "binary", n_spanning = 4L, fused_pos_bp = 100000,
               pair_edge_type = "gene_conversion", pair_partner_pos = 200000)
  pk_r <- list(best_edge_type = "binary", n_spanning = 4L, fused_pos_bp = 200000,
               pair_edge_type = "gene_conversion", pair_partner_pos = 100000)
  r <- classify_two_binary_junction(pk_l, pk_r, "ALT_fixed", params)
  r$call == "NCO_GC_LARGE"
})

expect("NULL left peak -> AMBIGUOUS", {
  pk_r <- list(best_edge_type = "binary", n_spanning = 5L)
  r    <- classify_two_binary_junction(NULL, pk_r, "REF_fixed", params)
  r$call == "AMBIGUOUS"
})

expect("low spanning reads across both peaks -> AMBIGUOUS(low_coverage)", {
  pk_l <- list(best_edge_type = "binary", n_spanning = 1L)
  pk_r <- list(best_edge_type = "binary", n_spanning = 1L)
  r <- classify_two_binary_junction(pk_l, pk_r, "REF_fixed", params)
  r$call == "AMBIGUOUS" && r$reason == "low_coverage"
})

# =============================================================================
#  SECTION 4 — CANONICALISE
# =============================================================================
section("canonicalise")

expect("F G F (same state, small gap) -> merged to single F", {
  toks <- list(
    make_token("TEL", start = 1L, end = 1L),
    make_token("F", state = "REF_fixed", start = 100L,   end = 10000L,  n_snps = 20L, depth_ratio = 1.0),
    make_token("G", start = 10001L, end = 11000L),   # 999 bp < 5000 merge_gap
    make_token("F", state = "REF_fixed", start = 11001L, end = 50000L, n_snps = 30L, depth_ratio = 1.0),
    make_token("TEL", start = 500000L, end = 500000L)
  )
  ch  <- make_chain(toks, "chrI", 500000L)
  ch2 <- canonicalise(ch, params)
  types <- sapply(ch2$tokens, `[[`, "type")
  sum(types == "F") == 1L
})

expect("merged token sets bridged_gap = TRUE", {
  toks <- list(
    make_token("TEL", start = 1L, end = 1L),
    make_token("F", state = "REF_fixed", start = 100L,   end = 10000L,  n_snps = 10L, depth_ratio = 1.0),
    make_token("G", start = 10001L, end = 10500L),
    make_token("F", state = "REF_fixed", start = 10501L, end = 50000L, n_snps = 20L, depth_ratio = 1.0),
    make_token("TEL", start = 500000L, end = 500000L)
  )
  ch  <- make_chain(toks, "chrI", 500000L)
  ch2 <- canonicalise(ch, params)
  f_tok <- Filter(function(t) t$type == "F", ch2$tokens)[[1]]
  isTRUE(f_tok$bridged_gap)
})

expect("F G F different states -> NOT merged", {
  toks <- list(
    make_token("TEL", start = 1L, end = 1L),
    make_token("F", state = "REF_fixed", start = 100L,   end = 10000L,  n_snps = 20L, depth_ratio = 1.0),
    make_token("G", start = 10001L, end = 11000L),
    make_token("F", state = "ALT_fixed", start = 11001L, end = 50000L, n_snps = 30L, depth_ratio = 1.0),
    make_token("TEL", start = 500000L, end = 500000L)
  )
  ch  <- make_chain(toks, "chrI", 500000L)
  ch2 <- canonicalise(ch, params)
  types <- sapply(ch2$tokens, `[[`, "type")
  sum(types == "F") == 2L   # stays as two separate F tokens
})

expect("F G F gap larger than merge_gap -> NOT merged", {
  p2 <- default_chain_params()
  p2$merge_gap_bp <- 500L    # tight limit
  toks <- list(
    make_token("TEL", start = 1L, end = 1L),
    make_token("F", state = "REF_fixed", start = 100L, end = 10000L, n_snps = 10L, depth_ratio = 1.0),
    make_token("G", start = 10001L, end = 12000L),   # 1999 bp > 500 limit
    make_token("F", state = "REF_fixed", start = 12001L, end = 50000L, n_snps = 20L, depth_ratio = 1.0),
    make_token("TEL", start = 500000L, end = 500000L)
  )
  ch  <- make_chain(toks, "chrI", 500000L)
  ch2 <- canonicalise(ch, p2)
  types <- sapply(ch2$tokens, `[[`, "type")
  sum(types == "F") == 2L
})

expect("H G H (small gap) -> merged to single H", {
  toks <- list(
    make_token("TEL", start = 1L, end = 1L),
    make_token("H", state = "HET", start = 1000L,  end = 20000L, n_snps = 50L),
    make_token("G", start = 20001L, end = 21000L),
    make_token("H", state = "HET", start = 21001L, end = 80000L, n_snps = 60L),
    make_token("TEL", start = 500000L, end = 500000L)
  )
  ch  <- make_chain(toks, "chrI", 500000L)
  ch2 <- canonicalise(ch, params)
  types <- sapply(ch2$tokens, `[[`, "type")
  sum(types == "H") == 1L
})

# =============================================================================
#  SECTION 5 — BUILD_RAW_CHAINS
# =============================================================================
section("build_raw_chains")

expect("produces one chain per chromosome", {
  segs <- loh_rbind(
    make_loh_seg("chrI",  50000L, 150000L, "REF_fixed"),
    make_loh_seg("chrII", 80000L, 200000L, "ALT_fixed")
  )
  cs   <- make_chr_span(c("chrI","chrII"), c(500000L, 600000L))
  chs  <- build_raw_chains(segs, cs, params)
  length(chs) == 2L && all(c("chrI","chrII") %in% names(chs))
})

expect("chain starts and ends with TEL sentinels", {
  segs <- loh_rbind(make_loh_seg("chrI", 50000L, 150000L, "REF_fixed"))
  cs   <- make_chr_span("chrI", 500000L)
  ch   <- build_raw_chains(segs, cs, params)$chrI
  types <- sapply(ch$tokens, `[[`, "type")
  types[1] == "TEL" && types[length(types)] == "TEL"
})

expect("left TEL absorbs leading unscored region (no G before first segment)", {
  # TEL sentinels are now placed at min/max SNP positions so that the first
  # and last called segments are directly adjacent to their sentinel.
  segs <- loh_rbind(make_loh_seg("chrI", 50000L, 150000L, "REF_fixed"))
  cs   <- make_chr_span("chrI", 500000L)
  ch   <- build_raw_chains(segs, cs, params)$chrI
  toks <- ch$tokens
  toks[[1]]$type == "TEL" && toks[[1]]$end == 49999L &&
  toks[[2]]$type == "F"   # F is directly adjacent, no G between them
})

expect("peak attached to overlapping F token as peak_over", {
  segs <- loh_rbind(make_loh_seg("chrI", 50000L, 150000L, "REF_fixed"))
  cs   <- make_chr_span("chrI", 500000L)
  pks  <- make_peak("chrI", pos = 100000L, edge_type = "gene_conversion")
  ch   <- build_raw_chains(segs, cs, params, fused_peaks = pks)$chrI
  f_toks <- Filter(function(t) t$type == "F", ch$tokens)
  !is.null(f_toks[[1]]$peak_over)
})

# =============================================================================
#  SECTION 6 — FULL PIPELINE: NCO_GC  (H F● H)
#
#  Chromosome layout (500 kb):
#    0─────5k : HET
#    5k──100k : HET   (one segment; peak at 75k inside it)
#    100k─200k: REF_fixed  ← event locus
#    200k─500k: HET
#
#  Peak at 150k (gene_conversion, 5 spanning reads) sits inside the LOH.
# =============================================================================
section("Full pipeline: NCO_GC (H [F●] H)")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)

  expect("NCO_GC event is called", {
    "NCO_GC" %in% event_classes(res)
  })

  expect("NCO_GC event spans the REF_fixed region", {
    ev <- res$event_table[event_class == "NCO_GC"]
    nrow(ev) >= 1L &&
      ev$start[1] <= 100000L && ev$end[1] >= 200000L
  })

  expect("NCO_GC has high confidence", {
    ev <- res$event_table[event_class == "NCO_GC"]
    ev$confidence[1] == "high"
  })

  expect("no unclaimed LOH after clean NCO_GC", {
    length(res$unclaimed_loh) == 0L
  })
}

# =============================================================================
#  SECTION 7 — FULL PIPELINE: CO_GC  (H [F●] H, crossover edge)
# =============================================================================
section("Full pipeline: CO_GC (H [F●] H, crossover)")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "ALT_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "crossover")
  cs  <- make_chr_span("chrI", 500000L)

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)

  expect("CO_GC event is called", {
    "CO_GC" %in% event_classes(res)
  })
}

# =============================================================================
#  SECTION 8 — FULL PIPELINE: TWO-BINARY-PEAK CO_GC
#
#  Layout:
#    0──100k:  HET
#    100k─200k: REF_fixed (the LOH tract)
#    200k─500k: HET
#
#  Two binary peaks, one at each junction (80k and 210k).
#  A peak_pair row links them with edge_type = "crossover".
# =============================================================================
section("Full pipeline: CO_GC from two binary peaks")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 95000L,  edge_type = "binary",
              start = 92000L, end = 98000L,  fusion_group = 1L),
    make_peak("chrI", pos = 205000L, edge_type = "binary",
              start = 202000L, end = 208000L, fusion_group = 2L)
  )
  prs <- make_pair("chrI", pos_a = 95000L, pos_b = 205000L,
                   edge_type = "crossover", n_spanning = 6L)
  cs  <- make_chr_span("chrI", 500000L)

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             peak_pairs = prs, chr_span = cs, params = params)

  expect("CO_GC called from two binary peaks", {
    any(grepl("CO_GC", event_classes(res)))
  })
}

# =============================================================================
#  SECTION 9 — FULL PIPELINE: CROSSOVER_NO_TRACT
#
#  Layout:
#    0─────200k: REF_fixed  (extends from near telomere)
#    200k──500k: ALT_fixed
#
#  Adjacent opposite-fixed blocks with a crossover peak at the junction.
#  No interstitial HET or fixed fragment — clean crossover, no conversion tract.
# =============================================================================
section("Full pipeline: CROSSOVER_NO_TRACT")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 200000L, "REF_fixed", n_snps = 400L),
    make_loh_seg("chrI", 200001L, 500000L, "ALT_fixed", n_snps = 600L)
  )
  # Peak sits right at the junction
  pks <- make_peak("chrI", pos = 200000L,
                   start = 197000L, end = 203000L,
                   edge_type = "crossover")
  cs  <- make_chr_span("chrI", 500000L)

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)

  expect("CROSSOVER_NO_TRACT event is called", {
    "CROSSOVER_NO_TRACT" %in% event_classes(res)
  })
}

# =============================================================================
#  SECTION 10 — FULL PIPELINE: TERMINAL_LOH
#
#  Layout:
#    0──150k:  REF_fixed  (starts at chromosome start -> telomeric)
#    150k─500k: HET
#
#  Peak at 148k (at the internal end of the fixed region).
#  Depth ratio = 1.0 (not a deletion).
# =============================================================================
section("Full pipeline: TERMINAL_LOH")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 150000L, "REF_fixed", n_snps = 300L),
    make_loh_seg("chrI", 150001L, 500000L, "HET",       n_snps = 700L)
  )
  # Terminal LOH boundary is a single haplotype switch -> "binary" peak
  pks <- make_peak("chrI", pos = 148000L,
                   start = 145000L, end = 151000L,
                   edge_type = "binary")
  cs  <- make_chr_span("chrI", 500000L)

  # Use a depth_drop threshold that the segment's normal depth won't trigger
  p2 <- default_chain_params()
  p2$depth_drop <- 0.3   # only flag very low depth as deletion

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = p2)

  expect("TERMINAL_LOH event is called", {
    "TERMINAL_LOH" %in% event_classes(res)
  })

  expect("gene_conversion peak at terminal position -> AMBIGUOUS", {
    pks_gc <- make_peak("chrI", pos = 148000L,
                        start = 145000L, end = 151000L,
                        edge_type = "gene_conversion")
    res_gc <- run_chain_analysis(loh_segments = loh, fused_peaks = pks_gc,
                                 chr_span = cs, params = p2)
    any(grepl("AMBIGUOUS", event_classes(res_gc)))
  })
}

# =============================================================================
#  SECTION 11 — FULL PIPELINE: TERMINAL_DELETION
#
#  Same layout as above but depth_ratio is below depth_drop threshold.
#  We simulate low depth by setting depth_drop high so the normal density
#  triggers it.
# =============================================================================
section("Full pipeline: TERMINAL_DELETION (low depth)")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 150000L, "REF_fixed", n_snps = 5L),   # very few SNPs -> low density
    make_loh_seg("chrI", 150001L, 500000L, "HET",       n_snps = 700L)
  )
  cs  <- make_chr_span("chrI", 500000L)

  # Set depth_drop = 0.99 so almost any segment density triggers deletion
  p2 <- default_chain_params()
  p2$depth_drop <- 0.99

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = NULL,
                             chr_span = cs, params = p2)

  expect("TERMINAL_DELETION event is called", {
    "TERMINAL_DELETION" %in% event_classes(res)
  })
}

# =============================================================================
#  SECTION 12 — FULL PIPELINE: NCO_GC_in_terminal  (R06 opp_sandwich)
#
#  Layout: a tiny ALT_fixed fragment sandwiched inside two large REF_fixed
#  blocks. The tiny fragment is < SMALL_FRAC * sum(flanks).
#  After the event fires, the flanks are merged into one REF_fixed token.
#
#    0────50k:  REF_fixed  (50 kb)
#    50k──60k:  ALT_fixed  (10 kb = tiny; 10 / (50+390) = 0.023 < 0.10)
#    60k─450k:  REF_fixed  (390 kb)
#    450k─500k: HET
# =============================================================================
section("Full pipeline: NCO_GC_in_terminal (opp_sandwich)")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  50000L, "REF_fixed", n_snps = 100L),
    make_loh_seg("chrI",  50001L,  60000L, "ALT_fixed", n_snps = 8L),
    make_loh_seg("chrI",  60001L, 450000L, "REF_fixed", n_snps = 780L),
    make_loh_seg("chrI", 450001L, 500000L, "HET",       n_snps = 100L)
  )
  pks <- make_peak("chrI", pos = 55000L,
                   start = 52000L, end = 58000L,
                   edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)

  expect("NCO_GC_in_terminal event is called", {
    "NCO_GC_in_terminal" %in% event_classes(res)
  })

  expect("flanking REF_fixed tokens merged (no unclaimed LOH)", {
    length(res$unclaimed_loh) == 0L
  })
}

# =============================================================================
#  SECTION 13 — FULL PIPELINE: DOUBLE_GC  (H F H F̄ H)
#
#  Layout:
#    0────50k:  HET
#    50k─100k:  REF_fixed  (first conversion)
#    100k─150k: HET
#    150k─200k: ALT_fixed  (second conversion)
#    200k─500k: HET
#
#  Two peaks, one per conversion tract, both gene_conversion.
# =============================================================================
section("Full pipeline: DOUBLE_GC")

{
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  49999L, "HET",       n_snps = 100L),
    make_loh_seg("chrI",  50000L, 100000L, "REF_fixed", n_snps = 50L),
    make_loh_seg("chrI", 100001L, 149999L, "HET",       n_snps = 100L),
    make_loh_seg("chrI", 150000L, 200000L, "ALT_fixed", n_snps = 50L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos =  75000L, edge_type = "gene_conversion", fusion_group = 1L),
    make_peak("chrI", pos = 175000L, edge_type = "gene_conversion", fusion_group = 2L)
  )
  cs  <- make_chr_span("chrI", 500000L)

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)

  expect("DOUBLE_GC event is called", {
    "DOUBLE_GC" %in% event_classes(res)
  })

  expect("no unclaimed LOH after DOUBLE_GC", {
    length(res$unclaimed_loh) == 0L
  })
}

# =============================================================================
#  SECTION 14 — CANONICALISE REWRITE FREES DOUBLE_GC
#
#  Two same-state F runs split by a small gap, with a tiny opposite-state
#  fragment also present. Canonicalise should merge the same-state runs first,
#  then the motif scanner sees the correct opp_sandwich pattern.
# =============================================================================
section("Canonicalise enables downstream motif match")

{
  # REF (50k) | gap (2k) | REF (300k) with a tiny ALT (8k) in the middle:
  # Before merge:  REF | G | REF | ALT_tiny | REF | HET
  # After merge:   REF (big) | ALT_tiny | REF (still there) | HET
  # -> should fire opp_sandwich
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  50000L, "REF_fixed", n_snps = 100L),
    # Gap 50001-52000 (no segment — will become G token in build_raw_chains)
    make_loh_seg("chrI",  52001L, 200000L, "REF_fixed", n_snps = 300L),
    make_loh_seg("chrI", 200001L, 210000L, "ALT_fixed", n_snps = 8L),
    make_loh_seg("chrI", 210001L, 450000L, "REF_fixed", n_snps = 480L),
    make_loh_seg("chrI", 450001L, 500000L, "HET",       n_snps = 100L)
  )
  pks <- make_peak("chrI", pos = 205000L,
                   start = 202000L, end = 208000L,
                   edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)

  # merge_gap must be >= 2000 to bridge the 50001-52000 gap
  p2 <- default_chain_params()
  p2$merge_gap_bp <- 3000L

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = p2)

  expect("opp_sandwich fires after gap merge", {
    "NCO_GC_in_terminal" %in% event_classes(res)
  })
}

# =============================================================================
#  SECTION 15 — MULTI-CHROMOSOME
# =============================================================================
section("Multi-chromosome run")

{
  loh <- loh_rbind(
    # chrI: NCO_GC
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L),
    # chrII: TERMINAL_LOH
    make_loh_seg("chrII",     1L, 150000L, "ALT_fixed", n_snps = 300L),
    make_loh_seg("chrII", 150001L, 600000L, "HET",      n_snps = 900L)
  )
  pks <- peaks_rbind(
    make_peak("chrI",  pos = 150000L, edge_type = "gene_conversion"),
    make_peak("chrII", pos = 148000L, start = 145000L, end = 151000L,
              edge_type = "binary")   # terminal boundary -> binary peak
  )
  cs <- make_chr_span(c("chrI","chrII"), c(500000L, 600000L))

  p2 <- default_chain_params()
  p2$depth_drop <- 0.3

  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = p2)

  expect("events found on both chromosomes", {
    chroms_with_events <- unique(res$event_table$chrom)
    "chrI" %in% chroms_with_events && "chrII" %in% chroms_with_events
  })

  expect("NCO_GC on chrI", {
    any(res$event_table$chrom == "chrI" &
        res$event_table$event_class == "NCO_GC")
  })

  expect("TERMINAL_LOH on chrII", {
    any(res$event_table$chrom == "chrII" &
        res$event_table$event_class == "TERMINAL_LOH")
  })
}

# =============================================================================
#  SECTION 16 — UNCLAIMED TOKENS
# =============================================================================
section("Unclaimed token reconciliation")

expect("LOH with no peak -> not an event, but reported as other-event LOH", {
  loh <- loh_rbind(
    make_loh_seg("chrI",     1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = NULL,
                             chr_span = cs, params = params)
  !("UNCATEGORIZED_LOH" %in% event_classes(res)) &&
    length(res$unclaimed_loh) == 1L &&
    res$unclaimed_loh[[1]]$chrom == "chrI"
})

expect("self-classifying peak with no LOH -> promoted to NCO_GC_subres event, not left unclaimed", {
  loh <- loh_rbind(
    make_loh_seg("chrI", 1L, 500000L, "HET", n_snps = 1000L)
  )
  pks <- make_peak("chrI", pos = 250000L, edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "NCO_GC_subres" %in% event_classes(res) &&
    length(res$unclaimed_peaks) == 0L
})

expect("peak at junction with gene_conversion but low spanning -> AMBIGUOUS", {
  loh <- loh_rbind(
    make_loh_seg("chrI",     1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  # n_spanning = 1 < min_span = 3
  pks <- make_peak("chrI", pos = 150000L, edge_type = "gene_conversion")
  pks[, n_spanning := 1L]    # override the default NA -> explicit low value
  # also need to pass it through peak_pairs so .get_chr_peaks picks it up
  prs <- make_pair("chrI", pos_a = 100000L, pos_b = 200000L,
                   edge_type = "gene_conversion", n_spanning = 1L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             peak_pairs = prs, chr_span = cs, params = params)
  any(grepl("AMBIGUOUS", event_classes(res)))
})

# =============================================================================
#  SECTION 17 — EXTENSIBILITY: meta field round-trip
# =============================================================================
section("Extensibility: meta field")

expect("meta content survives canonicalise merge", {
  toks <- list(
    make_token("TEL", start = 1L, end = 1L),
    make_token("F", state = "REF_fixed", start = 100L, end = 10000L,
               n_snps = 20L, depth_ratio = 1.0,
               meta = list(cn_state = 2L, note = "left")),
    make_token("G", start = 10001L, end = 11000L),
    make_token("F", state = "REF_fixed", start = 11001L, end = 50000L,
               n_snps = 20L, depth_ratio = 1.0,
               meta = list(cn_state = 2L, note = "right")),
    make_token("TEL", start = 500000L, end = 500000L)
  )
  ch  <- make_chain(toks, "chrI", 500000L)
  ch2 <- canonicalise(ch, params)
  f   <- Filter(function(t) t$type == "F", ch2$tokens)[[1]]
  # Both meta lists should have been combined
  "cn_state" %in% names(f$meta)
})

expect("new meta field added post-construction survives to event table", {
  # Simulate adding a ploidy call to a token after chain construction
  loh <- loh_rbind(
    make_loh_seg("chrI",     1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)

  # Add ploidy to every F token in every chain (as a future layer would)
  for (cname in names(res$chains)) {
    ch <- res$chains[[cname]]
    for (ti in seq_along(ch$tokens)) {
      if (ch$tokens[[ti]]$type == "F") {
        ch$tokens[[ti]]$meta$ploidy_call <- "diploid"
        ch$tokens[[ti]]$meta$cn_state    <- 2L
      }
    }
    res$chains[[cname]] <- ch
  }

  # Verify meta is accessible post-assignment
  f_toks <- Filter(function(t) t$type == "F",
                   res$chains[["chrI"]]$tokens)
  length(f_toks) > 0 &&
    isTRUE(f_toks[[1]]$meta$ploidy_call == "diploid")
})

# =============================================================================
#  SECTION 18 — BUILD_EVENT_TABLE
# =============================================================================
section("build_event_table")

expect("empty event list produces zero-row data.table", {
  et <- build_event_table(list())
  is.data.table(et) && nrow(et) == 0L
})

expect("event table has required columns", {
  et <- build_event_table(list())
  all(c("event_class","chrom","start","end","length_kb",
        "n_support","peak_edge_types","confidence","notes") %in% names(et))
})

expect("NCO_GC gets confidence=high", {
  ev <- list(list(
    event_class = "NCO_GC", chrom = "chrI",
    start = 100000L, end = 200000L, length_bp = 100000L,
    n_support = 5L, peak_edge_types = "gene_conversion", notes = ""
  ))
  et <- build_event_table(ev)
  et$confidence[1] == "high"
})

expect("AMBIGUOUS gets confidence=review", {
  ev <- list(list(
    event_class = "AMBIGUOUS(low_coverage)", chrom = "chrI",
    start = 100000L, end = 200000L, length_bp = 100000L,
    n_support = 1L, peak_edge_types = NA_character_, notes = ""
  ))
  et <- build_event_table(ev)
  et$confidence[1] == "review"
})

# =============================================================================
#  SECTION 19 — REAL-DATA LAYOUT (no HET rows in loh_segments)
#
#  compute_loh_map() only emits REF_fixed and ALT_fixed rows; HET stretches
#  appear as coordinate gaps (G tokens) between them, not as explicit H rows.
#  All motif rules must fire on G-flanked F tokens, not just H-flanked ones.
# =============================================================================
section("Real-data layout: no HET rows in loh_segments")

expect("NCO_GC fires when F flanked by G tokens (no HET rows)", {
  # Only a REF_fixed segment; HET flanks are implicit coordinate gaps
  loh <- loh_rbind(
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "NCO_GC" %in% event_classes(res)
})

expect("CO_GC fires when F flanked by G tokens", {
  loh <- loh_rbind(
    make_loh_seg("chrI", 100000L, 200000L, "ALT_fixed", n_snps = 40L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "crossover")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "CO_GC" %in% event_classes(res)
})

expect("TERMINAL_LOH fires when F starts near telomere with non-fixed on internal side", {
  # TEL sentinels are placed at min/max SNP positions, so a HET (or G) segment
  # is needed on the non-terminal side to serve as the non-fixed context that
  # distinguishes TERMINAL_LOH from a full-chromosome LOH.
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 150000L, "REF_fixed", n_snps = 300L),
    make_loh_seg("chrI", 155000L, 400000L, "HET",       n_snps = 500L)
  )
  # Terminal boundary is a single haplotype switch -> binary peak
  pks <- make_peak("chrI", pos = 148000L, start = 145000L, end = 151000L,
                   edge_type = "binary")
  cs  <- make_chr_span("chrI", 500000L)
  p2  <- default_chain_params(); p2$depth_drop <- 0.3
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = p2)
  "TERMINAL_LOH" %in% event_classes(res)
})

expect("DOUBLE_GC fires when both F tokens flanked only by G tokens", {
  loh <- loh_rbind(
    make_loh_seg("chrI",  50000L, 100000L, "REF_fixed", n_snps = 50L),
    make_loh_seg("chrI", 150000L, 200000L, "ALT_fixed", n_snps = 50L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos =  75000L, edge_type = "gene_conversion", fusion_group = 1L),
    make_peak("chrI", pos = 175000L, edge_type = "gene_conversion", fusion_group = 2L)
  )
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "DOUBLE_GC" %in% event_classes(res)
})

expect("CROSSOVER_NO_TRACT fires on adjacent F tokens with no intervening HET", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 200000L, "REF_fixed", n_snps = 400L),
    make_loh_seg("chrI", 200001L, 500000L, "ALT_fixed", n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 200000L, start = 197000L, end = 203000L,
                   edge_type = "crossover")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "CROSSOVER_NO_TRACT" %in% event_classes(res)
})

# =============================================================================
#  SECTION 20 — R10: PEAK-DIRECT CLASSIFICATION
#
#  gene_conversion / crossover / internal_crossover peaks directly over F
#  should classify it as NCO_GC or CO_GC without needing flanking context.
# =============================================================================
section("R10: peak-direct classification")

expect("gene_conversion peak over F -> NCO_GC", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "NCO_GC" %in% event_classes(res)
})

expect("crossover peak over F -> CO_GC", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "ALT_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "crossover")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "CO_GC" %in% event_classes(res)
})

expect("internal_crossover peak over F -> CO_GC", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 100000L, "REF_fixed", n_snps = 1L),
    make_loh_seg("chrI", 100001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 100000L, start = 99800L, end = 100200L,
                   edge_type = "internal_crossover")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "CO_GC" %in% event_classes(res)
})

expect("binary peak over F does NOT trigger R10 (falls through to R11 or uncategorized)", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "binary")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  # R10 should NOT fire; event should be UNCATEGORIZED or AMBIGUOUS, not NCO_GC/CO_GC
  !any(c("NCO_GC","CO_GC") %in% event_classes(res))
})

expect("gene_conversion peak over small F (1 SNP) upgrades from POSSIBLE_GC to NCO_GC", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 100000L, "REF_fixed", n_snps = 1L),
    make_loh_seg("chrI", 100001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 100000L, start = 99800L, end = 100200L,
                   edge_type = "gene_conversion")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "NCO_GC" %in% event_classes(res)
})

expect("gene_conversion peak with low spanning reads -> AMBIGUOUS (coverage gate preserved)", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 150000L, edge_type = "gene_conversion")
  prs <- make_pair("chrI", pos_a = 100000L, pos_b = 200000L,
                   edge_type = "gene_conversion", n_spanning = 1L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             peak_pairs = prs, chr_span = cs, params = params)
  any(grepl("AMBIGUOUS", event_classes(res)))
})

# =============================================================================
#  SECTION 21 — R11: TWO BINARY FLANKING PEAKS
#
#  An F token with binary peaks at both the H→F and F→H junctions, each
#  found by scanning transparent G gaps.  pair_edge_type drives NCO/CO.
# =============================================================================
section("R11: two binary flanking peaks")

expect("two binary flanking peaks + crossover pair -> CO_GC", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 95000L,  edge_type = "binary",
              start = 92000L, end = 98000L, fusion_group = 1L),
    make_peak("chrI", pos = 205000L, edge_type = "binary",
              start = 202000L, end = 208000L, fusion_group = 2L)
  )
  prs <- make_pair("chrI", pos_a = 95000L, pos_b = 205000L,
                   edge_type = "crossover", n_spanning = 6L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             peak_pairs = prs, chr_span = cs, params = params)
  any(grepl("CO_GC", event_classes(res)))
})

expect("two binary flanking peaks + gene_conversion pair -> NCO_GC", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "ALT_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 95000L,  edge_type = "binary",
              start = 92000L, end = 98000L, fusion_group = 1L),
    make_peak("chrI", pos = 205000L, edge_type = "binary",
              start = 202000L, end = 208000L, fusion_group = 2L)
  )
  prs <- make_pair("chrI", pos_a = 95000L, pos_b = 205000L,
                   edge_type = "gene_conversion", n_spanning = 6L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             peak_pairs = prs, chr_span = cs, params = params)
  any(grepl("NCO_GC", event_classes(res)))
})

expect("two binary peaks with no pair data -> GC_UNRESOLVED (confirmed GC, NCO/CO undetermined)", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 95000L,  edge_type = "binary",
              start = 92000L, end = 98000L, fusion_group = 1L),
    make_peak("chrI", pos = 205000L, edge_type = "binary",
              start = 202000L, end = 208000L, fusion_group = 2L)
  )
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)  # no peak_pairs
  any(grepl("GC_UNRESOLVED", event_classes(res)))
})

expect("R11 finds junction peaks across G gaps (real-data layout)", {
  # HET segments flank the F token with coordinate G gaps between them.
  # Peaks attach to H (not G) because .attach_peaks is called on H/F tokens
  # but not on gap tokens.  Peaks should be inside H but close to the H-G
  # boundary (within merge_gap_bp of the F token's edge — .left_junction_peak/
  # .right_junction_peak now bound their *_over fallback to that tolerance,
  # so a peak buried deep inside a long neighbouring H run isn't mistaken for
  # a different, far-away junction's peak) so they are found as left/right
  # junction peaks by R11.
  loh <- loh_rbind(
    make_loh_seg("chrI",  50000L,  89999L, "HET",       n_snps = 100L),
    make_loh_seg("chrI",  92000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 202000L, 500000L, "HET",       n_snps = 300L)
  )
  # Peaks near the H boundaries — inside H, reflecting the haplotype switch
  # position at the H→F junction that real chimera reads detect.
  pks <- peaks_rbind(
    make_peak("chrI", pos = 89500L, edge_type = "binary",
              start = 87000L, end = 91000L, fusion_group = 1L),   # inside H(50000-89999)
    make_peak("chrI", pos = 202500L, edge_type = "binary",
              start = 200500L, end = 204000L, fusion_group = 2L)  # inside H(202000-500000)
  )
  prs <- make_pair("chrI", pos_a = 89500L, pos_b = 202500L,
                   edge_type = "gene_conversion", n_spanning = 5L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             peak_pairs = prs, chr_span = cs, params = params)
  "NCO_GC" %in% event_classes(res)
})

expect("single binary peak (only left junction) does NOT trigger R11", {
  # F with a peak on the left side but nothing on the right
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 95000L, edge_type = "binary",
                   start = 92000L, end = 98000L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  !any(c("NCO_GC","CO_GC") %in% event_classes(res))
})

expect("single binary peak, other side genuinely peak-less -> GC_ONE_SIDED (R11b)", {
  # Same layout as the R11 negative test above: F flanked by a confirmed
  # binary peak on the left and NO peak at all on the right (e.g. the
  # reference is unmappable/repetitive there, so no SNP calls exist near
  # that boundary). R11 correctly declines (needs both sides); R11b should
  # report the one confirmed peak instead of silently dropping it to
  # UNCATEGORIZED_LOH.
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 200000L, "REF_fixed", n_snps = 40L),
    make_loh_seg("chrI", 200001L, 500000L, "HET",       n_snps = 600L)
  )
  pks <- make_peak("chrI", pos = 95000L, edge_type = "binary",
                   start = 92000L, end = 98000L)
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "GC_ONE_SIDED" %in% event_classes(res)
})

expect("a peak buried far inside a long neighbouring H run is NOT misattributed as this F's junction peak", {
  # Regression test for the RAD5_6 chrXV bug: a small REF_fixed island
  # (100000-103000) sits next to a long HET run (110000-800000) whose only
  # detected peak (750000) is far past merge_gap_bp from the island's right
  # edge (103000) -- it actually marks a DIFFERENT junction much further
  # along the same H run (a separate small F island further out, so the long
  # H run's far end is NOT itself telomeric -- that would pull in R02c
  # instead of the R11/R11b junction-peak logic this test targets).
  # R11/R11b must not borrow that distant peak for the first island's right
  # junction: with only a real left peak and no legitimate right peak, this
  # should resolve as GC_ONE_SIDED (from the left peak alone), never
  # NCO_GC/CO_GC (which would require pairing with the wrongly-borrowed
  # distant peak).
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L,  99999L, "HET",       n_snps = 200L),
    make_loh_seg("chrI", 100000L, 103000L, "REF_fixed", n_snps = 20L),
    make_loh_seg("chrI", 110000L, 800000L, "HET",       n_snps = 900L),
    make_loh_seg("chrI", 900000L, 905000L, "REF_fixed", n_snps = 20L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 95000L, edge_type = "binary",
              start = 92000L, end = 98000L, fusion_group = 1L),
    make_peak("chrI", pos = 750000L, edge_type = "binary",
              start = 748000L, end = 752000L, fusion_group = 2L)
  )
  cs  <- make_chr_span("chrI", 1000000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  ec <- event_classes(res)
  "GC_ONE_SIDED" %in% ec && !any(c("NCO_GC", "CO_GC", "GC_UNRESOLVED") %in% ec)
})

# =============================================================================
#  SECTION 22 — R03: TCO_CAPTURED_TCO
#
#  Layout: H [F] F̄→TEL  — a terminal crossover (ALT_fixed, reaching the
#  chromosome end) that itself captured an earlier terminal crossover
#  (REF_fixed), with a chimeric-read peak at each junction.
#
#    1────300k:    HET
#    300k─400k:    REF_fixed  (first, now-buried, terminal CO)
#    400k─500k:    ALT_fixed  (second, capturing terminal CO; reaches TEL)
# =============================================================================
section("R03: TCO_CAPTURED_TCO")

expect("TCO_CAPTURED_TCO fires for H-F-F-bar(telomeric) with peaks at both junctions", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 300000L, "HET",       n_snps = 300L),
    make_loh_seg("chrI", 300001L, 400000L, "REF_fixed", n_snps = 50L),
    make_loh_seg("chrI", 400001L, 500000L, "ALT_fixed", n_snps = 50L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 300000L, start = 298000L, end = 302000L,
              edge_type = "binary", fusion_group = 1L),
    make_peak("chrI", pos = 400000L, start = 398000L, end = 402000L,
              edge_type = "binary", fusion_group = 2L)
  )
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "TCO_CAPTURED_TCO" %in% event_classes(res)
})

expect("ordinary single terminal F (no captured TCO) still calls CO_TERM/TERMINAL_LOH, not TCO_CAPTURED_TCO", {
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 300000L, "HET",       n_snps = 300L),
    make_loh_seg("chrI", 300001L, 500000L, "ALT_fixed", n_snps = 100L)
  )
  pks <- make_peak("chrI", pos = 300000L, start = 298000L, end = 302000L,
                   edge_type = "binary")
  cs  <- make_chr_span("chrI", 500000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  !("TCO_CAPTURED_TCO" %in% event_classes(res))
})

expect("TCO_CAPTURED_TCO fires once the outer F is extended to the telomere by an R06 merge (real-data layout: G gaps + boundary-exact peaks)", {
  # Mirrors a real chrXV chain: H, then a captured REF_fixed island, then an
  # ALT_fixed run that only reaches the telomere after R06 (opp_sandwich)
  # peels a tiny embedded REF fragment further out and merges its flanks.
  # Both junction peaks sit exactly AT a segment boundary (not within
  # peak_pad_bp of it) and are separated from their nearest token by a G
  # gap wider than peak_pad_bp, so they only ever attach as peak_over —
  # exercising both the R06 peak_over-preservation fix and the
  # boundary-inclusive (<=/>=) junction-peak fix.
  loh <- loh_rbind(
    make_loh_seg("chrI",      1L, 300000L, "HET",       n_snps = 300L),
    make_loh_seg("chrI", 301001L, 400000L, "REF_fixed", n_snps = 50L),
    make_loh_seg("chrI", 401001L, 500000L, "ALT_fixed", n_snps = 50L),
    make_loh_seg("chrI", 501001L, 501001L, "REF_fixed", n_snps = 1L),
    make_loh_seg("chrI", 502001L, 600000L, "ALT_fixed", n_snps = 50L)
  )
  pks <- peaks_rbind(
    make_peak("chrI", pos = 301001L, start = 296000L, end = 306000L,
              edge_type = "binary", fusion_group = 1L),
    make_peak("chrI", pos = 400000L, start = 397000L, end = 403000L,
              edge_type = "binary", fusion_group = 2L)
  )
  cs  <- make_chr_span("chrI", 600000L)
  res <- run_chain_analysis(loh_segments = loh, fused_peaks = pks,
                             chr_span = cs, params = params)
  "TCO_CAPTURED_TCO" %in% event_classes(res)
})

# =============================================================================
#  SUMMARY
# =============================================================================
cat(sprintf(
  "\n══════════════════════════════════════\n  Results: %d passed, %d failed\n══════════════════════════════════════\n",
  .passed, .failed
))
if (.failed > 0L) {
  message("Some tests failed. See FAIL lines above for details.")
  quit(status = 1L, save = "no")
} else {
  message("All tests passed.")
}
