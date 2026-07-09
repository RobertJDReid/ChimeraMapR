#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# undefined_phase_diag.R  —  PROTOTYPE diagnostic pass (option B)
#
# For every peak that classify_peak_haplotype() labels "undefined", locate a
# FIX (REF- or ALT-fixed) island bounded by HET inside the classification
# window, then RE-CLASSIFY it from per-read PHASING across that island rather
# than population consensus. Consensus sees HET on both flanks (the two
# reciprocal crossover orientations average to ~0.5) and can never resolve
# the event; individual reads do.
#
# For each undefined peak it reports the left-flank vs right-flank majority-
# allele cross-tab of all reads spanning the island and a proposed call:
#   crossover (reciprocal) — flanks switch per-read, both directions, ~0 non-switch
#   crossover (one-sided)  — flanks switch per-read, single direction
#   gene_conversion (NCO)  — flanks AGREE per-read, fixed patch in middle
#   ambiguous / low_support / no_het_bounded_fix
#
# This is a READ-ONLY diagnostic. It does not modify any pipeline output.
#
# Usage:
#   Rscript undefined_phase_diag.R <read_data> <snp_data> <chr_size.fai>
# Defaults to the bundled test data if no args are given.
# ---------------------------------------------------------------------------

suppressMessages(library(data.table))
source("chimera_functions.R")

args <- commandArgs(trailingOnly = TRUE)
read_path <- if (length(args) >= 1) args[[1]] else "test_data/R187E_13.csv.gz"
snp_path  <- if (length(args) >= 2) args[[2]] else "test_data/full_SNPs.vcf.gz"
fai_path  <- if (length(args) >= 3) args[[3]] else "test_data/ref_genome.fasta.fai"

MIN_RUN     <- 2L    # zone_min_snps (matches CLI --min-run default)
MIN_SPAN    <- 5L    # min reads spanning the island to make a call
SWITCH_HI   <- 0.80  # >= this fraction switching per-read  -> crossover
SWITCH_LO   <- 0.20  # <= this fraction switching per-read  -> gene conversion

# ── Run the standard pipeline far enough to get labelled peaks ──────────────
message("Running analysis ...")
res <- run_chimera_analysis(read_path, snp_path, fai_path,
                            mapq_cutoff = 20L, baseq_cutoff = 10L,
                            min_run = MIN_RUN, min_peak_height = 10L, lambda = 1)
res$snp_peaks <- label_snp_peaks_haplotypes(
  res$snp_peaks, res$rt_df, res$transition_pos, zone_min_snps = MIN_RUN)

frl  <- res$full_read_loh          # ALL reads (MAPQ+is_del filtered): chrom,pos,read_id,IS_REF,ALLELE
peaks <- res$snp_peaks[!is.na(snp_pos)]
undef <- peaks[haplotype_label == "undefined"]
message(sprintf("Peaks: %d total, %d undefined", nrow(peaks), nrow(undef)))

# ── Locate a HET-bounded FIX island inside the classifier's seg_data ────────
# seg_data runs are in kb; xmax is the lead of the next run's xmin, so each
# run i spans [xmin_i, xmax_i) bp = [xmin_i*1000, xmin_{i+1}*1000).
find_het_bounded_fix <- function(seg) {
  if (is.null(seg) || nrow(seg) < 3) return(NULL)
  s <- as.character(seg$SNP_call)
  cand <- which(s %in% c("REF", "ALT") &
                shift(s, 1) == "HET" &
                shift(s, -1) == "HET")
  if (length(cand) == 0) return(NULL)
  # If several, take the widest island (most likely the real conversion tract)
  widths <- (seg$xmax[cand] - seg$xmin[cand])
  i <- cand[which.max(widths)]
  list(state = s[i],
       start_bp = seg$xmin[i] * 1000,
       end_bp   = seg$xmax[i] * 1000)
}

# ── Per-read L/R phasing across the island (all spanning reads) ─────────────
phase_across_island <- function(chr_name, win_start, win_end,
                                 isl_start, isl_end) {
  fr <- frl[as.character(chrom) == chr_name & pos >= win_start & pos <= win_end]
  if (nrow(fr) == 0) return(NULL)
  ph <- fr[, .(
    L = classify_zone_state(pos, IS_REF, win_start,    isl_start - 1L, MIN_RUN),
    R = classify_zone_state(pos, IS_REF, isl_end + 1L, win_end,        MIN_RUN)
  ), by = read_id][!is.na(L) & !is.na(R)]
  ph
}

classify_from_phase <- function(ph) {
  n <- nrow(ph)
  if (n < MIN_SPAN) return(list(call = "low_support", n = n))
  n_switch <- sum(ph$L != ph$R)
  frac     <- n_switch / n
  n_RtoA   <- sum(ph$L == "REF" & ph$R == "ALT")
  n_AtoR   <- sum(ph$L == "ALT" & ph$R == "REF")
  call <-
    if (frac >= SWITCH_HI) {
      if (n_RtoA > 0 && n_AtoR > 0) "crossover (reciprocal)" else "crossover (one-sided)"
    } else if (frac <= SWITCH_LO) {
      "gene_conversion (NCO)"
    } else "ambiguous"
  list(call = call, n = n, n_switch = n_switch, frac = frac,
       n_RtoA = n_RtoA, n_AtoR = n_AtoR)
}

# ── Diagnostic loop ─────────────────────────────────────────────────────────
rows <- list()
for (i in seq_len(nrow(undef))) {
  pk       <- undef[i]
  chr_name <- as.character(pk$chrom)
  touching <- res$transition_pos[
    as.character(chrom) == chr_name & pos >= pk$peak_start & pos <= pk$peak_end,
    unique(read_id)]
  hap <- classify_peak_haplotype(pk, chr_name, res$rt_df, touching, MIN_RUN)
  isl <- find_het_bounded_fix(hap$seg_data)

  base <- data.table(chrom = chr_name, snp_pos = pk$snp_pos,
                     win_start = hap$win_start, win_end = hap$win_end)
  if (is.null(isl)) {
    rows[[length(rows) + 1L]] <- cbind(base, data.table(
      island = "-", island_bp = "-", call = "no_het_bounded_fix",
      n_span = NA_integer_, n_switch = NA_integer_, n_RtoA = NA_integer_,
      n_AtoR = NA_integer_))
    next
  }
  ph <- phase_across_island(chr_name, hap$win_start, hap$win_end,
                            isl$start_bp, isl$end_bp)
  cl <- if (is.null(ph)) list(call = "low_support", n = 0) else classify_from_phase(ph)
  rows[[length(rows) + 1L]] <- cbind(base, data.table(
    island    = isl$state,
    island_bp = sprintf("%.1f-%.1fk", isl$start_bp/1000, isl$end_bp/1000),
    call      = cl$call,
    n_span    = cl$n,
    n_switch  = cl$n_switch %||% NA_integer_,
    n_RtoA    = cl$n_RtoA   %||% NA_integer_,
    n_AtoR    = cl$n_AtoR   %||% NA_integer_))
}

out <- rbindlist(rows, fill = TRUE)
cat("\n================ UNDEFINED-PEAK PHASE DIAGNOSTIC ================\n\n")
if (nrow(out) == 0) {
  cat("No undefined peaks.\n")
} else {
  setorder(out, chrom, snp_pos)
  print(out, nrow = 200)
  cat("\n--- proposed-call summary ---\n")
  print(out[, .N, by = call][order(-N)])
}
