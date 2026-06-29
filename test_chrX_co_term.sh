#!/usr/bin/env bash
# End-to-end regression test: S288C_chrX terminal crossover (CO_TERM)
#
# Reproduces the full Shiny-app pipeline on test_data/full.csv.gz and checks
# that the chimeric peak at ~510.6 KB + terminal LOH on chrX is labelled CO_TERM.
#
# Run from the repo root:
#   bash test_chrX_co_term.sh
# or with a different Rscript binary:
#   RSCRIPT=/path/to/Rscript bash test_chrX_co_term.sh

set -euo pipefail
RSCRIPT="${RSCRIPT:-Rscript}"
cd "$(dirname "$0")"

"$RSCRIPT" - <<'EOF'
suppressPackageStartupMessages({
  library(data.table)
})

source("chimera_functions.R")
source("loh_chain_analysis.R")

READ_DATA <- "test_data/full.csv.gz"
SNP_DATA  <- "test_data/full_SNPs.vcf.gz"
CHR_SIZE  <- "test_data/ref_genome.fasta.fai"
TARGET_CHR <- "S288C_chrX"

message("=== Step 1: run_chimera_analysis ===")
res <- run_chimera_analysis(
  read_data_path  = READ_DATA,
  snp_data_path   = SNP_DATA,
  chr_size_path   = CHR_SIZE,
  sample_name     = "test_chrX"
)

message("=== Step 2: compute_loh_map ===")
loh_out <- compute_loh_map(res$full_read_loh)

message("=== Step 3: compute_coverage_map ===")
cov_out <- compute_coverage_map(res$full_read_loh)

message("=== Step 4: compute_peak_pairs ===")
fusion_res <- compute_peak_pairs(
  snp_peaks      = res$snp_peaks,
  rt_df          = res$rt_df,
  transition_pos = res$transition_pos,
  loh_segments   = loh_out$loh_segments
)

message("=== Step 5: run_chain_analysis ===")
params <- default_chain_params()
chain_res <- run_chain_analysis(
  loh_segments      = loh_out$loh_segments,
  fused_peaks       = fusion_res$fused_peaks,
  peak_pairs        = fusion_res$peak_pairs,
  snp_peaks         = res$snp_peaks,
  rt_df             = res$rt_df,
  chr_span          = res$chr_span,
  coverage_segments = cov_out$coverage_segments,
  coverage_table    = cov_out$coverage_table,
  params            = params
)

# ── Report ────────────────────────────────────────────────────────────────────
et <- chain_res$event_table
chrX_events <- et[chrom == TARGET_CHR]

message("\n=== Events on ", TARGET_CHR, " ===")
if (nrow(chrX_events) == 0) {
  message("  (none)")
} else {
  print(chrX_events[, .(event_class, start, end, confidence, notes)])
}

has_co_term <- nrow(chrX_events[event_class == "CO_TERM"]) > 0

if (has_co_term) {
  message("\nPASS: CO_TERM found on ", TARGET_CHR)
  quit(status = 0)
} else {
  message("\nFAIL: no CO_TERM on ", TARGET_CHR,
          " (events: ", paste(chrX_events$event_class, collapse = ", "), ")")
  quit(status = 1)
}
EOF
