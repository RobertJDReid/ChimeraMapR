#!/usr/bin/env Rscript
# =============================================================================
#  chimera_cli.R  —  Command-line interface for ChimeraMapR
#
#  Usage:
#    Rscript chimera_cli.R [options] <read_data> <snp_data> <chr_size.fai>
#
#  Examples:
#    # Default: save overview PNG
#    Rscript chimera_cli.R reads.csv.gz snps.vcf.gz genome.fa.fai
#
#    # Named output path
#    Rscript chimera_cli.R --output results/my_overview.png \
#                          reads.csv.gz snps.vcf.gz genome.fa.fai
#
#    # Peak list CSV
#    Rscript chimera_cli.R --peak-list reads.csv.gz snps.vcf.gz genome.fa.fai
#
#    # Overview plot saved as an RDS object (for re-plotting in R)
#    Rscript chimera_cli.R --overview-rds reads.csv.gz snps.vcf.gz genome.fa.fai
#
#    # Per-position coverage table + collapsed coverage segments as CSV
#    Rscript chimera_cli.R --coverage-map reads.csv.gz snps.vcf.gz genome.fa.fai
#
#    # Only the final chain-caller events table as CSV
#    Rscript chimera_cli.R --events-table reads.csv.gz snps.vcf.gz genome.fa.fai
#
#    # Full parameter set
#    Rscript chimera_cli.R \
#      --sample-name My_Sample \
#      --mapq-cutoff 30 \
#      --baseq-cutoff 15 \
#      --min-run 3 \
#      --min-peak-height 10 \
#      --lambda 1 \
#      --output results/ \
#      reads.csv.gz snps.vcf.gz genome.fa.fai
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# ── Source shared analysis functions ─────────────────────────────────────────
# Resolve symlinks manually — normalizePath() does NOT follow symlinks on macOS
resolve_symlink <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  for (i in seq_len(20)) {          # guard against circular links
    target <- Sys.readlink(path)
    if (target == "") break         # not a symlink (or fully resolved)
    if (!startsWith(target, "/"))   # relative symlink → make absolute
      target <- file.path(dirname(path), target)
    path <- normalizePath(target, mustWork = FALSE)
  }
  path
}

script_dir <- tryCatch({
  argv     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", argv, value = TRUE)
  if (length(file_arg)) {
    real_path <- resolve_symlink(sub("^--file=", "", file_arg))
    dirname(real_path)
  } else {
    "."
  }
}, error = function(e) ".")

func_paths   <- c(
  file.path(script_dir, "R", "chimera_functions.R"),
  file.path(script_dir, "chimera_functions.R"),
  "R/chimera_functions.R",
  "chimera_functions.R"
)
found <- Find(file.exists, func_paths)
if (is.null(found))
  stop("Cannot find chimera_functions.R. Expected at: ",
       paste(func_paths, collapse = "\n  "))
source(found)

chain_paths <- c(
  file.path(script_dir, "R", "loh_chain_analysis.R"),
  file.path(script_dir, "loh_chain_analysis.R"),
  "R/loh_chain_analysis.R",
  "loh_chain_analysis.R"
)
chain_found <- Find(file.exists, chain_paths)
if (!is.null(chain_found)) source(chain_found)


# ── Option parser ─────────────────────────────────────────────────────────────
option_list <- list(

  # Analysis parameters
  make_option(c("-n", "--sample-name"),
              type    = "character",
              default = "Sample_01",
              metavar = "NAME",
              help    = "Sample name used in output filenames and plot titles [default: %default]"),

  make_option("--mapq-cutoff",
              type    = "integer",
              default = 20L,
              metavar = "INT",
              help    = "Minimum MAPQ value; reads below this are excluded [default: %default]"),

  make_option("--baseq-cutoff",
              type    = "integer",
              default = 10L,
              metavar = "INT",
              help    = "Minimum base quality at SNP position [default: %default]"),

  make_option("--del-rate-cutoff",
              type    = "double",
              default = 0.10,
              metavar = "FLOAT",
              help    = paste("SNPs where more than this fraction of confidently-mapped",
                              "reads register a deletion are excluded [default: %default]")),

  make_option("--min-run",
              type    = "integer",
              default = 2L,
              metavar = "INT",
              help    = "Minimum consecutive same-allele calls to count as a run [default: %default]"),

  make_option("--min-peak-height",
              type    = "integer",
              default = 10L,
              metavar = "INT",
              help    = "Minimum transition count at a SNP to qualify as a peak [default: %default]"),

  make_option(c("-l", "--lambda"),
              type    = "double",
              default = 1,
              metavar = "FLOAT",
              help    = "Whittaker smoothing penalty lambda (lower = tighter fit) [default: %default]"),

  # Output mode — mutually exclusive; default is PNG overview plot
  make_option("--peak-list",
              action  = "store_true",
              default = FALSE,
              help    = "Write a CSV of detected peaks instead of the overview plot"),

  make_option("--overview-rds",
              action  = "store_true",
              default = FALSE,
              help    = "Save the overview plot as an RDS object for re-plotting in R"),

  make_option("--events-table",
              action  = "store_true",
              default = FALSE,
              help    = paste("Run the chain-based caller and write ONLY the final events",
                              "table as a CSV (no intermediate step CSVs, no plot)")),

  # Coverage map (sequencing depth)
  make_option("--coverage-map",
              action  = "store_true",
              default = FALSE,
              help    = paste("Also write the per-position coverage table and collapsed",
                              "coverage segments from compute_coverage_map() as CSVs",
                              "(in addition to the primary output mode)")),

  # Chain-based LOH event calling
  make_option("--chain-all",
              action  = "store_true",
              default = FALSE,
              help    = paste("Run the chain-based LOH event caller and write",
                              "one CSV per pass (step 0–4).")),

  make_option("--tel-tol",
              type    = "integer",
              default = 5L,
              metavar = "KB",
              help    = "Telomere tolerance in kb for chain analysis [default: %default]"),

  make_option("--merge-gap",
              type    = "integer",
              default = 5L,
              metavar = "KB",
              help    = "Same-state NA-gap merge threshold in kb [default: %default]"),

  make_option("--min-span",
              type    = "integer",
              default = 3L,
              metavar = "INT",
              help    = "Minimum spanning reads required for a read-based chain call [default: %default]"),

  make_option("--peak-pad",
              type    = "integer",
              default = 200L,
              metavar = "BP",
              help    = paste("Peak association padding in bp: a peak is attached to a token",
                              "junction if its SNP position falls within this distance of the",
                              "token boundary [default: %default]")),

  # Output path
  make_option(c("-o", "--output"),
              type    = "character",
              default = NULL,
              metavar = "PATH",
              help    = paste("Output file path or directory.",
                              "If a directory (or omitted), a dated filename is auto-generated.",
                              "Default extension: .png / .csv / .rds depending on output mode."))
)

parser <- OptionParser(
  usage       = "%prog [options] <read_data> <snp_data> <chr_size.fai>",
  option_list = option_list,
  description = paste(
    "\nChimeraMapR command-line interface.",
    "Detects chimeric haplotypes in long-read sequencing data.",
    "\nThree positional arguments are required (in order):",
    "  1. Read data file  (.csv or .csv.gz)",
    "  2. SNP data file   (.csv, .vcf, or .vcf.gz)",
    "  3. Chromosome size file (.fai)",
    sep = "\n"
  ),
  epilogue = paste(
    "Output modes (pick one; default is PNG overview plot):",
    "  [default]       Genome-wide overview plot saved as PNG",
    "  --peak-list     Detected peaks as CSV",
    "  --overview-rds  Overview plot object saved as RDS (re-plot with readRDS + print)",
    "  --events-table  Final chain-caller events table as CSV (only the events)",
    sep = "\n"
  )
)


# ── Parse args ────────────────────────────────────────────────────────────────
parsed <- parse_args(parser, positional_arguments = 3)
opts   <- parsed$options
args   <- parsed$args

read_path <- args[1]
snp_path  <- args[2]
fai_path  <- args[3]

# Validate input files exist
for (f in c(read_path, snp_path, fai_path)) {
  if (!file.exists(f))
    stop("Input file not found: ", f)
}

# Validate output mode — the primary output modes are mutually exclusive
n_modes <- sum(opts[["peak-list"]], opts[["overview-rds"]], opts[["events-table"]])
if (n_modes > 1)
  stop("--peak-list, --overview-rds and --events-table are mutually exclusive. Choose one.")

# --events-table produces the chain caller's final events table only; asking for
# --chain-all at the same time (which writes every step CSV + the annotated PNG)
# is contradictory.
if (opts[["events-table"]] && opts[["chain-all"]])
  stop("--events-table and --chain-all are mutually exclusive. Use --chain-all for the full output set.")

# Both the full chain run and the events-only mode need the chain analysis code.
run_chain <- opts[["chain-all"]] || opts[["events-table"]]
if (run_chain && is.null(chain_found))
  stop("--chain-all / --events-table requires loh_chain_analysis.R but it was not found next to chimera_functions.R")


# ── Resolve output path ───────────────────────────────────────────────────────
resolve_output <- function(opts_output, sample_name, mode) {
  ext <- switch(mode,
    png    = ".png",
    csv    = "_peaks.csv",
    rds    = "_overview.rds",
    events = "_events.csv"
  )
  stem <- paste0(sample_name, "_chimera_", Sys.Date())

  if (is.null(opts_output)) {
    # Default: write to current directory
    return(paste0(stem, ext))
  }

  if (dir.exists(opts_output) || grepl("/$", opts_output)) {
    # It's a directory: auto-name the file inside it
    dir.create(opts_output, recursive = TRUE, showWarnings = FALSE)
    return(file.path(opts_output, paste0(stem, ext)))
  }

  # Treat as a full file path — create parent directory if needed
  dir.create(dirname(opts_output), recursive = TRUE, showWarnings = FALSE)
  opts_output
}

# Resolve a directory for "extra" per-run outputs (coverage map, chain CSVs)
# that sit alongside the main --output, independent of its mode/extension.
resolve_extra_dir <- function(opts_output, out_path) {
  if (!is.null(opts_output) && (dir.exists(opts_output) || grepl("/$", opts_output))) {
    opts_output
  } else {
    dirname(normalizePath(out_path, mustWork = FALSE))
  }
}

output_mode <- if (opts[["peak-list"]])    "csv"    else
               if (opts[["overview-rds"]]) "rds"    else
               if (opts[["events-table"]]) "events" else
                                           "png"

out_path <- resolve_output(opts[["output"]], opts[["sample-name"]], output_mode)


# ── Run analysis ──────────────────────────────────────────────────────────────
cat("ChimeraMapR", APP_VERSION, "—", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Sample:", opts[["sample-name"]], "\n")
cat("Inputs:\n")
cat("  Read data : ", read_path, "\n")
cat("  SNP data  : ", snp_path,  "\n")
cat("  FAI       : ", fai_path,  "\n")
cat("Parameters:\n")
cat("  MAPQ cutoff     :", opts[["mapq-cutoff"]],     "\n")
cat("  BaseQ cutoff    :", opts[["baseq-cutoff"]],    "\n")
cat("  Del rate cutoff :", opts[["del-rate-cutoff"]], "\n")
cat("  Min run         :", opts[["min-run"]],         "\n")
cat("  Min peak height :", opts[["min-peak-height"]], "\n")
cat("  Lambda (λ)      :", opts[["lambda"]],          "\n")
cat("Output mode       :", output_mode, "→", out_path, "\n")
if (opts[["coverage-map"]]) cat("Coverage map      : enabled\n")
if (opts[["chain-all"]]) {
  cat("Chain analysis    : enabled\n")
  cat("  Telomere tol.   :", opts[["tel-tol"]],   "kb\n")
  cat("  Merge gap       :", opts[["merge-gap"]], "kb\n")
  cat("  Min spanning    :", opts[["min-span"]],  "reads\n")
  cat("  Peak pad        :", opts[["peak-pad"]],  "bp\n")
}
cat("\n")

results <- run_chimera_analysis(
  read_data_path  = read_path,
  snp_data_path   = snp_path,
  chr_size_path   = fai_path,
  sample_name     = opts[["sample-name"]],
  mapq_cutoff     = opts[["mapq-cutoff"]],
  baseq_cutoff    = opts[["baseq-cutoff"]],
  del_rate_cutoff = opts[["del-rate-cutoff"]],
  min_run         = opts[["min-run"]],
  min_peak_height = opts[["min-peak-height"]],
  lambda          = opts[["lambda"]]
)

ploidy_map          <- get_chromosome_ploidy(results$full_read_loh)
results$ploidy_map  <- ploidy_map

# Classify each peak's haplotype run-pattern (binary / gene_conversion /
# internal_crossover / undefined) from its own read evidence, before any
# fusion or event calling -- the same per-peak classification shown on the
# app's individual peak plots. Included in --peak-list / --overview-rds, and
# reused (not recomputed) by --chain-all's compute_peak_pairs() below.
results$snp_peaks <- label_snp_peaks_haplotypes(
  results$snp_peaks, results$rt_df, results$transition_pos,
  zone_min_snps = opts[["min-run"]]
)


# ── Print run summary ─────────────────────────────────────────────────────────
n_chimeric  <- length(results$chimeric_read_ids)
n_peaks     <- nrow(results$peaks_genomic)
n_chr       <- uniqueN(results$snp_coverage$chrom)
aneuploid   <- ploidy_map[estimated_ploidy != 2L]
cat("\nResults:\n")
cat("  Chromosomes analyzed :", n_chr, "\n")
if (nrow(aneuploid) > 0) {
  cat("  Ploidy anomalies     :\n")
  for (i in seq_len(nrow(aneuploid))) {
    cat(sprintf("    %-20s %dN  (depth ratio %.2f)\n",
                as.character(aneuploid$chrom[i]),
                aneuploid$estimated_ploidy[i],
                aneuploid$depth_ratio[i]))
  }
} else {
  cat("  Ploidy               : all chromosomes 2N\n")
}
cat("  Chimeric reads       :", n_chimeric, "\n")
cat("  Peaks detected       :", n_peaks, "\n")
if (!is.null(results$snp_peaks) && "haplotype_label" %in% names(results$snp_peaks)) {
  label_counts <- results$snp_peaks[!is.na(snp_pos), .N, by = haplotype_label]
  setorder(label_counts, -N)
  cat("  Peak types           :\n")
  for (i in seq_len(nrow(label_counts))) {
    cat(sprintf("    %-20s %d\n", label_counts$haplotype_label[i], label_counts$N[i]))
  }
}
cat("\n")


# ── Write output ──────────────────────────────────────────────────────────────
if (output_mode == "csv") {

  # Peak list
  message("Writing peak list → ", out_path)
  data.table::fwrite(results$snp_peaks, out_path)
  cat("Peak list written to:", out_path, "\n")

} else if (output_mode == "rds") {

  # Overview plot RDS — same payload as the Shiny download_plot_rds handler
  message("Saving overview RDS → ", out_path)
  saveRDS(
    list(
      snp_coverage    = data.table::copy(results$snp_coverage),
      chromosome_fits = data.table::copy(results$chromosome_fits),
      peaks_genomic   = results$peaks_genomic,
      snp_peaks       = results$snp_peaks,
      sample_name     = opts[["sample-name"]],
      app_version     = APP_VERSION
    ),
    out_path
  )
  cat("Overview RDS saved to:", out_path, "\n")
  cat("Re-plot with:\n")
  cat('  library(ggplot2); source("R/chimera_functions.R")\n')
  cat('  d <- readRDS("', out_path, '")\n', sep = "")
  cat('  print(build_overview_plot(d))\n')

} else if (output_mode == "events") {

  # Events-only mode: the events table is produced by the chain pipeline below,
  # so defer writing until the chain block has built final_events.
  message("Events table will be written after chain analysis ...")

} else if (opts[["chain-all"]]) {

  # PNG overview plot, but --chain-all is also requested: defer the render.
  # The chain block below populates results$loh_segments and
  # results$event_table, which build_overview_plot() overlays as the LOH band
  # and event symbols. Building here would produce a peaks-only plot, so we
  # wait and render the annotated overview at the end of the chain block.
  message("Overview plot will be rendered after chain analysis (annotated) ...")

} else {

  # Default: PNG overview plot
  message("Building overview plot ...")
  p       <- build_overview_plot(results)
  n_chr   <- length(unique(results$snp_coverage$chrom))
  png_h   <- max(3, min(16, n_chr * 1.2))

  message("Saving PNG → ", out_path)
  ggplot2::ggsave(out_path, plot = p, width = 12, height = png_h, dpi = 300)
  cat("Overview plot saved to:", out_path, "\n")
}

# When output is produced by the chain block below (deferred annotated PNG for
# --chain-all, or the events table for --events-table), defer "Done." too so it
# prints after that output is actually written.
if (!(output_mode == "png" && opts[["chain-all"]]) && output_mode != "events")
  cat("Done.\n")


# ── Coverage map (sequencing depth) ─────────────────────────────────────────────
# Real per-position read depth, modeled with the same EM + Viterbi HMM approach
# as the LOH map but on total depth instead of allele balance
# (compute_coverage_map() in chimera_functions.R). Computed here whenever either
# --coverage-map or --chain-all is requested, since chain analysis's Step 0a
# needs the same result and would otherwise recompute it.
coverage_result <- NULL
if (opts[["coverage-map"]] || run_chain) {
  message("Computing coverage map ...")
  coverage_result <- compute_coverage_map(results$full_read_loh)
}

if (opts[["coverage-map"]]) {
  cov_dir  <- resolve_extra_dir(opts[["output"]], out_path)
  dir.create(cov_dir, recursive = TRUE, showWarnings = FALSE)
  cov_stem <- file.path(cov_dir, paste0(opts[["sample-name"]], "_", Sys.Date()))
  cov_table_path <- paste0(cov_stem, "_coverage_table.csv")
  cov_segs_path  <- paste0(cov_stem, "_coverage_segments.csv")

  message("Writing coverage table → ", cov_table_path)
  data.table::fwrite(coverage_result$coverage_table, cov_table_path)
  message("Writing coverage segments → ", cov_segs_path)
  data.table::fwrite(coverage_result$coverage_segments, cov_segs_path)

  cat(sprintf("Coverage table     → %s (%d positions)\n",
              cov_table_path, nrow(coverage_result$coverage_table)))
  cat(sprintf("Coverage segments  → %s (%d segments; baseline depth = %.1f)\n",
              cov_segs_path, nrow(coverage_result$coverage_segments),
              coverage_result$baseline_depth))
}


# ── Chain-based LOH event calling (--chain-all / --events-table) ──────────────
# The full pipeline (steps 0–4) always runs; `full_chain` decides how much is
# written. Under --chain-all we emit every intermediate step CSV and the
# annotated overview PNG. Under --events-table we run the same pipeline but write
# only the final events table (to out_path).
if (run_chain) {

  full_chain <- opts[["chain-all"]]

  # ── Helpers: flatten chain structures to data.tables ─────────────────────────

  # One row per token per chromosome
  tokens_to_dt <- function(chains) {
    rows <- list()
    for (cname in names(chains)) {
      ch <- chains[[cname]]
      for (idx in seq_along(ch$tokens)) {
        tok <- ch$tokens[[idx]]
        rows[[length(rows) + 1L]] <- data.table::data.table(
          chrom          = cname,
          token_idx      = idx,
          type           = tok$type,
          state          = tok$state,
          start          = tok$start,
          end            = tok$end,
          length_bp      = tok$length_bp,
          n_snps         = tok$n_snps,
          depth_ratio    = tok$depth_ratio,
          bridged_gap    = tok$bridged_gap,
          peak_over_pos  = if (!is.null(tok$peak_over))
                             tok$peak_over$fused_pos_bp %||% tok$peak_over$snp_pos
                           else NA_integer_,
          peak_left_pos  = if (!is.null(tok$peak_left))
                             tok$peak_left$fused_pos_bp  %||% tok$peak_left$snp_pos
                           else NA_integer_,
          peak_right_pos = if (!is.null(tok$peak_right))
                             tok$peak_right$fused_pos_bp %||% tok$peak_right$snp_pos
                           else NA_integer_
        )
      }
    }
    data.table::rbindlist(rows, fill = TRUE)
  }

  # One row per event — used for Step 3 (pre-reconcile) events
  scan_events_to_dt <- function(scan_results) {
    rows <- list()
    for (cname in names(scan_results)) {
      for (ev in scan_results[[cname]]$events) {
        rows[[length(rows) + 1L]] <- data.table::data.table(
          event_class     = ev$event_class,
          chrom           = ev$chrom,
          start           = as.integer(ev$start),
          end             = as.integer(ev$end),
          length_bp       = as.integer(ev$length_bp),
          n_support       = as.integer(ev$n_support),
          peak_edge_types = ev$peak_edge_types %||% NA_character_,
          notes           = ev$notes %||% ""
        )
      }
    }
    if (length(rows) == 0)
      return(data.table::data.table(
        event_class=character(), chrom=character(),
        start=integer(), end=integer(), length_bp=integer(),
        n_support=integer(), peak_edge_types=character(), notes=character()
      ))
    data.table::rbindlist(rows, fill = TRUE)
  }

  # Resolve the output directory — place chain CSVs alongside the main output
  out_dir <- resolve_extra_dir(opts[["output"]], out_path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(out_dir,
                    paste0(opts[["sample-name"]], "_", Sys.Date(), "_chain"))

  # Build chain params from CLI opts (remaining tuning uses defaults)
  cp <- default_chain_params()
  cp$tel_tol_bp        <- as.integer(opts[["tel-tol"]]   * 1000L)
  cp$merge_gap_bp      <- as.integer(opts[["merge-gap"]] * 1000L)
  cp$min_span          <- as.integer(opts[["min-span"]])
  cp$peak_pad_bp       <- as.integer(opts[["peak-pad"]])
  cp$min_snps_for_peak <- as.integer(opts[["min-run"]])

  cat("\n── Chain analysis ───────────────────────────────────────────────────────────\n")

  # ── Step 0: LOH map ───────────────────────────────────────────────────────────
  cat("[chain] Step 0: Computing LOH map ...\n")
  loh_result <- compute_loh_map(results$full_read_loh)
  loh_segs   <- loh_result$loh_segments

  if (full_chain) {
    step0_snp <- paste0(stem, "_step0_loh_snps.csv")
    step0_seg <- paste0(stem, "_step0_loh_segments.csv")
    data.table::fwrite(loh_result$snp_table, step0_snp)
    data.table::fwrite(loh_segs,             step0_seg)
    cat(sprintf("  LOH SNP table    → %s (%d rows)\n",  step0_snp, nrow(loh_result$snp_table)))
    cat(sprintf("  LOH segments     → %s (%d segments)\n\n", step0_seg, nrow(loh_segs)))
  }

  # ── Step 0a: Coverage map ────────────────────────────────────────────────────
  # Real per-position read depth, modeled with the same EM+HMM approach as
  # the LOH map above but on total depth instead of allele balance. Used by
  # R01/R02b to tell a true terminal/arm deletion (real depth drop) apart
  # from a terminal LOH/crossover that is simply missing its junction peak
  # (depth consistent with the rest of the chromosome). Already computed above
  # if --coverage-map was also requested; otherwise computed here on demand.
  cat("[chain] Step 0a: Computing coverage map ...\n")
  if (is.null(coverage_result)) coverage_result <- compute_coverage_map(results$full_read_loh)
  coverage_segs    <- coverage_result$coverage_segments

  if (full_chain) {
    step0a_cov <- paste0(stem, "_step0a_coverage_segments.csv")
    data.table::fwrite(coverage_segs, step0a_cov)
    cat(sprintf("  Coverage segments → %s (%d segments; baseline depth = %.1f)\n\n",
                step0a_cov, nrow(coverage_segs), coverage_result$baseline_depth))
  }

  # ── Step 0b: Haplotype labeling + peak fusion ────────────────────────────────
  # compute_peak_pairs() now runs classify_peak_haplotype() internally for any
  # unlabeled peak before evaluating pairs, so no separate labeling step is
  # needed here. Labels and labeled snp_peaks are returned in fusion_res.
  cat("[chain] Step 0b: Labeling peak haplotypes + computing peak pairs ...\n")
  fusion_res <- compute_peak_pairs(
    snp_peaks      = results$snp_peaks,
    rt_df          = results$rt_df,
    transition_pos = results$transition_pos,
    loh_segments   = loh_segs,
    zone_min_snps  = as.integer(opts[["min-run"]]),
    homog_frac     = cp$homog_frac
  )
  # Update snp_peaks with the labels compute_peak_pairs assigned so downstream
  # steps (chain build, reconcile, CSV export) see them.
  if (!is.null(fusion_res$snp_peaks)) results$snp_peaks <- fusion_res$snp_peaks

  n_labeled <- if (!is.null(results$snp_peaks))
    sum(!is.na(results$snp_peaks$haplotype_label)) else 0L
  n_pairs   <- if (!is.null(fusion_res$peak_pairs)) nrow(fusion_res$peak_pairs) else 0L
  cat(sprintf("  Labeled %d / %d peaks\n", n_labeled,
              if (!is.null(results$snp_peaks)) nrow(results$snp_peaks) else 0L))
  cat(sprintf("  Peak pairs: %d candidate pairs evaluated\n\n", n_pairs))

  # ── Step 1: Build raw chains ──────────────────────────────────────────────────
  cat("[chain] Step 1: Building raw chains ...\n")
  raw_chains <- build_raw_chains(
    loh_segments      = loh_segs,
    chr_span          = results$chr_span,
    params            = cp,
    fused_peaks       = fusion_res$fused_peaks,
    peak_pairs        = fusion_res$peak_pairs,
    snp_peaks         = results$snp_peaks,
    coverage_segments = coverage_segs,
    coverage_table    = coverage_result$coverage_table
  )

  step1_dt <- tokens_to_dt(raw_chains)
  if (full_chain) {
    step1_out <- paste0(stem, "_step1_raw_tokens.csv")
    data.table::fwrite(step1_dt, step1_out)
    cat(sprintf("  Raw tokens       → %s (%d tokens across %d chromosomes)\n\n",
                step1_out, nrow(step1_dt), length(raw_chains)))
  }

  # ── Step 2: Canonicalise ──────────────────────────────────────────────────────
  cat("[chain] Step 2: Canonicalising (merging same-state gaps) ...\n")
  canonical_chains <- lapply(raw_chains, canonicalise, params = cp)

  if (full_chain) {
    step2_out <- paste0(stem, "_step2_canonical_tokens.csv")
    step2_dt  <- tokens_to_dt(canonical_chains)
    data.table::fwrite(step2_dt, step2_out)
    cat(sprintf("  Canonical tokens → %s (%d tokens; %d merged from step 1)\n\n",
                step2_out, nrow(step2_dt), nrow(step1_dt) - nrow(step2_dt)))
  }

  # ── Step 3: Motif scan ────────────────────────────────────────────────────────
  cat("[chain] Step 3: Scanning for recombination motifs ...\n")
  scan_results <- lapply(names(canonical_chains), function(cname) {
    scan_chain(canonical_chains[[cname]], cp)
  })
  names(scan_results) <- names(canonical_chains)

  # Update chains with any token rewrites produced during scanning
  for (cname in names(canonical_chains)) {
    canonical_chains[[cname]] <- scan_results[[cname]]$chain
  }

  if (full_chain) {
    step3_ev  <- paste0(stem, "_step3_events.csv")
    step3_tok <- paste0(stem, "_step3_tokens_post_scan.csv")
    step3_dt  <- scan_events_to_dt(scan_results)
    step3_tok_dt <- tokens_to_dt(canonical_chains)
    data.table::fwrite(step3_dt,     step3_ev)
    data.table::fwrite(step3_tok_dt, step3_tok)
    cat(sprintf("  Pre-reconcile events → %s (%d events)\n",    step3_ev,  nrow(step3_dt)))
    cat(sprintf("  Post-scan tokens     → %s (%d tokens)\n\n",  step3_tok, nrow(step3_tok_dt)))
  }

  # ── Step 4: Reconcile ─────────────────────────────────────────────────────────
  cat("[chain] Step 4: Reconciling other events (unresolved LOH/peaks) ...\n")
  rec <- reconcile(
    scan_results = scan_results,
    chains       = canonical_chains,
    fused_peaks  = fusion_res$fused_peaks,
    peak_pairs   = fusion_res$peak_pairs,
    snp_peaks    = results$snp_peaks
  )

  # Under --chain-all the events table is one of many step CSVs (stem-named);
  # under --events-table it is the single requested output, written to out_path.
  step4_ev <- if (full_chain) paste0(stem, "_step4_final_events.csv") else out_path
  final_events <- build_event_table(rec$events)
  data.table::fwrite(final_events, step4_ev)
  cat(sprintf("  Final events     → %s (%d total; %d high, %d review)\n",
              step4_ev, nrow(final_events),
              sum(final_events$confidence == "high"),
              sum(final_events$confidence == "review")))

  if (full_chain && length(rec$unclaimed_loh) > 0) {
    step4_ul <- paste0(stem, "_step4_other_events_loh.csv")
    uncl_loh <- data.table::rbindlist(lapply(rec$unclaimed_loh, function(u)
      data.table::data.table(
        chrom     = u$chrom,
        start     = u$start,
        end       = u$end,
        length_bp = u$end - u$start,
        state     = u$state,
        n_snps    = u$n_snps,
        reason    = u$reason
      )))
    data.table::fwrite(uncl_loh, step4_ul)
    cat(sprintf("  Other events (LOH)   → %s (%d segments)\n", step4_ul, nrow(uncl_loh)))
  }

  if (full_chain && length(rec$unclaimed_peaks) > 0) {
    step4_up <- paste0(stem, "_step4_other_events_peaks.csv")
    uncl_pk <- data.table::rbindlist(lapply(rec$unclaimed_peaks, function(u)
      data.table::data.table(
        chrom     = u$chrom,
        snp_pos   = u$snp_pos,
        edge_type = u$edge_type,
        reason    = u$reason
      )))
    data.table::fwrite(uncl_pk, step4_up)
    cat(sprintf("  Other events (peaks) → %s (%d peaks)\n", step4_up, nrow(uncl_pk)))
  }

  # ── Annotated overview PNG ────────────────────────────────────────────────────
  # In the default (PNG) output mode the overview render was deferred above so
  # it could be annotated. Feed the chain results back into `results` — the LOH
  # segment table and the final event table — so build_overview_plot() draws the
  # LOH band and event symbols, matching the app's annotated overview. Skipped
  # for --peak-list / --overview-rds, where out_path is a CSV/RDS file.
  if (output_mode == "png") {
    results$loh_segments <- loh_segs
    results$event_table  <- final_events

    message("Building annotated overview plot ...")
    p     <- build_overview_plot(results)
    n_chr <- length(unique(results$snp_coverage$chrom))
    png_h <- max(3, min(16, n_chr * 1.2))

    message("Saving PNG → ", out_path)
    ggplot2::ggsave(out_path, plot = p, width = 12, height = png_h, dpi = 300)
    cat(sprintf("Annotated overview plot (LOH + events) saved to: %s\n", out_path))
  }

  cat("\n── Chain analysis complete ───────────────────────────────────────────────────\n")

  # Deferred "Done." — the non-chain path prints this right after writing the
  # main output; here the PNG (deferred --chain-all) or the events table
  # (--events-table) is written above, so print it now.
  if (output_mode %in% c("png", "events")) cat("Done.\n")
}
