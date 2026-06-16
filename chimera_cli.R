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

# Validate output mode
if (opts[["peak-list"]] && opts[["overview-rds"]])
  stop("--peak-list and --overview-rds are mutually exclusive. Choose one.")

if (opts[["chain-all"]] && is.null(chain_found))
  stop("--chain-all requires loh_chain_analysis.R but it was not found next to chimera_functions.R")


# ── Resolve output path ───────────────────────────────────────────────────────
resolve_output <- function(opts_output, sample_name, mode) {
  ext <- switch(mode,
    png = ".png",
    csv = "_peaks.csv",
    rds = "_overview.rds"
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

output_mode <- if (opts[["peak-list"]])    "csv" else
               if (opts[["overview-rds"]]) "rds" else
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
cat("  Min run         :", opts[["min-run"]],         "\n")
cat("  Min peak height :", opts[["min-peak-height"]], "\n")
cat("  Lambda (λ)      :", opts[["lambda"]],          "\n")
cat("Output mode       :", output_mode, "→", out_path, "\n")
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
  min_run         = opts[["min-run"]],
  min_peak_height = opts[["min-peak-height"]],
  lambda          = opts[["lambda"]]
)


# ── Print run summary ─────────────────────────────────────────────────────────
n_chimeric  <- length(results$chimeric_read_ids)
n_peaks     <- nrow(results$peaks_genomic)
n_chr       <- uniqueN(results$snp_coverage$chrom)
cat("\nResults:\n")
cat("  Chromosomes analyzed :", n_chr, "\n")
cat("  Chimeric reads       :", n_chimeric, "\n")
cat("  Peaks detected       :", n_peaks, "\n\n")


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

cat("Done.\n")


# ── Chain-based LOH event calling (--chain-all) ───────────────────────────────
if (opts[["chain-all"]]) {

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
  out_dir <- if (!is.null(opts[["output"]]) &&
                 (dir.exists(opts[["output"]]) || grepl("/$", opts[["output"]]))) {
    opts[["output"]]
  } else {
    dirname(normalizePath(out_path, mustWork = FALSE))
  }
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

  step0_snp <- paste0(stem, "_step0_loh_snps.csv")
  step0_seg <- paste0(stem, "_step0_loh_segments.csv")
  data.table::fwrite(loh_result$snp_table, step0_snp)
  data.table::fwrite(loh_segs,             step0_seg)
  cat(sprintf("  LOH SNP table    → %s (%d rows)\n",  step0_snp, nrow(loh_result$snp_table)))
  cat(sprintf("  LOH segments     → %s (%d segments)\n\n", step0_seg, nrow(loh_segs)))

  # ── Step 0b: Haplotype labeling ───────────────────────────────────────────────
  # classify_peak_haplotype() is defined in chimera_functions.R.
  # The chain rules require haplotype_label on each snp_peaks row so that
  # best_edge_type can be resolved; without this all events come out AMBIGUOUS.
  cat("[chain] Step 0b: Labeling peak haplotypes ...\n")
  sp <- results$snp_peaks
  if (!is.null(sp) && nrow(sp) > 0 && !is.null(results$transition_pos)) {
    sp[, .row_idx := .I]
    hap_rows <- list()
    for (chr_name in unique(as.character(sp$chrom))) {
      chr_peaks <- sp[as.character(chrom) == chr_name & !is.na(snp_pos)][order(snp_pos)]
      for (pk_i in seq_len(nrow(chr_peaks))) {
        pk <- chr_peaks[pk_i]
        touching_ids <- results$transition_pos[
          as.character(chrom) == chr_name &
            pos >= pk$peak_start & pos <= pk$peak_end,
          unique(read_id)
        ]
        if (length(touching_ids) == 0) next
        hap <- classify_peak_haplotype(
          pk            = pk,
          chr_name      = chr_name,
          rt_df         = results$rt_df,
          touching_ids  = touching_ids,
          zone_min_snps = as.integer(opts[["min-run"]])
        )
        hap_rows[[length(hap_rows) + 1L]] <- data.table::data.table(
          .row_idx        = pk$.row_idx,
          haplotype_label = hap$label
        )
      }
    }
    if (length(hap_rows) > 0) {
      hap_dt <- data.table::rbindlist(hap_rows)
      sp <- merge(sp, hap_dt, by = ".row_idx", all.x = TRUE)
    }
    sp[, .row_idx := NULL]
    results$snp_peaks <- sp
    n_labeled <- if (length(hap_rows) > 0) nrow(data.table::rbindlist(hap_rows)) else 0L
    cat(sprintf("  Labeled %d / %d peaks\n\n", n_labeled, nrow(sp)))
  }

  # ── Step 0c: Peak fusion ─────────────────────────────────────────────────────
  cat("[chain] Step 0c: Computing peak pairs (fusion analysis) ...\n")
  fusion_res <- compute_peak_pairs(
    snp_peaks      = results$snp_peaks,
    rt_df          = results$rt_df,
    transition_pos = results$transition_pos,
    loh_segments   = loh_segs,
    zone_min_snps  = as.integer(opts[["min-run"]])
  )
  n_pairs <- if (!is.null(fusion_res$peak_pairs)) nrow(fusion_res$peak_pairs) else 0L
  cat(sprintf("  Peak pairs: %d candidate pairs evaluated\n\n", n_pairs))

  # ── Step 1: Build raw chains ──────────────────────────────────────────────────
  cat("[chain] Step 1: Building raw chains ...\n")
  raw_chains <- build_raw_chains(
    loh_segments = loh_segs,
    chr_span     = results$chr_span,
    params       = cp,
    fused_peaks  = fusion_res$fused_peaks,
    peak_pairs   = fusion_res$peak_pairs,
    snp_peaks    = results$snp_peaks
  )

  step1_out  <- paste0(stem, "_step1_raw_tokens.csv")
  step1_dt   <- tokens_to_dt(raw_chains)
  data.table::fwrite(step1_dt, step1_out)
  cat(sprintf("  Raw tokens       → %s (%d tokens across %d chromosomes)\n\n",
              step1_out, nrow(step1_dt), length(raw_chains)))

  # ── Step 2: Canonicalise ──────────────────────────────────────────────────────
  cat("[chain] Step 2: Canonicalising (merging same-state gaps) ...\n")
  canonical_chains <- lapply(raw_chains, canonicalise, params = cp)

  step2_out <- paste0(stem, "_step2_canonical_tokens.csv")
  step2_dt  <- tokens_to_dt(canonical_chains)
  data.table::fwrite(step2_dt, step2_out)
  cat(sprintf("  Canonical tokens → %s (%d tokens; %d merged from step 1)\n\n",
              step2_out, nrow(step2_dt), nrow(step1_dt) - nrow(step2_dt)))

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

  step3_ev  <- paste0(stem, "_step3_events.csv")
  step3_tok <- paste0(stem, "_step3_tokens_post_scan.csv")
  step3_dt  <- scan_events_to_dt(scan_results)
  step3_tok_dt <- tokens_to_dt(canonical_chains)
  data.table::fwrite(step3_dt,     step3_ev)
  data.table::fwrite(step3_tok_dt, step3_tok)
  cat(sprintf("  Pre-reconcile events → %s (%d events)\n",    step3_ev,  nrow(step3_dt)))
  cat(sprintf("  Post-scan tokens     → %s (%d tokens)\n\n",  step3_tok, nrow(step3_tok_dt)))

  # ── Step 4: Reconcile ─────────────────────────────────────────────────────────
  cat("[chain] Step 4: Reconciling unclaimed tokens and peaks ...\n")
  rec <- reconcile(
    scan_results = scan_results,
    chains       = canonical_chains,
    fused_peaks  = fusion_res$fused_peaks,
    peak_pairs   = fusion_res$peak_pairs,
    snp_peaks    = results$snp_peaks
  )

  step4_ev <- paste0(stem, "_step4_final_events.csv")
  final_events <- build_event_table(rec$events)
  data.table::fwrite(final_events, step4_ev)
  cat(sprintf("  Final events     → %s (%d total; %d high, %d review)\n",
              step4_ev, nrow(final_events),
              sum(final_events$confidence == "high"),
              sum(final_events$confidence == "review")))

  if (length(rec$unclaimed_loh) > 0) {
    step4_ul <- paste0(stem, "_step4_unclaimed_loh.csv")
    uncl_loh <- data.table::rbindlist(lapply(rec$unclaimed_loh, function(u)
      data.table::data.table(
        chrom     = u$chrom,
        start     = u$start,
        end       = u$end,
        length_bp = u$end - u$start,
        state     = u$state,
        n_snps    = u$n_snps
      )))
    data.table::fwrite(uncl_loh, step4_ul)
    cat(sprintf("  Unclaimed LOH    → %s (%d segments)\n", step4_ul, nrow(uncl_loh)))
  }

  if (length(rec$unclaimed_peaks) > 0) {
    step4_up <- paste0(stem, "_step4_unclaimed_peaks.csv")
    uncl_pk <- data.table::rbindlist(lapply(rec$unclaimed_peaks, function(u)
      data.table::data.table(
        chrom     = u$chrom,
        snp_pos   = u$snp_pos,
        edge_type = u$edge_type
      )))
    data.table::fwrite(uncl_pk, step4_up)
    cat(sprintf("  Unclaimed peaks  → %s (%d peaks)\n", step4_up, nrow(uncl_pk)))
  }

  cat("\n── Chain analysis complete ───────────────────────────────────────────────────\n")
}
