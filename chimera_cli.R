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
#      --min-peak-height 15 \
#      --lambda 0.5 \
#      --output results/ \
#      reads.csv.gz snps.vcf.gz genome.fa.fai
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# ── Source shared analysis functions ─────────────────────────────────────────
# Looks for chimera_functions.R in R/ relative to this script's location,
# then falls back to the current working directory.
script_dir   <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),
                          error = function(e) ".")
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
              default = 1.0,
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
cat("Output mode       :", output_mode, "→", out_path, "\n\n")

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
