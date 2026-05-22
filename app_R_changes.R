# =============================================================================
#  Changes required in app.R to use chimera_functions.R
#
#  Three edits — all small.
# =============================================================================


# ── EDIT 1: top of file ───────────────────────────────────────────────────────
# Replace the inline whittaker() definition and library calls with a source().
# Keep the shiny library; remove data.table/pracma/ggplot2/Matrix (now in functions file).

# REMOVE (lines 1-23 of current app.R):
library(shiny)
library(data.table)
library(pracma)
library(ggplot2)
library(Matrix)

whittaker <- function(y, lambda = 1, d = 2) { ... }   # entire function

APP_VERSION <- "0.5.0"

# REPLACE WITH:
library(shiny)
source("R/chimera_functions.R")   # loads all packages + whittaker + run_chimera_analysis etc.
# APP_VERSION is now defined inside chimera_functions.R


# ── EDIT 2: inside observeEvent(input$run_analysis, { ... }) ─────────────────
# Replace the entire analysis body (steps 1–13, roughly lines 251–726)
# with a single call to run_chimera_analysis(), then unpack results.

# REMOVE (everything inside withProgress({ ... }) except the incProgress calls):
#   Steps 1–11 (load data → build chromosome_fits)

# REPLACE WITH:
withProgress(message = "Processing data...", value = 0, {

  incProgress(0.1, detail = "Running analysis")

  res <- run_chimera_analysis(
    read_data_path  = input$read_data_file$datapath,
    snp_data_path   = input$snp_data_file$datapath,
    chr_size_path   = input$chr_size_file$datapath,
    sample_name     = input$sample_name,
    mapq_cutoff     = input$mapq_cutoff,
    baseq_cutoff    = input$baseq_cutoff,
    min_run         = input$min_run,
    min_peak_height = input$min_peak_height,
    lambda          = input$lambda,
    # Redirect warnings to Shiny notifications instead of console messages
    warn_fn = function(msg) showNotification(msg, type = "warning", duration = 10)
  )

  incProgress(0.85, detail = "Storing results")

  results$rt_df             <- res$rt_df
  results$chimeric_read_ids <- res$chimeric_read_ids
  results$transition_pos    <- res$transition_pos
  results$snp_coverage      <- res$snp_coverage
  results$peaks_genomic     <- res$peaks_genomic
  results$snp_peaks         <- res$snp_peaks
  results$chromosome_fits   <- res$chromosome_fits
  results$chr_span          <- res$chr_span

  # Peak plots are still built here — they depend on Shiny output IDs and
  # are not worth moving into chimera_functions.R.
  # Keep the existing lapply(peak_chrs, ...) block unchanged.
  incProgress(0.9, detail = "Creating individual peak plots")
  # ... (existing peak_plots_list code unchanged) ...

  incProgress(1, detail = "Complete")
})


# ── EDIT 3: output$chr_plot renderPlot and output$download_plot ──────────────
# Both now call build_overview_plot() instead of duplicating the ggplot code.

# output$chr_plot — REPLACE the ggplot block with:
output$chr_plot <- renderPlot({
  req(results$snp_coverage, results$chromosome_fits)
  build_overview_plot(results)
}, height = function() {
  req(results$snp_coverage)
  n_chr <- length(unique(results$snp_coverage$chrom))
  min(1600, max(400, n_chr * 120))
})

# output$download_plot content function — REPLACE the ggplot block with:
output$download_plot <- downloadHandler(
  filename = function() paste0(input$sample_name, "_chromosome_tracking_", Sys.Date(), ".png"),
  content  = function(file) {
    req(results$snp_coverage, results$chromosome_fits)
    p     <- build_overview_plot(results)
    n_chr <- length(unique(results$snp_coverage$chrom))
    png_h <- max(3, min(16, n_chr * 1.2))
    ggsave(file, plot = p, width = 12, height = png_h, dpi = 300)
  }
)
