# =============================================================================
#  app_chain_patch.R
#
#  This file documents the MINIMAL changes needed to add chain-based LOH event
#  calling to the existing app.R.  It is NOT a standalone file — apply the
#  changes described below to app.R in the locations indicated.
#
#  No existing code is removed or modified; all additions are purely additive.
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
#  CHANGE 1 — Source the new module (add after the existing source() call)
#  Location: near line 6, after:  source("chimera_functions.R")
# ─────────────────────────────────────────────────────────────────────────────

source("loh_chain_analysis.R")


# ─────────────────────────────────────────────────────────────────────────────
#  CHANGE 2 — Add chain parameters to the sidebar
#  Location: in sidebarPanel(), after the existing Peak Fusion Parameters block
#  (after the jaccard_threshold numericInput, before the actionButton block)
# ─────────────────────────────────────────────────────────────────────────────

hr(),
h3("Event Caller Parameters"),

numericInput("tel_tol_kb",
             "Telomere Tolerance (Kb):",
             value = 5, min = 0, step = 1),
helpText("LOH regions within this distance of a chromosome end are treated as terminal."),

numericInput("merge_gap_kb",
             "Same-State Merge Gap (Kb):",
             value = 5, min = 0, step = 1),
helpText("NA gaps smaller than this between two same-state LOH runs will be merged."),

numericInput("min_span_reads",
             "Min Spanning Reads for Read Call:",
             value = 3, min = 1, step = 1),
helpText("Minimum spanning reads required to make a read-based event call. Below this, result is Ambiguous."),

numericInput("homog_frac",
             "NCO-GC Homogeneity Fraction:",
             value = 0.80, min = 0.5, max = 1.0, step = 0.05),
helpText("Fraction of spanning reads that must share the return pattern (AAABBBAAA) to call NCO-GC."),

actionButton("run_chain",
             "Run Event Caller",
             class = "btn-success btn-lg"),
helpText("Run after Peak Fusion. Calls recombination events from LOH + peak chain."),


# ─────────────────────────────────────────────────────────────────────────────
#  CHANGE 3 — Add chain results to reactiveValues
#  Location: inside results <- reactiveValues(...), after the existing fields
# ─────────────────────────────────────────────────────────────────────────────

chain_result   = NULL,   # full output of run_chain_analysis()
event_table    = NULL,   # flat data.table of called events
chain_params   = NULL,   # params list used for the last chain run


# ─────────────────────────────────────────────────────────────────────────────
#  CHANGE 4 — Add chain analysis observer
#  Location: in server(), after the existing observeEvent(input$run_fusion, {...})
#  block (around line 1751)
# ─────────────────────────────────────────────────────────────────────────────

observeEvent(input$run_chain, {
  req(results$loh_segments, results$snp_peaks, results$chr_span)

  # Build params from UI inputs
  cp <- default_chain_params()
  cp$tel_tol_bp   <- as.integer(round(input$tel_tol_kb   * 1000))
  cp$merge_gap_bp <- as.integer(round(input$merge_gap_kb * 1000))
  cp$min_span     <- as.integer(input$min_span_reads)
  cp$homog_frac   <- input$homog_frac

  results$chain_params <- cp

  withProgress(message = "Running chain event caller...", value = 0.2, {

    chain_res <- run_chain_analysis(
      loh_segments = results$loh_segments,
      fused_peaks  = results$fused_peaks,   # NULL if fusion not yet run — ok
      peak_pairs   = results$peak_pairs,
      snp_peaks    = results$snp_peaks,     # fallback when fused_peaks is NULL
      rt_df        = results$rt_df,
      chr_span     = results$chr_span,
      params       = cp
    )

    incProgress(0.8, detail = "Storing event table")
    results$chain_result <- chain_res
    results$event_table  <- chain_res$event_table
  })

  n_high   <- sum(results$event_table$confidence == "high",   na.rm = TRUE)
  n_review <- sum(results$event_table$confidence == "review", na.rm = TRUE)
  showNotification(
    paste0("Event calling complete. ",
           n_high, " high-confidence events, ",
           n_review, " requiring review."),
    type = "message", duration = 6
  )
})


# ─────────────────────────────────────────────────────────────────────────────
#  CHANGE 5 — Add Event Table tab to the UI tabsetPanel
#  Location: in mainPanel tabsetPanel(), after the "LOH Regions" tabPanel
# ─────────────────────────────────────────────────────────────────────────────

tabPanel("Recombination Events",
  h4("Called Recombination Events"),
  helpText(
    "Events called by the chain-based LOH + peak motif scanner. ",
    "Run Analysis, then Run Peak Fusion, then Run Event Caller. ",
    "High-confidence events have unambiguous read and LOH support. ",
    "'Review' events need manual inspection."
  ),
  br(),
  uiOutput("event_table_legend"),
  br(),
  tableOutput("event_table_output"),
  br(),
  fluidRow(
    column(3, downloadButton("download_event_table",
                             "Download Event Table (.csv)")),
    column(3, downloadButton("download_chain_rds",
                             "Download Chain Objects (.rds)"))
  ),
  br(),
  h5("Unclaimed LOH / Peaks"),
  helpText("LOH regions and peaks not assigned to any event. These may need parameter adjustment or manual review."),
  tableOutput("unclaimed_table")
),


# ─────────────────────────────────────────────────────────────────────────────
#  CHANGE 6 — Server outputs for the new tab
#  Location: in server(), after the existing loh_regions_table output block
# ─────────────────────────────────────────────────────────────────────────────

# Event table legend
output$event_table_legend <- renderUI({
  tags$div(
    style = "font-size:0.85em; color:#444; margin-bottom:8px;",
    strong("Event classes: "),
    tags$span("NCO_GC",             style = "background:#d4edda; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("NCO_GC_LARGE",       style = "background:#c3e6cb; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("NCO_GC_subres",      style = "background:#b8daff; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("CO_GC",              style = "background:#ffeeba; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("CROSSOVER_NO_TRACT", style = "background:#ffe8a1; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("TERMINAL_LOH",       style = "background:#d6d8db; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("TERMINAL_DELETION",  style = "background:#f5c6cb; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("TCO_CAPTURED_TCO",   style = "background:#e2d9f3; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("DOUBLE_GC",          style = "background:#bee5eb; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("AMBIGUOUS",          style = "background:#fff3cd; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("UNCATEGORIZED",      style = "background:#f8d7da; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    br(), br(),
    strong("Confidence: "),
    tags$span("high",   style = "background:#d4edda; padding:2px 6px; border-radius:3px; margin-right:4px;"),
    tags$span("review", style = "background:#fff3cd; padding:2px 6px; border-radius:3px; margin-right:4px;")
  )
})

# Main event table
output$event_table_output <- renderTable({
  req(results$event_table)
  et <- copy(results$event_table)

  # Resolve strain names in event notes / display
  s_ref <- if (!is.null(results$strain_ref) && nzchar(results$strain_ref))
    results$strain_ref else "REF"
  s_alt <- if (!is.null(results$strain_alt) && nzchar(results$strain_alt))
    results$strain_alt else "ALT"

  et[, `Start (bp)`  := format(start, big.mark = ",", scientific = FALSE)]
  et[, `End (bp)`    := format(end,   big.mark = ",", scientific = FALSE)]
  et[, `Length (kb)` := round(length_kb, 2)]
  et[, start := NULL][, end := NULL][, length_kb := NULL]

  setcolorder(et, c("event_class", "chrom", "Start (bp)", "End (bp)",
                    "Length (kb)", "n_support", "peak_edge_types",
                    "confidence", "notes"))
  setnames(et,
    c("event_class", "chrom", "n_support", "peak_edge_types", "confidence"),
    c("Event",       "Chr",   "N Support", "Peak Evidence",   "Confidence")
  )
  et
}, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s", na = "\u2014")

# Unclaimed LOH + peaks summary
output$unclaimed_table <- renderTable({
  req(results$chain_result)
  cr <- results$chain_result

  uncl_loh <- if (length(cr$unclaimed_loh) > 0) {
    rbindlist(lapply(cr$unclaimed_loh, function(u) data.table(
      Type  = "LOH",
      Chr   = u$chrom,
      Start = format(u$start, big.mark = ","),
      End   = format(u$end,   big.mark = ","),
      `Length (kb)` = round((u$end - u$start) / 1000, 2),
      Detail = paste0(u$state, "; snps=", u$n_snps)
    )), fill = TRUE)
  } else NULL

  uncl_pk <- if (length(cr$unclaimed_peaks) > 0) {
    rbindlist(lapply(cr$unclaimed_peaks, function(u) data.table(
      Type  = "Peak",
      Chr   = u$chrom,
      Start = format(as.integer(u$snp_pos), big.mark = ","),
      End   = format(as.integer(u$snp_pos), big.mark = ","),
      `Length (kb)` = 0,
      Detail = paste0("edge=", u$edge_type %||% "?")
    )), fill = TRUE)
  } else NULL

  out <- rbindlist(list(uncl_loh, uncl_pk), fill = TRUE)
  if (nrow(out) == 0)
    return(data.frame(Message = "All LOH segments and peaks accounted for."))
  setorder(out, Chr, Start)
  out
}, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s")

# Download event table
output$download_event_table <- downloadHandler(
  filename = function()
    paste0(input$sample_name, "_recombination_events_", Sys.Date(), ".csv"),
  content = function(file) {
    req(results$event_table)
    fwrite(results$event_table, file)
  }
)

# Download chain objects (RDS — for downstream analysis / ploidy layer)
output$download_chain_rds <- downloadHandler(
  filename = function()
    paste0(input$sample_name, "_loh_chains_", Sys.Date(), ".rds"),
  content = function(file) {
    req(results$chain_result)
    saveRDS(
      list(
        chains       = results$chain_result$chains,
        events       = results$chain_result$events,
        event_table  = results$event_table,
        chain_params = results$chain_params,
        sample_name  = input$sample_name,
        schema_version = CHAIN_SCHEMA_VERSION
      ),
      file
    )
  }
)
