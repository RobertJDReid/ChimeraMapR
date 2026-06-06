library(shiny)
if (!requireNamespace("igraph", quietly = TRUE))
  stop("Package 'igraph' is required. Install with: install.packages('igraph')")
library(igraph)

source("chimera_functions.R")   # loads all packages + whittaker + run_chimera_analysis etc.
# APP_VERSION is now defined inside chimera_functions.R

# ─────────────────────────────────────────────────────────────────────────────
#   PEAK FUSION FUNCTIONS
#   All fusion logic is self-contained here so it can be sourced / tested
#   independently of the Shiny server.
# ─────────────────────────────────────────────────────────────────────────────

# ---------------------------------------------------------------------------
# classify_zone_state()
#   For a single read's RLE-smoothed IS_REF vector and positions, return the
#   majority-vote haplotype ("REF" / "ALT") for a defined genomic interval.
#   Returns NA if fewer than min_snps positions fall in the interval.
# ---------------------------------------------------------------------------
classify_zone_state <- function(pos_vec, is_ref_vec, zone_start, zone_end, min_snps) {
  idx <- pos_vec >= zone_start & pos_vec <= zone_end
  if (sum(idx) < min_snps) return(NA_character_)
  majority <- mean(is_ref_vec[idx], na.rm = TRUE)
  if (is.na(majority)) return(NA_character_)
  if (majority >= 0.5) "REF" else "ALT"
}

# ---------------------------------------------------------------------------
# classify_edge_type()
#   Given the set of spanning reads (reads with data in all three zones),
#   determine whether the pair represents a gene_conversion, crossover,
#   independent_events, or ambiguous relationship.
#
#   state_df: data.frame with columns read_id, state_L, state_M, state_R
#             (each is "REF", "ALT", or NA)
#   Returns a character scalar.
# ---------------------------------------------------------------------------
classify_edge_type <- function(state_df) {

  # Drop reads where any zone is NA (insufficient SNPs for classification)
  classifiable <- state_df[!is.na(state_df$state_L) &
                           !is.na(state_df$state_R), ]

  # Two-zone fallback: if all middle zones are NA (no SNPs between peaks),
  # use L and R only to test for bimodality.
  no_middle <- all(is.na(state_df$state_M))

  if (nrow(classifiable) == 0) return("unresolvable")

  if (no_middle) {
    # Crossover with no LOH region: expect two complementary classes
    # e.g. REF-ALT and ALT-REF at the two flanks
    lr_patterns <- paste(classifiable$state_L, classifiable$state_R, sep = "-")
    unique_pats <- unique(lr_patterns)
    if (length(unique_pats) == 2) {
      # Check complementarity: REF-ALT paired with ALT-REF
      sorted <- sort(unique_pats)
      if (sorted[1] == "ALT-REF" && sorted[2] == "REF-ALT")
        return("crossover")
    }
    # Homogeneous: could be gene_conversion but no middle to confirm return
    if (length(unique_pats) == 1) return("ambiguous")
    return("ambiguous")
  }

  # Full three-zone classification
  classifiable_3 <- classifiable[!is.na(classifiable$state_M), ]

  if (nrow(classifiable_3) == 0) {
    # Fall back to two-zone
    lr_patterns <- paste(classifiable$state_L, classifiable$state_R, sep = "-")
    unique_pats <- unique(lr_patterns)
    if (length(unique_pats) == 2) {
      sorted <- sort(unique_pats)
      if (sorted[1] == "ALT-REF" && sorted[2] == "REF-ALT")
        return("crossover")
    }
    return("ambiguous")
  }

  patterns <- paste(classifiable_3$state_L,
                    classifiable_3$state_M,
                    classifiable_3$state_R, sep = "-")
  unique_pats <- unique(patterns)

  # Gene conversion: all reads show XYX (return to origin)
  is_return <- function(p) {
    parts <- strsplit(p, "-")[[1]]
    parts[1] == parts[3] && parts[1] != parts[2]
  }
  all_return <- all(sapply(unique_pats, is_return))
  if (all_return && length(unique_pats) <= 2) return("gene_conversion")

  # Crossover: two complementary classes, each switches at one peak only
  # REF-ALT-ALT paired with ALT-REF-REF  (or mirror)
  crossover_pairs <- list(
    c("REF-ALT-ALT", "ALT-REF-REF"),
    c("ALT-REF-REF", "REF-ALT-ALT")
  )
  if (length(unique_pats) == 2) {
    sorted <- sort(unique_pats)
    if ((sorted[1] == "ALT-REF-REF" && sorted[2] == "REF-ALT-ALT"))
      return("crossover")
  }

  # Independent events: no read returns to origin AND patterns are not
  # complementary crossover classes.  e.g. AAAAAAABBB + BBBAAAABBB
  no_return <- all(!sapply(unique_pats, is_return))
  if (no_return) return("independent_events")

  return("ambiguous")
}

# ---------------------------------------------------------------------------
# compute_peak_pairs()
#   Main fusion analysis function.  Takes the existing snp_peaks table and
#   rt_df (raw read data with RLE-smoothed IS_REF), returns:
#     - peak_pairs:  one row per candidate pair with Jaccard, edge type, etc.
#     - fused_peaks: original peaks table with fusion group membership added
# ---------------------------------------------------------------------------
compute_peak_pairs <- function(snp_peaks,
                               rt_df,
                               transition_pos,
                               jaccard_threshold  = 0.20,
                               zone_min_snps      = 2L,
                               supervised_override = NULL) {
  # supervised_override: character vector of "chrX_peakA_peakB" IDs that the
  # user has manually approved for fusion in the review UI.

  if (is.null(snp_peaks) || nrow(snp_peaks) == 0)
    return(list(peak_pairs = NULL, fused_peaks = NULL))

  peaks_dt <- copy(snp_peaks)
  peaks_dt[, peak_id := .I]                          # integer row ID
  peaks_dt[, chrom   := as.character(chrom)]

  all_pairs <- list()

  for (chr_name in unique(peaks_dt$chrom)) {

    chr_peaks <- peaks_dt[chrom == chr_name][order(snp_pos)]
    if (nrow(chr_peaks) < 2) next

    # Reads associated with each peak: transition boundary overlaps peak position
    # (using peak_start..peak_end window consistent with peak plot code)
    get_peak_reads <- function(pk) {
      transition_pos[
        as.character(chrom) == chr_name &
          pos >= pk$peak_start & pos <= pk$peak_end,
        unique(read_id)
      ]
    }

    # Pre-compute read sets for all peaks on this chromosome
    read_sets <- lapply(seq_len(nrow(chr_peaks)), function(i) get_peak_reads(chr_peaks[i]))
    names(read_sets) <- chr_peaks$peak_id

    # Evaluate adjacent pairs only
    for (j in seq_len(nrow(chr_peaks) - 1L)) {

      pk_a <- chr_peaks[j]
      pk_b <- chr_peaks[j + 1L]

      reads_a <- read_sets[[as.character(pk_a$peak_id)]]
      reads_b <- read_sets[[as.character(pk_b$peak_id)]]

      # Proximity threshold: median read length of reads associated with either peak
      all_pair_reads <- union(reads_a, reads_b)
      if (length(all_pair_reads) == 0) next

      read_lengths <- rt_df[read_id %in% all_pair_reads,
                            .(len = max(pos) - min(pos)), by = read_id]$len
      median_read_len <- median(read_lengths, na.rm = TRUE)
      if (is.na(median_read_len) || median_read_len <= 0) next

      # Skip pairs where either peak has no qualifying SNP position
      if (is.na(pk_a$snp_pos) || is.na(pk_b$snp_pos)) next

      gap_bp <- pk_b$snp_pos - pk_a$snp_pos
      if (is.na(gap_bp) || gap_bp > median_read_len) next     # too far apart

      # Jaccard on read sets
      n_shared <- length(intersect(reads_a, reads_b))
      n_union  <- length(union(reads_a, reads_b))
      jaccard  <- if (n_union > 0) n_shared / n_union else 0

      # Spanning reads: in the intersection AND have RLE data on both sides
      # of BOTH peak positions.
      spanning_ids <- intersect(reads_a, reads_b)

      # Zone boundaries
      zone_l_start <- -Inf
      zone_l_end   <- pk_a$snp_pos
      zone_m_start <- pk_a$snp_pos
      zone_m_end   <- pk_b$snp_pos
      zone_r_start <- pk_b$snp_pos
      zone_r_end   <- Inf

      n_spanning <- 0L
      edge_type  <- "unresolvable"

      if (length(spanning_ids) > 0) {
        span_df <- rt_df[read_id %in% spanning_ids &
                           as.character(chrom) == chr_name]

        if (nrow(span_df) > 0) {
          # Per-read zone classification
          state_list <- span_df[, {
            sL <- classify_zone_state(pos, IS_REF, zone_l_start, zone_l_end, zone_min_snps)
            sM <- classify_zone_state(pos, IS_REF, zone_m_start, zone_m_end, zone_min_snps)
            sR <- classify_zone_state(pos, IS_REF, zone_r_start, zone_r_end, zone_min_snps)
            .(state_L = sL, state_M = sM, state_R = sR)
          }, by = read_id]

          n_spanning <- nrow(state_list)
          edge_type  <- classify_edge_type(as.data.frame(state_list))
        }
      }

      # Fusion mode decision
      pair_key <- paste0(chr_name, "_", pk_a$peak_id, "_", pk_b$peak_id)

      fusion_mode <- if (edge_type %in% c("independent_events", "unresolvable")) {
        "none"
      } else if (jaccard >= jaccard_threshold &&
                 edge_type %in% c("gene_conversion", "crossover")) {
        "automatic"
      } else if (jaccard > 0 &&
                 edge_type %in% c("gene_conversion", "crossover", "ambiguous")) {
        "supervised"
      } else {
        "none"
      }

      # User-approved supervised fusions override to automatic
      if (!is.null(supervised_override) && pair_key %in% supervised_override) {
        if (fusion_mode == "supervised") fusion_mode <- "automatic"
      }

      all_pairs[[length(all_pairs) + 1L]] <- data.table(
        pair_key     = pair_key,
        chrom        = chr_name,
        peak_id_a    = pk_a$peak_id,
        peak_id_b    = pk_b$peak_id,
        snp_pos_a    = pk_a$snp_pos,
        snp_pos_b    = pk_b$snp_pos,
        gap_bp       = gap_bp,
        median_read_len_bp = round(median_read_len),
        n_reads_a    = length(reads_a),
        n_reads_b    = length(reads_b),
        n_shared     = n_shared,
        n_spanning   = n_spanning,
        jaccard      = round(jaccard, 4),
        edge_type    = edge_type,
        fusion_mode  = fusion_mode
      )
    }
  }

  if (length(all_pairs) == 0)
    return(list(peak_pairs = data.table(), fused_peaks = peaks_dt))

  pairs_dt <- rbindlist(all_pairs)

  # ── Graph-based connected components ──────────────────────────────────────
  auto_edges <- pairs_dt[fusion_mode == "automatic"]

  peaks_dt[, fusion_group := peak_id]   # default: each peak is its own group

  if (nrow(auto_edges) > 0) {
    g <- igraph::graph_from_data_frame(
      d        = auto_edges[, .(from = peak_id_a, to = peak_id_b)],
      directed = FALSE,
      vertices = data.frame(name = peaks_dt$peak_id)
    )
    comps <- igraph::components(g)
    # Map component membership back to peak IDs
    mem   <- comps$membership              # named by peak_id (as character)
    peaks_dt[, fusion_group := mem[as.character(peak_id)]]
  }

  # ── Compute fused peak coordinates per group ───────────────────────────────
  fused_coords <- peaks_dt[, .(
    fused_pos_bp    = round(mean(snp_pos,    na.rm = TRUE)),
    fused_start_bp  = min(peak_start, na.rm = TRUE),
    fused_end_bp    = max(peak_end,   na.rm = TRUE),
    n_sub_peaks     = .N,
    constituent_ids = paste(peak_id, collapse = ",")
  ), by = fusion_group]

  peaks_dt <- merge(peaks_dt, fused_coords, by = "fusion_group", all.x = TRUE)

  # Annotate with edge_type and fusion_mode from pairs where this peak is
  # involved (take the "strongest" relationship if multiple edges)
  edge_priority <- c(gene_conversion = 1, crossover = 2,
                     ambiguous = 3, independent_events = 4,
                     unresolvable = 5)

  get_peak_edge_info <- function(pid) {
    rel_pairs <- pairs_dt[peak_id_a == pid | peak_id_b == pid]
    if (nrow(rel_pairs) == 0)
      return(data.table(peak_id = pid, best_edge_type = NA_character_,
                        best_fusion_mode = NA_character_,
                        adjacent_pair_keys = NA_character_))
    rel_pairs[, priority := edge_priority[edge_type]]
    best <- rel_pairs[which.min(priority)]
    data.table(
      peak_id            = pid,
      best_edge_type     = best$edge_type,
      best_fusion_mode   = best$fusion_mode,
      adjacent_pair_keys = paste(rel_pairs$pair_key, collapse = ";")
    )
  }

  edge_info <- rbindlist(lapply(peaks_dt$peak_id, get_peak_edge_info))
  peaks_dt  <- merge(peaks_dt, edge_info, by = "peak_id", all.x = TRUE)

  list(
    peak_pairs  = pairs_dt,
    fused_peaks = peaks_dt
  )
}

# ─────────────────────────────────────────────
#                   UI
# ─────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well .form-group              { margin-top: 2px; margin-bottom: 2px; }
    .well .shiny-input-container   { margin-top: 2px; margin-bottom: 2px; }
    .well hr                       { margin-top: 2px; margin-bottom: 2px; }
    .well h3                       { margin-top: 2px; margin-bottom: 2px; }
    .fusion-auto   { background-color: #d4edda !important; }
    .fusion-super  { background-color: #fff3cd !important; }
    .fusion-none   { background-color: #f8d7da !important; }
    .fusion-unresolv { background-color: #e2e3e5 !important; }
  "))),
  titlePanel("ChimeraMapR: Detect Chimeric Haplotypes in Long Sequence Reads"),

  sidebarLayout(
    sidebarPanel(
      h3("Data Files"),
      fileInput("read_data_file",
                "Read Data File (CSV/GZ):",
                accept = c(".csv", ".gz", ".csv.gz")),

      fileInput("snp_data_file",
                "SNP Data File (CSV or VCF):",
                accept = c(".csv", ".vcf", ".vcf.gz", ".gz")),

      fileInput("chr_size_file",
                "Chromosome Size File (FAI):",
                accept = c(".fai", ".txt")),

      h3("Analysis Parameters"),
      textInput("sample_name",
                "Sample Name:",
                value = "Sample_01"),

      numericInput("mapq_cutoff",
                   "Minimum MAPQ Value:",
                   value = 20,
                   min = 0,
                   step = 1),

      numericInput("baseq_cutoff",
                   "Base Quality Minimum:",
                   value = 10,
                   min = 0,
                   step = 1),

      numericInput("min_run",
                   "Minimum Run Length:",
                   value = 2,
                   min = 1,
                   step = 1),
      helpText("Minimum consecutive same-allele calls to count as a run; increase for noisier data"),

      numericInput("min_peak_height",
                   "Minimum Peak Height:",
                   value = 10,
                   min = 1,
                   step = 1),

      numericInput("lambda",
                   "Whittaker Lambda (\u03bb):",
                   value = 1,
                   min   = 0.01,
                   step  = 0.5),
      helpText("Smoothness penalty for Whittaker smoother. Lower = tighter fit (preserves sharp peaks); higher = smoother curve"),

      hr(),
      h3("Peak Fusion Parameters"),

      numericInput("jaccard_threshold",
                   "Jaccard Threshold for Auto-Fusion:",
                   value = 0.20,
                   min   = 0.01,
                   max   = 1.0,
                   step  = 0.01),
      helpText("Minimum Jaccard index of peak-associated reads to trigger automatic peak fusion"),

      actionButton("run_analysis",
                   "Run Analysis",
                   class = "btn-primary btn-lg"),

      actionButton("run_fusion",
                   "Run Peak Fusion",
                   class = "btn-warning btn-lg"),
      helpText("Run after analysis. Re-run after approving supervised fusions below."),

      br(),
      tags$small(
        style = "color: gray;",
        paste("ChimeraMapR version", APP_VERSION)
      ),
      width = 3
    ),

    mainPanel(
      tabsetPanel(
        id = "main_tabs",

        tabPanel("Overview Plot",
                 helpText("Genome-wide overview of chimeric read coverage and Whittaker-smoothed signal across chromosomes."),
                 plotOutput(
                   "chr_plot",
                   height = "1200px"
                 ),
                 fluidRow(
                   column(3, downloadButton("download_plot",     "Download Plot (.png)")),
                   column(4, downloadButton("download_plot_rds", "Download R Object (.rds)"))
                 )
        ),

        tabPanel("Chromosome Plots",
                 h4("Per-Chromosome Coverage Plots"),
                 helpText("Select a chromosome tab, brush an x region, then click 'Plot Selected Region' to view chimeric reads in that interval."),
                 br(),
                 uiOutput("chr_plots_tabs")
        ),

        tabPanel("Peak Summary",
                 h4("Detected Peaks"),
                 tableOutput("peaks_table"),
                 downloadButton("download_peaks", "Download Peak Data")
        ),

        tabPanel("Peak Fusion",
                 h4("Peak Pair Analysis"),
                 helpText(
                   "Adjacent peak pairs are evaluated for shared reads and haplotype signatures. ",
                   "Green rows = automatic fusion. Yellow rows = supervised (check box to approve). ",
                   "Red rows = not fused (independent events). Grey rows = unresolvable."
                 ),
                 br(),
                 h5("Candidate Peak Pairs"),
                 uiOutput("fusion_review_ui"),
                 br(),
                 actionButton("apply_supervised", "Apply Checked Fusions", class = "btn-warning"),
                 helpText("After checking boxes above, click to apply and re-run fusion."),
                 br(), br(),
                 h5("Fused Peak Table"),
                 tableOutput("fused_peaks_table"),
                 downloadButton("download_fused_peaks",  "Download Fused Peak Data"),
                 downloadButton("download_peak_pairs",   "Download Pair Analysis")
        ),

        tabPanel("Individual Peak Plots",
                 h4("Detailed Peak Visualizations by Chromosome"),
                 helpText("Each peak shows all chimeric reads that intersect that position.
                          Colors indicate REF (blue) vs ALT (red) alleles."),
                 br(),
                 uiOutput("peak_plots_tabs")
        ),

        tabPanel("Curve Fits",
                 h4("Whittaker Smoother Parameters"),
                 tableOutput("span_table"),
                 helpText("Shows the lambda (\u03bb) value used for each chromosome in the analysis"),
                 br(),
                 h5("Export Whittaker fitted curves"),
                 helpText("Downloads the fitted Whittaker curves for all chromosomes from this run, including run parameters for later comparison across runs."),
                 downloadButton("download_curve_fits", "Download Curve Fits")
        ),

        tabPanel("Read Statistics",
                 h4("Analysis Summary"),
                 verbatimTextOutput("summary_stats"),
                 br(),
                 downloadButton("download_read_ids", "Download Chimeric Read IDs")
        ),

        tabPanel("About",
                 h4("About This Analysis"),
                 p("This application identifies chimeric reads in sequencing data by tracking
                   allele changes across chromosomes. Chimeric reads contain sequence from
                   multiple parental chromosomes and can indicate recombination events."),
                 p(strong("Version:"), APP_VERSION),
                 p(
                   "View the source code and README on",
                   tags$a("GitHub", href = "https://github.com/RobertJDReid/ChimeraMapR", target = "_blank"),
                   "."
                 ),
                 h5("Method:"),
                 tags$ul(
                   tags$li("Classifies base calls at SNP positions as REF or ALT alleles"),
                   tags$li("Uses run-length encoding to detect consecutive allele switches"),
                   tags$li("Identifies reads with multiple sustained allele changes"),
                   tags$li("Counts chimeric reads at each SNP position"),
                   tags$li("Applies Whittaker smoothing to identify peaks of per-read haplotype switches"),
                   tags$li("Groups adjacent peaks by shared reads and haplotype signature (gene conversion, crossover, independent events)")
                 ),
                 h5("Peak Fusion Edge Types:"),
                 tags$ul(
                   tags$li(strong("gene_conversion:"), "All spanning reads show A\u2192B\u2192A pattern (return to original haplotype). Fused automatically if Jaccard \u2265 threshold."),
                   tags$li(strong("crossover:"), "Spanning reads split into two complementary classes (A\u2192B and B\u2192A). Fused automatically if Jaccard \u2265 threshold."),
                   tags$li(strong("independent_events:"), "No spanning read returns to origin; divergent terminal haplotypes. Never fused."),
                   tags$li(strong("ambiguous:"), "Read patterns do not fit a clean model. Supervised review only."),
                   tags$li(strong("unresolvable:"), "No spanning reads with sufficient SNPs in all three zones. Not fused.")
                 ),
                 h5("Input Files:"),
                 tags$ul(
                   tags$li(strong("Read Data:"), "CSV file with SNP position information from BAM file (columns: chrom, pos, read_id, call, is_del, etc.). For csv files > 200 Mb, compress with ", em("gzip"), " or ", em("pigz"), " prior to upload."),
                   tags$li(strong("SNP Data:"), "CSV file with SNP positions (columns: CHROM, POS, REF, ALT), or a VCF file (plain or gzipped). Multi-allelic VCF sites are split into one row per ALT allele."),
                   tags$li(strong("Chromosome Size:"), "FAI index file with chromosome lengths")
                 ),
                 h5("Whittaker Lambda (\u03bb):"),
                 tags$ul(
                   tags$li(strong("Low \u03bb (0.01\u20131):"), "Tight fit \u2014 preserves sharp, narrow peaks well. Recommended for impulsive boundary-count signals."),
                   tags$li(strong("High \u03bb (10\u20131000):"), "Heavy smoothing \u2014 useful for broad signal or very noisy data, but will attenuate sharp peaks.")
                 )
        )
      ),
      width = 9
    )
  )
)

# ─────────────────────────────────────────────
#               SERVER
# ─────────────────────────────────────────────
server <- function(input, output, session) {

  # Set max upload size to 200 MB
  options(shiny.maxRequestSize = 200 * 1024^2)

  # Shared trigger for region plot building
  rv_trigger_region <- reactiveVal(0L)

  # Reactive values to store analysis results
  results <- reactiveValues(
    rt_df                 = NULL,
    transition_pos        = NULL,
    snp_coverage          = NULL,
    peaks_genomic         = NULL,
    snp_peaks             = NULL,
    chromosome_fits       = NULL,
    chr_span              = NULL,
    plot                  = NULL,
    chimeric_read_ids     = NULL,
    peak_plots_by_chr     = NULL,
    selected_region       = NULL,
    selected_regions      = list(),
    selected_region_plot  = NULL,
    selected_region_data  = NULL,
    # Fusion results (non-destructive — sit alongside originals)
    peak_pairs            = NULL,
    fused_peaks           = NULL
  )

  # Tracks which supervised pair_keys the user has checked
  supervised_approved <- reactiveVal(character(0))

  # ── RLE helper (operates on plain vectors) ──────────────────────
  rle_helper <- function(x) {
    r  <- rle(x)[[1]]
    rn <- rep(r, r)
    return(rn)
  }

  # ── Overview plot (genome-wide) ───────────────────────────────────────────────
  output$chr_plot <- renderPlot({
    req(results$snp_coverage, results$chromosome_fits)
    build_overview_plot(results)
  }, height = function() {
    req(results$snp_coverage)
    n_chr <- length(unique(results$snp_coverage$chrom))
    min(1600, max(400, n_chr * 120))
  })

  # ── Main analysis ────────────────────────────────────────────────────────────
  observeEvent(input$run_analysis, {

    req(input$read_data_file, input$snp_data_file, input$chr_size_file)

    # Reset outputs on re-run
    results$selected_region      <- NULL
    results$selected_region_plot <- NULL
    results$selected_region_data <- NULL
    results$peak_pairs           <- NULL
    results$fused_peaks          <- NULL
    supervised_approved(character(0))

    try(removeTab(inputId = "main_tabs", target = "Selected Region"), silent = TRUE)

    withProgress(message = "Processing data...", value = 0, {

      incProgress(0.1, detail = "Running analysis")

      restore_ext <- function(datapath, original_name) {
        suffix <- sub("^[^.]+", "", original_name)
        if (nchar(suffix) == 0) return(datapath)
        new_path <- paste0(datapath, suffix)
        file.copy(datapath, new_path, overwrite = TRUE)
        new_path
      }

      read_path <- restore_ext(input$read_data_file$datapath, input$read_data_file$name)
      snp_path  <- restore_ext(input$snp_data_file$datapath,  input$snp_data_file$name)
      chr_path  <- restore_ext(input$chr_size_file$datapath,  input$chr_size_file$name)

      res <- run_chimera_analysis(
        read_data_path  = read_path,
        snp_data_path   = snp_path,
        chr_size_path   = chr_path,
        sample_name     = input$sample_name,
        mapq_cutoff     = input$mapq_cutoff,
        baseq_cutoff    = input$baseq_cutoff,
        min_run         = input$min_run,
        min_peak_height = input$min_peak_height,
        lambda          = input$lambda,
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

      incProgress(0.9, detail = "Creating individual peak plots")

      peak_chrs <- character(0)
      if (!is.null(res$snp_peaks) && nrow(res$snp_peaks) > 0) {
        peak_chrs <- unique(as.character(res$snp_peaks$chrom))
      }

      peak_plots_list <- lapply(peak_chrs, function(chr_name) {

        chr_peaks <- res$snp_peaks[as.character(chrom) == chr_name]
        chr_peaks <- chr_peaks[!is.na(snp_pos)]
        chr_peaks <- chr_peaks[order(snp_pos)]

        if (nrow(chr_peaks) == 0) return(list(chromosome = chr_name, plots = list()))

        plots <- lapply(seq_len(nrow(chr_peaks)), function(pk_i) {
          pk       <- chr_peaks[pk_i]
          snp_p    <- pk$snp_pos
          pk_start <- pk$peak_start
          pk_end   <- pk$peak_end

          touching_ids <- res$transition_pos[
            as.character(chrom) == chr_name & pos >= pk_start & pos <= pk_end,
            unique(read_id)
          ]
          if (length(touching_ids) == 0) return(NULL)

          pad_bp     <- 5000L
          chr_pos    <- res$snp_coverage[as.character(chrom) == chr_name, pos]
          chr_min    <- min(chr_pos, na.rm = TRUE)
          chr_max    <- max(chr_pos, na.rm = TRUE)
          plot_start <- max(chr_min, pk_start - pad_bp)
          plot_end   <- min(chr_max, pk_end   + pad_bp)

          plot_df <- res$rt_df[
            as.character(chrom) == chr_name &
              read_id %in% touching_ids &
              pos >= plot_start &
              pos <= plot_end
          ]
          if (nrow(plot_df) == 0) return(NULL)
          setorder(plot_df, read_id, pos)

          peak_point_df <- data.frame(
            pos_kb      = snp_p / 1000,
            peak_start  = pk_start / 1000,
            peak_end    = pk_end   / 1000,
            peak_height = pk$peak_height
          )

          read_peak_df <- plot_df[pos >= pk_start & pos <= pk_end]

          seg_data <- NULL
          if (nrow(read_peak_df) > 0) {
            peak_summary <- read_peak_df[, .(
              REF = sum(IS_REF),
              num = .N
            ), by = pos][, allele_balance := REF / num][
              , SNP_call := data.table::fcase(
                  allele_balance < 0.2, "ALT",
                  allele_balance > 0.8, "REF",
                  default = "HET"
              )
            ]
            setorder(peak_summary, pos)
            peak_summary[, run := data.table::rleid(SNP_call)]
            seg_data <- peak_summary[, .(
              xmin     = min(pos) / 1000,
              SNP_call = SNP_call[1]
            ), by = run][order(xmin)][
              , xmax := data.table::shift(xmin, type = "lead",
                                          fill = pk_end / 1000)
            ][, SNP_call := factor(SNP_call, levels = c("ALT", "HET", "REF"))]
          }

          x_lims <- range(plot_df$pos / 1000)

          p_reads <- ggplot(plot_df, aes(x = pos / 1000, y = 1, colour = IS_REF)) +
            geom_point(size = 0.8) +
            geom_vline(
              data        = peak_point_df,
              aes(xintercept = pos_kb),
              colour      = "dodgerblue",
              linewidth   = 1,
              linetype    = 1,
              inherit.aes = FALSE
            ) +
            geom_vline(xintercept = pk_start / 1000,
                       color = "grey60", linewidth = 0.6, linetype = 2) +
            geom_vline(xintercept = pk_end / 1000,
                       color = "grey60", linewidth = 0.6, linetype = 2) +
            facet_grid(read_id ~ .) +
            scale_color_viridis_d(option = "turbo", begin = 0.87, end = 0.2) +
            scale_x_continuous(limits = x_lims, expand = expansion(mult = 0)) +
            theme_bw() +
            theme(
              axis.text.y      = element_blank(),
              axis.ticks.y     = element_blank(),
              axis.title.y     = element_blank(),
              axis.text.x      = element_blank(),
              axis.ticks.x     = element_blank(),
              axis.title.x     = element_blank(),
              legend.position  = "none",
              plot.background  = element_blank(),
              strip.background = element_blank(),
              panel.border     = element_rect(fill = NA, linewidth = 0.1, linetype = 3),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.text.y     = element_blank(),
              panel.spacing    = unit(0, "mm")
            ) +
            ggtitle(paste0(
              "Chr ", chr_name, "  \u2014  Peak ", pk_i,
              "  (SNP: ", round(snp_p / 1000, 2), " Kb;  Display: ",
              round(plot_start / 1000, 2), " \u2013 ",
              round(plot_end   / 1000, 2), " Kb)"
            ))

          if (!is.null(seg_data) && nrow(seg_data) > 0) {
            p_seg <- ggplot(seg_data) +
              geom_rect(
                aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
                    fill = SNP_call),
                alpha = 0.45
              ) +
              scale_fill_manual(
                values   = c(ALT = "firebrick", HET = "gray60", REF = "dodgerblue"),
                drop     = FALSE,
                na.value = "black"
              ) +
              coord_cartesian(
                xlim   = x_lims,
                ylim   = c(0, 1),
                expand = FALSE
              ) +
              labs(x = "Position (Kb)", y = "Haplotype\nregion") +
              theme_bw() +
              theme(
                panel.grid      = element_blank(),
                axis.text.y     = element_blank(),
                axis.ticks.y    = element_blank(),
                panel.border    = element_blank(),
                legend.position = "none"
              )
            p <- patchwork::wrap_plots(p_reads, p_seg, ncol = 1,
                                       heights = c(8, 1))
          } else {
            p <- p_reads + xlab("Position (Kbp)")
          }

          attr(p, "seg_data") <- seg_data
          p
        })

        plots <- Filter(Negate(is.null), plots)
        list(chromosome = chr_name, plots = plots)
      })

      results$peak_plots_by_chr <- Filter(
        function(x) length(x$plots) > 0,
        peak_plots_list
      )

      incProgress(1, detail = "Complete")
    })

    showNotification("Analysis complete. Click 'Run Peak Fusion' to evaluate adjacent peaks.", type = "message", duration = 5)
  })


  # ── Peak Fusion ───────────────────────────────────────────────────────────────
  observeEvent(input$run_fusion, {
    req(results$snp_peaks, results$rt_df, results$transition_pos)

    # zone_min_snps mirrors min_run (parameterizable later)
    zone_min_snps <- input$min_run

    withProgress(message = "Running peak fusion analysis...", value = 0.3, {

      fusion_res <- compute_peak_pairs(
        snp_peaks          = results$snp_peaks,
        rt_df              = results$rt_df,
        transition_pos     = results$transition_pos,
        jaccard_threshold  = input$jaccard_threshold,
        zone_min_snps      = zone_min_snps,
        supervised_override = supervised_approved()
      )

      results$peak_pairs  <- fusion_res$peak_pairs
      results$fused_peaks <- fusion_res$fused_peaks
    })

    showNotification("Peak fusion complete. Review results in the 'Peak Fusion' tab.", type = "message", duration = 4)
  })

  # Apply supervised fusions: collect checked boxes, store, re-run fusion
  observeEvent(input$apply_supervised, {
    req(results$peak_pairs)

    pairs <- results$peak_pairs
    super_pairs <- pairs[fusion_mode == "supervised", pair_key]

    approved <- Filter(function(pk) {
      cb_id <- paste0("fuse_cb_", gsub("[^A-Za-z0-9]", "_", pk))
      isTRUE(input[[cb_id]])
    }, super_pairs)

    supervised_approved(approved)

    # Re-run fusion with updated approvals
    zone_min_snps <- input$min_run

    withProgress(message = "Re-running fusion with approved pairs...", value = 0.3, {
      fusion_res <- compute_peak_pairs(
        snp_peaks           = results$snp_peaks,
        rt_df               = results$rt_df,
        transition_pos      = results$transition_pos,
        jaccard_threshold   = input$jaccard_threshold,
        zone_min_snps       = zone_min_snps,
        supervised_override = supervised_approved()
      )
      results$peak_pairs  <- fusion_res$peak_pairs
      results$fused_peaks <- fusion_res$fused_peaks
    })

    n_approved <- length(approved)
    showNotification(
      paste0(n_approved, " supervised fusion(s) applied."),
      type = "message", duration = 4
    )
  })

  # ── Fusion Review UI ─────────────────────────────────────────────────────────
  output$fusion_review_ui <- renderUI({
    if (is.null(results$peak_pairs) || nrow(results$peak_pairs) == 0) {
      return(p("No candidate pairs found, or fusion has not been run yet. Click 'Run Peak Fusion'."))
    }

    pairs <- results$peak_pairs

    row_tags <- lapply(seq_len(nrow(pairs)), function(i) {
      # Extract plain scalars from the data.table row
      fm    <- pairs$fusion_mode[[i]]
      et    <- pairs$edge_type[[i]]
      pk    <- pairs$pair_key[[i]]
      chr_i <- pairs$chrom[[i]]
      pos_a <- pairs$snp_pos_a[[i]]
      pos_b <- pairs$snp_pos_b[[i]]
      gap   <- pairs$gap_bp[[i]]
      jac   <- pairs$jaccard[[i]]
      n_sh  <- pairs$n_shared[[i]]
      n_sp  <- pairs$n_spanning[[i]]

      # Row background colour
      bg <- if (fm == "automatic") {
        "#d4edda"
      } else if (fm == "supervised") {
        "#fff3cd"
      } else if (et == "independent_events") {
        "#f8d7da"
      } else {
        "#e2e3e5"
      }

      cb_id   <- paste0("fuse_cb_", gsub("[^A-Za-z0-9]", "_", pk))
      show_cb <- fm == "supervised"

      decision_label <- if (fm == "automatic") {
        "\u2714 Auto-fused"
      } else if (fm == "supervised") {
        "Supervised"
      } else {
        "Not fused"
      }

      tags$div(
        style = paste0("background-color:", bg, "; padding:6px 10px; margin-bottom:4px; border-radius:4px; display:flex; align-items:center; gap:16px;"),
        tags$span(style = "min-width:60px;",  strong(chr_i)),
        tags$span(style = "min-width:120px;",
          paste0(round(pos_a / 1000, 1), " \u2013 ",
                 round(pos_b / 1000, 1), " Kb")),
        tags$span(style = "min-width:80px;",
          paste0("Gap: ", round(gap / 1000, 1), " Kb")),
        tags$span(style = "min-width:90px;",
          paste0("Jaccard: ", round(jac, 3))),
        tags$span(style = "min-width:60px;",
          paste0("Shared: ", n_sh)),
        tags$span(style = "min-width:70px;",
          paste0("Spanning: ", n_sp)),
        tags$span(style = "min-width:160px; font-style:italic;", et),
        tags$span(style = "min-width:100px;",
          strong(decision_label)),
        if (show_cb) {
          checkboxInput(cb_id, label = "Approve fusion",
                        value = pk %in% supervised_approved())
        } else {
          tags$span("")
        }
      )
    })

    tagList(
      # Legend
      tags$div(
        style = "margin-bottom:10px; font-size:0.85em; color:#555;",
        tags$span(style = "background:#d4edda; padding:2px 8px; border-radius:3px; margin-right:6px;", "Automatic"),
        tags$span(style = "background:#fff3cd; padding:2px 8px; border-radius:3px; margin-right:6px;", "Supervised (check to approve)"),
        tags$span(style = "background:#f8d7da; padding:2px 8px; border-radius:3px; margin-right:6px;", "Independent events"),
        tags$span(style = "background:#e2e3e5; padding:2px 8px; border-radius:3px;", "Unresolvable")
      ),
      # Column headers
      tags$div(
        style = "display:flex; gap:16px; font-weight:bold; padding:4px 10px; font-size:0.85em; color:#333;",
        tags$span(style="min-width:60px;",  "Chr"),
        tags$span(style="min-width:120px;", "Position"),
        tags$span(style="min-width:80px;",  "Gap"),
        tags$span(style="min-width:90px;",  "Jaccard"),
        tags$span(style="min-width:60px;",  "Shared"),
        tags$span(style="min-width:70px;",  "Spanning"),
        tags$span(style="min-width:160px;", "Edge Type"),
        tags$span(style="min-width:100px;", "Decision"),
        tags$span("Approve?")
      ),
      do.call(tagList, row_tags)
    )
  })

  # ── Fused peaks table ────────────────────────────────────────────────────────
  output$fused_peaks_table <- renderTable({
    req(results$fused_peaks)
    fp <- copy(results$fused_peaks)

    # Show one row per fusion group (the representative / fused coordinates)
    group_rep <- fp[, .(
      Chromosome          = chrom[1],
      `Fused Pos (Kb)`    = round(fused_pos_bp[1]   / 1000, 2),
      `Fused Start (Kb)`  = round(fused_start_bp[1] / 1000, 2),
      `Fused End (Kb)`    = round(fused_end_bp[1]   / 1000, 2),
      `Sub-peaks`         = n_sub_peaks[1],
      `Constituent IDs`   = constituent_ids[1],
      `Edge Type`         = ifelse(n_sub_peaks[1] == 1L, "singleton",
                                   ifelse(is.na(best_edge_type[1]), "—", best_edge_type[1])),
      `Fusion Mode`       = ifelse(n_sub_peaks[1] == 1L, "—",
                                   ifelse(is.na(best_fusion_mode[1]), "—", best_fusion_mode[1]))
    ), by = fusion_group]

    setorder(group_rep, Chromosome, `Fused Pos (Kb)`)
    group_rep[, fusion_group := NULL]
    group_rep
  })

  # ── Peak Summary table (original, unchanged) ─────────────────────────────────
  output$peaks_table <- renderTable({
    req(results$peaks_genomic)
    out <- copy(results$snp_peaks)
    out[, `:=`(
      peak_pos_kb   = round(peak_pos   / 1000, 2),
      peak_start_kb = round(peak_start / 1000, 2),
      peak_end_kb   = round(peak_end   / 1000, 2),
      peak_height   = round(peak_height, 2)
    )]
    out[, snp_pos_kb_str := ifelse(
      is.na(snp_pos),
      "None above cutoff",
      as.character(round(snp_pos / 1000, 2))
    )]
    out[, snp_n_str := ifelse(is.na(snp_n), "", as.character(snp_n))]
    out[, chimeric_reads_str := ifelse(
      is.na(chimeric_reads_at_snp), "", as.character(chimeric_reads_at_snp)
    )]
    out <- out[, .(
      Chromosome              = chrom,
      `Peak Position (Kb)`    = peak_pos_kb,
      `Qualifying SNP (Kb)`   = snp_pos_kb_str,
      `Raw Count at SNP`      = snp_n_str,
      `Peak Start (Kb)`       = peak_start_kb,
      `Peak End (Kb)`         = peak_end_kb,
      `Peak Height`           = peak_height,
      `Chimeric Reads at SNP` = chimeric_reads_str
    )]
    setorder(out, Chromosome, `Peak Position (Kb)`)
    out
  })

  # ── Dynamic UI: per-chromosome tabs of individual peak plots ─────────────────
  output$peak_plots_tabs <- renderUI({
    req(results$peak_plots_by_chr)

    if (is.null(results$peak_plots_by_chr) || length(results$peak_plots_by_chr) == 0)
      return(h4("No peaks detected in the analysis."))

    chr_tabs <- lapply(results$peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots    <- chr_data$plots

      if (length(plots) == 0) return(NULL)

      plot_outputs <- lapply(seq_along(plots), function(i) {
        tagList(
          plotOutput(paste0("peak_plot_", chr_name, "_", i), height = "600px"),
          fluidRow(
            column(3,
              downloadButton(
                paste0("dl_peak_png_", chr_name, "_", i),
                "Download Plot Image (.png)",
                class = "btn-sm btn-default"
              )
            ),
            column(4,
              downloadButton(
                paste0("dl_peak_rds_", chr_name, "_", i),
                "Download Plot Data (.rds)",
                class = "btn-sm btn-default"
              )
            )
          ),
          hr()
        )
      })

      tabPanel(
        paste("Chr", chr_name),
        h5(paste0("Chromosome ", chr_name, " - ", length(plots), " peak(s) detected")),
        hr(),
        do.call(tagList, plot_outputs)
      )
    })

    chr_tabs <- Filter(Negate(is.null), chr_tabs)
    if (length(chr_tabs) == 0) return(h4("No peaks detected in the analysis."))
    do.call(tabsetPanel, chr_tabs)
  })

  # Render individual peak plots with download handlers
  observe({
    req(results$peak_plots_by_chr)

    lapply(results$peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots    <- chr_data$plots

      lapply(seq_along(plots), function(i) {
        output_name <- paste0("peak_plot_", chr_name, "_", i)
        png_dl_id   <- paste0("dl_peak_png_", chr_name, "_", i)
        rds_dl_id   <- paste0("dl_peak_rds_", chr_name, "_", i)
        local({
          p    <- plots[[i]]
          .chr <- chr_name
          .i   <- i

          output[[output_name]] <- renderPlot({ p })

          output[[png_dl_id]] <- downloadHandler(
            filename = function() {
              paste0(input$sample_name, "_peak_chr", .chr, "_peak", .i,
                     "_", Sys.Date(), ".png")
            },
            content = function(file) {
              reads_plot <- if (inherits(p, "patchwork")) p[[1]] else p
              n_reads    <- length(unique(ggplot_build(reads_plot)$data[[1]]$group))
              has_seg    <- inherits(p, "patchwork")
              seg_extra  <- if (has_seg) 1.5 else 0
              plot_h     <- max(4, min(22, n_reads * 0.4 + 2)) + seg_extra
              ggsave(file, plot = p, width = 10, height = plot_h, dpi = 300)
            }
          )

          output[[rds_dl_id]] <- downloadHandler(
            filename = function() {
              paste0(input$sample_name, "_peak_chr", .chr, "_peak", .i,
                     "_", Sys.Date(), ".rds")
            },
            content = function(file) {
              reads_plot      <- if (inherits(p, "patchwork")) p[[1]] else p
              plot_data       <- as.data.table(reads_plot$data)
              plot_data       <- plot_data[, .(chrom, pos, read_id, IS_REF)]
              peak_points_dt  <- as.data.table(reads_plot$layers[[2]]$data)
              seg_dt          <- attr(p, "seg_data")
              saveRDS(
                list(
                  plot_data   = plot_data,
                  peak_points = peak_points_dt,
                  seg_data    = seg_dt,
                  chromosome  = .chr,
                  peak_number = .i,
                  sample_name = input$sample_name,
                  app_version = APP_VERSION
                ),
                file
              )
            }
          )
        })
      })
    })
  })

  # ── Chromosome plots tab ──────────────────────────────────────────────────────
  output$chr_plots_tabs <- renderUI({
    req(results$snp_coverage, results$chromosome_fits)

    chr_levels <- levels(results$snp_coverage$chrom)
    if (length(chr_levels) == 0) return(h4("Run analysis first."))

    tabs <- lapply(chr_levels, function(chr_name) {
      plot_id  <- paste0("chr_cov_plot_", chr_name)
      brush_id <- paste0("chr_cov_brush_", chr_name)
      btn_id   <- paste0("chr_plot_btn_", chr_name)
      txt_id   <- paste0("chr_region_text_", chr_name)

      tabPanel(
        chr_name,
        br(),
        plotOutput(
          plot_id,
          height = "400px",
          brush  = brushOpts(id = brush_id, direction = "x", resetOnNew = TRUE, clip = FALSE)
        ),
        br(),
        fluidRow(
          column(3, actionButton(btn_id, "Plot Selected Region", class = "btn-primary")),
          column(9, verbatimTextOutput(txt_id))
        ),
        br()
      )
    })

    do.call(tabsetPanel, c(list(id = "chr_plots_subtabs"), tabs))
  })

  # Render per-chromosome coverage plots and wire up brush observers
  observe({
    req(results$snp_coverage, results$chromosome_fits)

    chr_levels    <- levels(results$snp_coverage$chrom)
    snp_cov_all   <- copy(results$snp_coverage)
    fits_all      <- copy(results$chromosome_fits)
    peaks_all     <- results$peaks_genomic
    snp_peaks_all <- if (!is.null(results$snp_peaks) && nrow(results$snp_peaks) > 0)
      copy(results$snp_peaks)[, chrom := as.character(chrom)]
    else NULL
    snp_cov_all[, chrom := as.character(chrom)]
    fits_all[,    chrom := as.character(chrom)]
    if (!is.null(peaks_all) && nrow(peaks_all) > 0)
      peaks_all <- copy(peaks_all)[, chrom := as.character(chrom)]

    lapply(chr_levels, function(chr_name) {
      chr_c    <- as.character(chr_name)
      plot_id  <- paste0("chr_cov_plot_", chr_c)
      brush_id <- paste0("chr_cov_brush_", chr_c)
      btn_id   <- paste0("chr_plot_btn_",  chr_c)
      txt_id   <- paste0("chr_region_text_", chr_c)
      local({
        .chr           <- chr_c
        .plot_id       <- plot_id
        .brush_id      <- brush_id
        .btn_id        <- btn_id
        .txt_id        <- txt_id
        .snp           <- snp_cov_all[chrom == .chr]
        .fits          <- fits_all[chrom == .chr]
        .peaks         <- if (!is.null(peaks_all) && nrow(peaks_all) > 0)
                            peaks_all[chrom == .chr] else NULL
        .snp_peaks_chr <- if (!is.null(snp_peaks_all) && nrow(snp_peaks_all) > 0)
                            snp_peaks_all[chrom == .chr]
                          else NULL

        output[[.plot_id]] <- renderPlot({
          p <- ggplot(.snp, aes(x = pos_kb, y = n)) +
            geom_line(
              data  = .fits,
              aes(x = uniform_pos / 1000, y = uniform_fit),
              color = "firebrick", linewidth = 0.8, alpha = 0.7
            ) +
            geom_point(color = "black", alpha = 0.5, size = 0.8, shape = 21) +
            scale_x_continuous(minor_breaks = seq(0, 1600, 100)) +
            xlab("Position (Kbp)") +
            ylab("Number of Reads") +
            ylim(0, max(30, max(.snp$n))) +
            ggtitle(paste("Chromosome", .chr)) +
            theme_bw() +
            theme(
              panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
              panel.grid.major.x = element_line(linewidth = 0.05, color = "red")
            )

          if (!is.null(.snp_peaks_chr) && nrow(.snp_peaks_chr) > 0) {
            peak_highlight <- merge(
              .snp_peaks_chr[, .(pos = snp_pos)],
              .snp[, .(pos, pos_kb, n)],
              by = "pos"
            )
            if (nrow(peak_highlight) > 0)
              p <- p + geom_point(
                data        = peak_highlight,
                aes(x = pos_kb, y = n),
                color       = "black",
                fill        = "dodgerblue",
                size        = 3,
                shape       = 21,
                alpha       = 0.9,
                inherit.aes = FALSE
              )
          }
          p
        })

        observeEvent(input[[.brush_id]], {
          brush <- input[[.brush_id]]
          if (is.null(brush)) return(NULL)

          chr_max_kb <- max(.snp$pos_kb, na.rm = TRUE)
          if (brush$xmin > chr_max_kb * 2 || brush$xmax > chr_max_kb * 10) return(NULL)

          xmin_kb <- min(brush$xmin, brush$xmax)
          xmax_kb <- max(brush$xmin, brush$xmax)
          xmin_bp <- as.integer(round(xmin_kb * 1000))
          xmax_bp <- as.integer(round(xmax_kb * 1000))

          results$selected_regions[[.chr]] <- list(
            chrom = .chr,
            start = xmin_bp,
            end   = xmax_bp
          )
        }, ignoreNULL = TRUE)

        output[[.txt_id]] <- renderText({
          reg <- results$selected_regions[[.chr]]
          if (is.null(reg)) return("No region selected.")
          paste0(
            "Selected: ", format(reg$start, big.mark = ","),
            " - ", format(reg$end, big.mark = ","), " bp  (",
            round((reg$end - reg$start) / 1000, 2), " kb)"
          )
        })

        observeEvent(input[[.btn_id]], {
          reg <- results$selected_regions[[.chr]]
          req(!is.null(reg))
          results$selected_region <- reg
          rv_trigger_region(isolate(rv_trigger_region()) + 1L)
        }, ignoreNULL = TRUE)
      })
    })
  })

  # Whittaker lambda table
  output$span_table <- renderTable({
    req(results$chr_span)
    out <- copy(results$chr_span)
    out[, length_kb := round(length / 1000, 1)]
    out <- out[, .(
      Chromosome    = chrom,
      `Length (Kb)` = length_kb,
      `SNP Count`   = snp_count,
      `Lambda` = lambda
    )]
    setorder(out, Chromosome)
    out
  })

  # Summary statistics
  output$summary_stats <- renderText({
    req(results$rt_df, results$snp_coverage, results$peaks_genomic)

    n_chimeric_reads <- uniqueN(results$rt_df$read_id)
    n_chromosomes    <- uniqueN(results$snp_coverage$chrom)
    n_peaks          <- nrow(results$peaks_genomic)
    total_snps       <- nrow(results$snp_coverage)
    boundary_snps    <- results$snp_coverage[n > 0, .N]
    total_boundaries <- sum(results$snp_coverage$n)

    fusion_summary <- if (!is.null(results$peak_pairs) && nrow(results$peak_pairs) > 0) {
      pp <- results$peak_pairs
      paste0(
        "\nPeak Fusion Summary:\n",
        "  Candidate pairs evaluated: ", nrow(pp), "\n",
        "  Auto-fused pairs:          ", sum(pp$fusion_mode == "automatic"), "\n",
        "  Supervised (pending):      ", sum(pp$fusion_mode == "supervised"), "\n",
        "  Not fused (indep/unresolv):", sum(pp$fusion_mode == "none"), "\n"
      )
    } else {
      "\nPeak Fusion: not yet run (click 'Run Peak Fusion')\n"
    }

    paste0(
      "Sample: ", input$sample_name, "\n\n",
      "Chimeric Reads Detected: ",             n_chimeric_reads,  "\n",
      "Chromosomes Analyzed: ",                n_chromosomes,     "\n",
      "Total SNP Positions: ",                 total_snps,        "\n",
      "SNP Positions with Transition Boundaries: ", boundary_snps,
      " (", round(100 * boundary_snps / total_snps, 1), "%)\n",
      "Total Transition Boundary Events: ",    total_boundaries,  "\n",
      "Peaks Detected: ",                      n_peaks,           "\n",
      fusion_summary,
      "\nAnalysis Parameters:\n",
      "  MAPQ Cutoff: ",              input$mapq_cutoff,     "\n",
      "  Base Quality Cutoff: ",      input$baseq_cutoff,    "\n",
      "  Min Run Length: ",           input$min_run,         "\n",
      "  Min Peak Height: ",          input$min_peak_height, "\n",
      "  Whittaker Lambda (\u03bb): ",input$lambda,          "\n",
      "  Jaccard Threshold: ",        input$jaccard_threshold
    )
  })

  # ── Shared helper: build and display the selected-region read plot ───────────
  build_region_plot <- function() {
    req(results$selected_region, results$rt_df, results$transition_pos)

    reg <- results$selected_region

    touching_ids <- results$transition_pos[
      as.character(chrom) == reg$chrom & pos >= reg$start & pos <= reg$end,
      unique(read_id)
    ]

    if (length(touching_ids) == 0) {
      showNotification("No chimeric reads overlap the selected region.", type = "warning")
      return(NULL)
    }

    pad_bp <- 5000L

    chr_positions <- results$snp_coverage[as.character(chrom) == reg$chrom, pos]
    chr_min <- min(chr_positions, na.rm = TRUE)
    chr_max <- max(chr_positions, na.rm = TRUE)

    plot_start <- max(chr_min, reg$start - pad_bp)
    plot_end   <- min(chr_max, reg$end + pad_bp)

    plot_df <- results$rt_df[
      as.character(chrom) == reg$chrom &
        read_id %in% touching_ids &
        pos >= plot_start &
        pos <= plot_end
    ]

    if (nrow(plot_df) == 0) {
      showNotification("No plotted points available in the padded selected region.", type = "warning")
      return(NULL)
    }

    setorder(plot_df, read_id, pos)
    results$selected_region_data <- copy(plot_df)

    p <- ggplot(plot_df, aes(x = pos / 1000, y = 1, colour = IS_REF)) +
      geom_vline(xintercept = reg$start / 1000,
                 color = "grey60", linewidth = 0.8, linetype = 2) +
      geom_vline(xintercept = reg$end / 1000,
                 color = "grey60", linewidth = 0.8, linetype = 2) +
      geom_point() +
      facet_grid(read_id ~ .) +
      scale_color_viridis_d(option = "turbo", begin = 0.87, end = 0.2) +
      theme_bw() +
      theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.title.y     = element_blank(),
        legend.position  = "none",
        plot.background  = element_blank(),
        strip.background = element_blank(),
        panel.border     = element_rect(linewidth = 0.1, linetype = 3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y     = element_blank(),
        panel.spacing    = unit(0, "mm")
      ) +
      xlab("Position (Kbp)") +
      ggtitle(paste0(
        "Selected Region: Chr ", reg$chrom,
        "  (Selected: ", round(reg$start / 1000, 2), " - ",
        round(reg$end   / 1000, 2), " Kb;  Display: ",
        round(plot_start / 1000, 2), " - ",
        round(plot_end   / 1000, 2), " Kb)"
      ))

    results$selected_region_plot <- p

    if (!"Selected Region" %in% isolate(input$main_tabs)) {
      insertTab(
        inputId  = "main_tabs",
        target   = "Overview Plot",
        position = "after",
        select   = TRUE,
        tabPanel(
          title = "Selected Region",
          value = "Selected Region",
          h4("Regional Read Plot"),
          helpText("Shows all chimeric reads touching the selected interval, plotted in a padded display window."),
          plotOutput("selected_region_plot", height = "800px"),
          br(),
          fluidRow(
            column(3, downloadButton("download_selected_region_plot", "Download Plot Image (.png)")),
            column(4, downloadButton("download_selected_region_rds",  "Download Plot Data (.rds)")),
            column(3, actionButton("close_selected_region_tab", "Close Tab", class = "btn-secondary"))
          )
        )
      )
    } else {
      updateTabsetPanel(session, "main_tabs", selected = "Selected Region")
    }
  }

  observeEvent(rv_trigger_region(), {
    if (rv_trigger_region() > 0L) build_region_plot()
  }, ignoreInit = TRUE)

  output$selected_region_plot <- renderPlot({
    req(results$selected_region_plot)
    results$selected_region_plot
  })

  observeEvent(input$close_selected_region_tab, {
    try(removeTab(inputId = "main_tabs", target = "Selected Region"), silent = TRUE)
    results$selected_region_plot <- NULL
    results$selected_region_data <- NULL
  })

  # ── Download handlers ────────────────────────────────────────────────────────

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

  output$download_plot_rds <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chromosome_tracking_", Sys.Date(), ".rds"),
    content  = function(file) {
      req(results$snp_coverage, results$chromosome_fits)
      saveRDS(
        list(
          snp_coverage    = copy(results$snp_coverage),
          chromosome_fits = copy(results$chromosome_fits),
          peaks_genomic   = results$peaks_genomic,
          snp_peaks       = results$snp_peaks,
          sample_name     = input$sample_name,
          app_version     = APP_VERSION
        ),
        file
      )
    }
  )

  output$download_peaks <- downloadHandler(
    filename = function() paste0(input$sample_name, "_peaks_", Sys.Date(), ".csv"),
    content  = function(file) fwrite(results$snp_peaks, file)
  )

  output$download_fused_peaks <- downloadHandler(
    filename = function() paste0(input$sample_name, "_fused_peaks_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(results$fused_peaks)
      fwrite(results$fused_peaks, file)
    }
  )

  output$download_peak_pairs <- downloadHandler(
    filename = function() paste0(input$sample_name, "_peak_pairs_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(results$peak_pairs)
      fwrite(results$peak_pairs, file)
    }
  )

  output$download_read_ids <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chimeric_read_ids_", Sys.Date(), ".txt"),
    content  = function(file) writeLines(results$chimeric_read_ids, file)
  )

  output$download_selected_region_plot <- downloadHandler(
    filename = function() {
      req(results$selected_region)
      reg <- results$selected_region
      paste0(input$sample_name, "_selected_region_chr", reg$chrom,
             "_", reg$start, "_", reg$end, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(results$selected_region_plot)
      ggsave(file, plot = results$selected_region_plot, width = 10, height = 12, dpi = 300)
    }
  )

  output$download_selected_region_rds <- downloadHandler(
    filename = function() {
      req(results$selected_region)
      reg <- results$selected_region
      paste0(input$sample_name, "_selected_region_chr", reg$chrom,
             "_", reg$start, "_", reg$end, "_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(results$selected_region_data, results$selected_region)
      reg       <- results$selected_region
      plot_data <- as.data.table(results$selected_region_data)
      plot_data <- plot_data[, .(chrom, pos, read_id, IS_REF, ALLELE)]
      saveRDS(
        list(
          plot_data       = plot_data,
          selected_region = reg,
          sample_name     = input$sample_name,
          app_version     = APP_VERSION
        ),
        file
      )
    }
  )

  output$download_curve_fits <- downloadHandler(
    filename = function() paste0(input$sample_name, "_curve_fits_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(results$chromosome_fits)
      fwrite(results$chromosome_fits, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server)
