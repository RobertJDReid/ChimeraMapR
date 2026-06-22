library(shiny)
if (!requireNamespace("igraph", quietly = TRUE))
  stop("Package 'igraph' is required. Install with: install.packages('igraph')")
library(igraph)

source("chimera_functions.R")   # loads all packages + whittaker + run_chimera_analysis etc.
source("loh_chain_analysis.R")  # loads chain functions

# APP_VERSION is now defined inside chimera_functions.R

# Null-coalescing operator (base R does not provide one)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# build_fused_peak_plots()
#   Builds the per-chromosome list of fused peak group plots, mirroring the
#   structure of peak_plots_by_chr.  Called from both run_fusion and
#   apply_supervised observers so the plots stay in sync with fusion results.
#
#   Returns a list of: list(chromosome, plots)  — same shape as
#   results$peak_plots_by_chr so the renderUI / observe pattern is reusable.
# ---------------------------------------------------------------------------
build_fused_peak_plots <- function(fused_peaks, rt_df, transition_pos,
                                   snp_coverage, zone_min_snps) {

  if (is.null(fused_peaks) || nrow(fused_peaks) == 0) return(list())

  fp <- copy(fused_peaks)
  fp[, chrom := as.character(chrom)]

  # One representative row per fusion group (same as fused_peaks_table logic)
  group_rep <- fp[, .(
    chrom           = chrom[1],
    fused_pos_bp    = fused_pos_bp[1],
    fused_start_bp  = fused_start_bp[1],
    fused_end_bp    = fused_end_bp[1],
    n_sub_peaks     = n_sub_peaks[1],
    constituent_ids = constituent_ids[1],
    best_edge_type  = ifelse(n_sub_peaks[1] == 1L, "singleton",
                             ifelse(is.na(best_edge_type[1]), "—", best_edge_type[1])),
    best_fusion_mode = ifelse(n_sub_peaks[1] == 1L, "—",
                              ifelse(is.na(best_fusion_mode[1]), "—", best_fusion_mode[1])),
    # Collect all constituent peak_ids as integer vector for read lookup
    constituent_peak_ids = list(peak_id)
  ), by = fusion_group]

  setorder(group_rep, chrom, fused_pos_bp)
  
  # add to avoid NAs creeping into logic
  group_rep <- group_rep[!is.na(fused_pos_bp) &
                           !is.na(fused_start_bp) &
                           !is.na(fused_end_bp)]

  fused_chrs <- unique(group_rep$chrom)

  plots_by_chr <- lapply(fused_chrs, function(chr_name) {

    chr_groups <- group_rep[chrom == chr_name]

    plots <- lapply(seq_len(nrow(chr_groups)), function(gi) {

      grp         <- chr_groups[gi]
      snp_p       <- grp$fused_pos_bp
      fused_start <- grp$fused_start_bp
      fused_end   <- grp$fused_end_bp
      n_sub       <- grp$n_sub_peaks
      edge_lbl    <- grp$best_edge_type

      # Union of reads from all constituent sub-peaks (by their peak windows)
      # We use the original sub-peak boundaries stored in fused_peaks rows
      sub_peak_ids <- grp$constituent_peak_ids[[1]]
      sub_peak_rows <- fp[peak_id %in% sub_peak_ids & chrom == chr_name]

      touching_ids <- unique(unlist(lapply(seq_len(nrow(sub_peak_rows)), function(si) {
        spk <- sub_peak_rows[si]
        transition_pos[
          as.character(chrom) == chr_name &
            pos >= spk$peak_start & pos <= spk$peak_end,
          unique(read_id)
        ]
      })))

      if (length(touching_ids) == 0) return(NULL)

      # Re-classify using fused window boundaries
      hap <- classify_fused_peak_haplotype(
        fg            = grp,
        chr_name      = chr_name,
        rt_df         = rt_df,
        touching_ids  = touching_ids,
        zone_min_snps = zone_min_snps
      )

      seg_data <- hap$seg_data

      pad_bp     <- 5000L
      chr_pos    <- snp_coverage[as.character(chrom) == chr_name, pos]
      chr_min    <- min(chr_pos, na.rm = TRUE)
      chr_max    <- max(chr_pos, na.rm = TRUE)
      plot_start <- max(chr_min, fused_start - pad_bp)
      plot_end   <- min(chr_max, fused_end   + pad_bp)

      plot_df <- rt_df[
        as.character(chrom) == chr_name &
          read_id %in% touching_ids &
          pos >= plot_start &
          pos <= plot_end
      ]
      if (nrow(plot_df) == 0) return(NULL)
      setorder(plot_df, read_id, pos)

      x_lims <- range(plot_df$pos / 1000)

      hap_label_str <- gsub("_", " ", hap$label)
      expanded_note <- if (hap$expanded) " [window expanded]" else ""
      sub_note      <- if (n_sub > 1L)
        paste0("  \u2014  ", n_sub, " sub-peaks fused  [", edge_lbl, "]")
      else
        "  \u2014  singleton"

      plot_title <- paste0(
        "Chr ", chr_name, "  \u2014  Fused Group ", gi,
        "  (SNP: ", round(snp_p / 1000, 2), " Kb;  Display: ",
        round(plot_start / 1000, 2), " \u2013 ",
        round(plot_end   / 1000, 2), " Kb)",
        sub_note,
        "\n\u25b6 Haplotype classification: ", hap_label_str, expanded_note
      )

      # Sub-peak anchor lines (one vertical per constituent snp_pos)
      sub_snp_pos <- sub_peak_rows$snp_pos
      sub_snp_pos <- sub_snp_pos[!is.na(sub_snp_pos)]

      peak_point_df <- data.frame(
        pos_kb      = snp_p / 1000,
        peak_start  = fused_start / 1000,
        peak_end    = fused_end   / 1000
      )

      p_reads <- ggplot(plot_df, aes(x = pos / 1000, y = 1, colour = IS_REF)) +
        geom_point(size = 0.8) +
        # Fused centre position (blue solid)
        geom_vline(
          data        = peak_point_df,
          aes(xintercept = pos_kb),
          colour      = "dodgerblue",
          linewidth   = 1,
          linetype    = 1,
          inherit.aes = FALSE
        ) +
        # Fused boundary lines (grey dashed)
        geom_vline(xintercept = fused_start / 1000,
                   color = "grey60", linewidth = 0.6, linetype = 2) +
        geom_vline(xintercept = fused_end / 1000,
                   color = "grey60", linewidth = 0.6, linetype = 2)

      # Individual sub-peak SNP positions (purple dotted) when group is fused
      if (n_sub > 1L && length(sub_snp_pos) > 0) {
        p_reads <- p_reads +
          geom_vline(
            data        = data.frame(xint = sub_snp_pos / 1000),
            aes(xintercept = xint),
            colour      = "mediumpurple",
            linewidth   = 0.6,
            linetype    = 3,
            inherit.aes = FALSE
          )
      }

      p_reads <- p_reads +
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
        ggtitle(plot_title)

      if (hap$expanded) {
        p_reads <- p_reads +
          geom_vline(xintercept = hap$win_start / 1000,
                     color = "darkorange", linewidth = 0.7, linetype = 3) +
          geom_vline(xintercept = hap$win_end / 1000,
                     color = "darkorange", linewidth = 0.7, linetype = 3)
      }

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
          annotate("text",
                   x     = mean(c(hap$win_start, hap$win_end)) / 1000,
                   y     = 0.5,
                   label = hap_label_str,
                   size  = 3,
                   color = "black",
                   fontface = "bold") +
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

      attr(p, "seg_data")        <- seg_data
      attr(p, "haplotype_label") <- hap$label
      p
    })

    plots <- Filter(Negate(is.null), plots)
    list(chromosome = chr_name, plots = plots)
  })

  Filter(function(x) length(x$plots) > 0, plots_by_chr)
}

# ---------------------------------------------------------------------------
# extract_fused_haplotype_labels()
#   Walks the list returned by build_fused_peak_plots and extracts the
#   haplotype_label attr from each plot.  Returns a data.table with one row
#   per fusion group containing:
#     chrom, group_index (within chromosome), fusion_group (integer ID),
#     fused_classification (character label from the post-fusion window)
#
#   This is the authoritative post-fusion classification source for the
#   blended peak summary table.
# ---------------------------------------------------------------------------
extract_fused_haplotype_labels <- function(fused_peak_plots_by_chr, fused_peaks) {

  if (is.null(fused_peak_plots_by_chr) || length(fused_peak_plots_by_chr) == 0)
    return(NULL)
  if (is.null(fused_peaks) || nrow(fused_peaks) == 0)
    return(NULL)

  # Build a lookup: (chrom, within-chr group index) -> fusion_group integer
  # This mirrors the ordering used in build_fused_peak_plots (setorder chrom, fused_pos_bp)
  fp <- copy(fused_peaks)
  fp[, chrom := as.character(chrom)]

  group_rep <- fp[, .(
    chrom          = chrom[1],
    fused_pos_bp   = fused_pos_bp[1],
    n_sub_peaks    = n_sub_peaks[1]
  ), by = fusion_group]
  group_rep <- group_rep[!is.na(fused_pos_bp)]
  setorder(group_rep, chrom, fused_pos_bp)

  # Add a within-chromosome index (gi) matching build_fused_peak_plots loop order
  group_rep[, gi := seq_len(.N), by = chrom]

  rows <- list()
  for (chr_data in fused_peak_plots_by_chr) {
    chr_name <- chr_data$chromosome
    plots    <- chr_data$plots
    for (i in seq_along(plots)) {
      p   <- plots[[i]]
      lbl <- attr(p, "haplotype_label")
      if (is.null(lbl)) lbl <- NA_character_

      # Match back to the fusion_group integer via (chrom, gi)
      match_row <- group_rep[chrom == chr_name & gi == i]
      fg        <- if (nrow(match_row) == 1L) match_row$fusion_group else NA_integer_

      rows[[length(rows) + 1L]] <- data.table(
        chrom                 = chr_name,
        group_index           = i,
        fusion_group          = fg,
        fused_classification  = lbl
      )
    }
  }

  if (length(rows) == 0) return(NULL)
  rbindlist(rows, fill = TRUE)
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

      textInput("strain_ref",
                "REF Strain Name (optional):",
                value = "",
                placeholder = "e.g. S288C"),
      textInput("strain_alt",
                "ALT Strain Name (optional):",
                value = "",
                placeholder = "e.g. SK1"),
      helpText("Strain names appear in LOH band legends. Leave blank to use generic REF / ALT labels."),

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
      
      # LOH chain params
      
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

      numericInput("peak_pad_bp",
                   "Peak Association Padding (bp):",
                   value = 200, min = 0, step = 50),
      helpText("A peak is attached to a token junction if its SNP position falls within this distance of the token boundary."),

      numericInput("homog_frac",
                   "NCO-GC Homogeneity Fraction:",
                   value = 0.80, min = 0.5, max = 1.0, step = 0.05),
      helpText("Fraction of spanning reads that must share the return pattern (AAABBBAAA) to call NCO-GC."),
      
      actionButton("run_chain",
                   "Run Event Caller",
                   class = "btn-success btn-lg"),
      helpText("Run after Peak Fusion. Calls recombination events from LOH + peak chain."),

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
                   "Peaks already classified as gene conversion or internal crossover are excluded from fusion ",
                   "regardless of Jaccard score; only binary (or unlabelled) peaks are candidates. ",
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

        tabPanel("Fused Peak Plots",
                 h4("Post-Fusion Peak Visualizations by Chromosome"),
                 helpText(
                   "One plot per fusion group using the merged window boundaries and the union of reads ",
                   "from all constituent sub-peaks. Haplotype classification is re-run on the fused window. ",
                   "Blue solid line = fused group SNP centre; grey dashed lines = fused boundaries; ",
                   "purple dotted lines = individual sub-peak positions (fused groups only). ",
                   "Run 'Run Peak Fusion' first to populate this tab."
                 ),
                 br(),
                 uiOutput("fused_peak_plots_tabs")
        ),

        tabPanel("Post Fusion Peak Summary",
                 h4("Blended Peak Summary — Initial Peaks with Post-Fusion Classification"),
                 helpText(
                   "This table combines the initial peak positions and metrics with updated haplotype ",
                   "classifications after peak fusion. Singleton peaks retain their original classification. ",
                   "Fused groups show the re-classified haplotype from the merged window. ",
                   "Run 'Run Peak Fusion' first to populate the fused classifications."
                 ),
                 br(),
                 uiOutput("post_fusion_summary_legend"),
                 br(),
                 tableOutput("post_fusion_peak_summary_table"),
                 br(),
                 downloadButton("download_post_fusion_summary", "Download Post-Fusion Summary (.csv)")
        ),

        tabPanel("LOH Regions",
                 h4("Loss of Heterozygosity Regions"),
                 helpText("Fixed-haplotype LOH segments (REF_fixed and ALT_fixed) identified by a beta-binomial EM mixture + Viterbi HMM on per-position allele balance. Run analysis first to populate this table."),
                 br(),
                 tableOutput("loh_regions_table"),
                 br(),
                 downloadButton("download_loh_regions", "Download LOH Regions (.csv)")
        ),
        
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
                 h5("Peak-Class Fusion Guard:"),
                 p("Before edge-type scoring, each peak's haplotype label (assigned during the initial analysis) is checked.",
                   "Eligibility rules are defined in ", strong("FUSION_HEURISTICS"), " and enforced by ",
                   strong("peak_is_fusion_eligible()"), ".",
                   "Peaks classified as ", strong("gene_conversion"), " or ", strong("internal_crossover"),
                   " are treated as self-contained events and will ", strong("never"), " be merged with a neighbour,",
                   "regardless of Jaccard score or edge type.",
                   "Only peaks labelled ", strong("binary"), " (or still unlabelled / undefined) are eligible for fusion.",
                   "The fusion mode decision (automatic / supervised / none) is made by ",
                   strong("decide_fusion_mode()"), " — edit that function or ", strong("FUSION_HEURISTICS"),
                   " to change scoring criteria without touching the main analysis loop."),
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
    peak_pairs              = NULL,
    fused_peaks             = NULL,
    fused_peak_plots_by_chr = NULL,
    fused_hap_labels        = NULL,  # data.table: fusion_group, chrom, fused_classification
    loh_map                 = NULL,   # snp_table from compute_loh_map()$snp_table; per-SNP; used by gap_has_loh()
    loh_segments            = NULL,   # loh_segments from compute_loh_map()$loh_segments; used by plots
    coverage_segments       = NULL,   # coverage_segments from compute_coverage_map(); chromosome-wide depth_ratio
    coverage_table          = NULL,   # coverage_table from compute_coverage_map(); per-position depth for terminal-event rules
    strain_ref              = "",     # display name for REF haplotype (from UI input)
    strain_alt              = "",      # display name for ALT haplotype (from UI input)
    chain_result   = NULL,   # full output of run_chain_analysis()
    event_table    = NULL,   # flat data.table of called events
    chain_params   = NULL,   # params list used for the last chain run
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
    results$loh_segments   # explicit reactive dep — triggers re-render after LOH finishes
    results$event_table    # explicit reactive dep — triggers re-render after chain caller finishes
    build_overview_plot(results)
  }, height = function() {
    req(results$snp_coverage)
    n_chr   <- length(unique(results$snp_coverage$chrom))
    has_loh <- !is.null(results$loh_segments) &&
               nrow(results$loh_segments) > 0 &&
               any(results$loh_segments$loh_state %in% c("REF_fixed", "ALT_fixed"), na.rm = TRUE)
    # Add ~15px per chromosome for the LOH strip when present
    loh_bonus <- if (has_loh) n_chr * 15 else 0
    min(1800, max(400, n_chr * 120 + loh_bonus))
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
    results$fused_hap_labels     <- NULL
    results$loh_map              <- NULL
    results$loh_segments         <- NULL
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
      results$strain_ref        <- trimws(input$strain_ref)
      results$strain_alt        <- trimws(input$strain_alt)

      # Build LOH map using the beta-binomial EM + Viterbi HMM from chimera_functions.R.
      # full_read_loh applies only MAPQ + is_del filters (no base-quality filter)
      # so that allele-balance estimates reflect the deepest possible pileup.
      incProgress(0.87, detail = "Computing LOH map (HMM)")
      loh_out              <- compute_loh_map(
        res$full_read_loh,
        warn_fn = function(msg) showNotification(msg, type = "warning", duration = 10)
      )
      results$loh_map      <- loh_out$snp_table    # per-SNP; used by gap_has_loh()
      results$loh_segments <- loh_out$loh_segments # pre-collapsed runs; used by plots

      # Build a real per-position coverage map (same EM+HMM architecture as
      # the LOH map, but modeling total depth instead of allele balance).
      # Used by the chain analysis's terminal-event rules to tell a true
      # deletion (real depth drop) apart from a terminal LOH/crossover that
      # is simply missing its junction peak.
      incProgress(0.88, detail = "Computing coverage map (HMM)")
      coverage_out <- compute_coverage_map(
        res$full_read_loh,
        warn_fn = function(msg) showNotification(msg, type = "warning", duration = 10)
      )
      results$coverage_segments <- coverage_out$coverage_segments
      results$coverage_table    <- coverage_out$coverage_table

      incProgress(0.9, detail = "Creating individual peak plots")

      peak_chrs <- character(0)
      if (!is.null(res$snp_peaks) && nrow(res$snp_peaks) > 0) {
        peak_chrs <- unique(as.character(res$snp_peaks$chrom))
      }

      # ── Haplotype label lookup table (keyed by snp_peaks row index) ─────────
      # Tag snp_peaks with a stable row index before the loop so we can join
      # labels back even when some peaks are skipped (no touching reads).
      res$snp_peaks[, .row_idx := .I]

      hap_label_rows  <- list()

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

          # ── Haplotype classification (uses rt_df, may expand the window) ──
          hap <- classify_peak_haplotype(
            pk            = pk,
            chr_name      = chr_name,
            rt_df         = res$rt_df,
            touching_ids  = touching_ids,
            zone_min_snps = input$min_run
          )

          # Record label for writing back to snp_peaks later (keyed by .row_idx
          # so skipped peaks safely receive NA on the join)
          hap_label_rows[[length(hap_label_rows) + 1L]] <<- data.table(
            .row_idx         = pk$.row_idx,
            haplotype_label  = hap$label,
            hap_win_start    = hap$win_start,
            hap_win_end      = hap$win_end,
            hap_win_expanded = hap$expanded,
            n_read_support   = hap$n_support
          )

          # seg_data now comes from the classifier (window may be expanded)
          seg_data <- hap$seg_data

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

          x_lims <- range(plot_df$pos / 1000)

          # ── Build plot title with haplotype label ─────────────────────────
          hap_label_str <- gsub("_", " ", hap$label)
          expanded_note <- if (hap$expanded) " [window expanded]" else ""
          plot_title <- paste0(
            "Chr ", chr_name, "  \u2014  Peak ", pk_i,
            "  (SNP: ", round(snp_p / 1000, 2), " Kb;  Display: ",
            round(plot_start / 1000, 2), " \u2013 ",
            round(plot_end   / 1000, 2), " Kb)",
            "\n\u25b6 Haplotype classification: ", hap_label_str, expanded_note
          )

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
            # Original peak boundary lines (solid grey if not expanded,
            # or kept as reference even when the analysis window widened)
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
            ggtitle(plot_title)

          # If the analysis window was expanded, add orange dashed lines
          # showing the actual window used for classification
          if (hap$expanded) {
            p_reads <- p_reads +
              geom_vline(xintercept = hap$win_start / 1000,
                         color = "darkorange", linewidth = 0.7, linetype = 3) +
              geom_vline(xintercept = hap$win_end / 1000,
                         color = "darkorange", linewidth = 0.7, linetype = 3)
          }

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
              # Label the classification inside the segment bar
              annotate("text",
                       x     = mean(c(hap$win_start, hap$win_end)) / 1000,
                       y     = 0.5,
                       label = hap_label_str,
                       size  = 3,
                       color = "black",
                       fontface = "bold") +
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

          attr(p, "seg_data")       <- seg_data
          attr(p, "haplotype_label") <- hap$label
          p
        })

        plots <- Filter(Negate(is.null), plots)
        list(chromosome = chr_name, plots = plots)
      })

      # ── Write haplotype labels back to snp_peaks ─────────────────────────
      if (length(hap_label_rows) > 0) {
        hap_dt <- rbindlist(hap_label_rows, fill = TRUE)
        # .row_idx was added to res$snp_peaks before the loop; join on it.
        res$snp_peaks <- merge(
          res$snp_peaks, hap_dt,
          by = ".row_idx", all.x = TRUE
        )
      }
      res$snp_peaks[, .row_idx := NULL]
      results$snp_peaks <- res$snp_peaks

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
        loh_segments       = results$loh_segments,
        jaccard_threshold  = input$jaccard_threshold,
        zone_min_snps      = zone_min_snps,
        supervised_override = supervised_approved(),
        homog_frac         = input$homog_frac
      )

      results$peak_pairs  <- fusion_res$peak_pairs
      results$fused_peaks <- fusion_res$fused_peaks

      incProgress(0.7, detail = "Building fused peak plots")

      results$fused_peak_plots_by_chr <- build_fused_peak_plots(
        fused_peaks    = fusion_res$fused_peaks,
        rt_df          = results$rt_df,
        transition_pos = results$transition_pos,
        snp_coverage   = results$snp_coverage,
        zone_min_snps  = zone_min_snps
      )

      # Extract fused haplotype labels from plots and store for the blended summary table
      results$fused_hap_labels <- extract_fused_haplotype_labels(
        results$fused_peak_plots_by_chr,
        results$fused_peaks
      )
    })

    showNotification("Peak fusion complete. Review results in the 'Peak Fusion' and 'Fused Peak Plots' tabs.", type = "message", duration = 4)
  })
  
  # LOH chain events
  
  observeEvent(input$run_chain, {
    req(results$loh_segments, results$snp_peaks, results$chr_span)
    
    # Build params from UI inputs
    cp <- default_chain_params()
    cp$tel_tol_bp       <- as.integer(round(input$tel_tol_kb   * 1000))
    cp$merge_gap_bp     <- as.integer(round(input$merge_gap_kb * 1000))
    cp$min_span         <- as.integer(input$min_span_reads)
    cp$peak_pad_bp      <- as.integer(input$peak_pad_bp)
    cp$homog_frac       <- input$homog_frac
    cp$min_snps_for_peak <- as.integer(input$min_run)
    
    results$chain_params <- cp
    
    withProgress(message = "Running chain event caller...", value = 0.2, {
      
      chain_res <- run_chain_analysis(
        loh_segments      = results$loh_segments,
        fused_peaks       = results$fused_peaks, # NULL if fusion not run
        peak_pairs        = results$peak_pairs,
        snp_peaks         = results$snp_peaks,     # fallback when fused_peaks is NULL
        rt_df             = results$rt_df,
        chr_span          = results$chr_span,
        coverage_segments = results$coverage_segments,
        coverage_table    = results$coverage_table,
        params            = cp
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
        loh_segments        = results$loh_segments,
        jaccard_threshold   = input$jaccard_threshold,
        zone_min_snps       = zone_min_snps,
        supervised_override = supervised_approved(),
        homog_frac          = input$homog_frac
      )
      results$peak_pairs  <- fusion_res$peak_pairs
      results$fused_peaks <- fusion_res$fused_peaks

      incProgress(0.7, detail = "Rebuilding fused peak plots")

      results$fused_peak_plots_by_chr <- build_fused_peak_plots(
        fused_peaks    = fusion_res$fused_peaks,
        rt_df          = results$rt_df,
        transition_pos = results$transition_pos,
        snp_coverage   = results$snp_coverage,
        zone_min_snps  = zone_min_snps
      )

      results$fused_hap_labels <- extract_fused_haplotype_labels(
        results$fused_peak_plots_by_chr,
        results$fused_peaks
      )
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
      loh_g <- if ("loh_in_gap" %in% names(pairs)) isTRUE(pairs$loh_in_gap[[i]]) else FALSE

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
        tags$span(style = "min-width:55px;",
          if (loh_g)
            tags$span(style = "background:#6f42c1; color:white; padding:1px 5px; border-radius:3px; font-size:0.8em;", "LOH")
          else
            tags$span(style = "color:#aaa; font-size:0.8em;", "\u2014")),
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
        tags$span(style="min-width:55px;",  "LOH"),
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

  # ── Peak Summary table ───────────────────────────────────────────────────────
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
    out[, snp_n_str         := ifelse(is.na(snp_n), "", as.character(snp_n))]
    out[, chimeric_reads_str := ifelse(
      is.na(chimeric_reads_at_snp), "", as.character(chimeric_reads_at_snp)
    )]

    # Haplotype label columns (present only after analysis has run the plot loop)
    has_hap  <- "haplotype_label"  %in% names(out)
    has_wins <- "hap_win_start"    %in% names(out) && "hap_win_end" %in% names(out)

    base_cols <- list(
      Chromosome                = out$chrom,
      `Peak Position (Kb)`      = out$peak_pos_kb,
      `Qualifying SNP (Kb)`     = out$snp_pos_kb_str,
      `Raw Count at SNP`        = out$snp_n_str,
      `Peak Start (Kb)`         = out$peak_start_kb,
      `Peak End (Kb)`           = out$peak_end_kb,
      `Peak Height`             = out$peak_height,
      `Chimeric Reads at SNP`   = out$chimeric_reads_str
    )
    if (has_hap) {
      base_cols[["Haplotype Label"]] <-
        ifelse(is.na(out$haplotype_label), "\u2014",
               gsub("_", " ", out$haplotype_label))
    }
    if (has_hap && has_wins) {
      base_cols[["Hap Win Start (Kb)"]] <-
        ifelse(is.na(out$hap_win_start), "\u2014",
               as.character(round(out$hap_win_start / 1000, 2)))
      base_cols[["Hap Win End (Kb)"]] <-
        ifelse(is.na(out$hap_win_end), "\u2014",
               as.character(round(out$hap_win_end / 1000, 2)))
      base_cols[["Win Expanded"]] <-
        ifelse(is.na(out$hap_win_expanded), "\u2014",
               ifelse(out$hap_win_expanded, "yes", "no"))
    }

    display <- as.data.table(base_cols)
    setorder(display, Chromosome, `Peak Position (Kb)`)
    display
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

  # ── Fused Peak Plots tab (post-fusion) ───────────────────────────────────────
  output$fused_peak_plots_tabs <- renderUI({
    if (is.null(results$fused_peak_plots_by_chr) ||
        length(results$fused_peak_plots_by_chr) == 0) {
      return(p("Run 'Run Peak Fusion' first to generate fused peak plots."))
    }

    chr_tabs <- lapply(results$fused_peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots    <- chr_data$plots

      if (length(plots) == 0) return(NULL)

      plot_outputs <- lapply(seq_along(plots), function(i) {
        tagList(
          plotOutput(paste0("fused_peak_plot_", chr_name, "_", i), height = "600px"),
          fluidRow(
            column(3,
              downloadButton(
                paste0("dl_fused_peak_png_", chr_name, "_", i),
                "Download Plot Image (.png)",
                class = "btn-sm btn-default"
              )
            ),
            column(4,
              downloadButton(
                paste0("dl_fused_peak_rds_", chr_name, "_", i),
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
        h5(paste0("Chromosome ", chr_name, " — ", length(plots), " fused group(s)")),
        hr(),
        do.call(tagList, plot_outputs)
      )
    })

    chr_tabs <- Filter(Negate(is.null), chr_tabs)
    if (length(chr_tabs) == 0)
      return(p("No fused peak groups to display."))
    do.call(tabsetPanel, chr_tabs)
  })

  # Render fused peak plots with download handlers
  observe({
    req(results$fused_peak_plots_by_chr)

    lapply(results$fused_peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots    <- chr_data$plots

      lapply(seq_along(plots), function(i) {
        output_name <- paste0("fused_peak_plot_", chr_name, "_", i)
        png_dl_id   <- paste0("dl_fused_peak_png_", chr_name, "_", i)
        rds_dl_id   <- paste0("dl_fused_peak_rds_", chr_name, "_", i)
        local({
          p    <- plots[[i]]
          .chr <- chr_name
          .i   <- i

          output[[output_name]] <- renderPlot({ p })

          output[[png_dl_id]] <- downloadHandler(
            filename = function() {
              paste0(input$sample_name, "_fused_peak_chr", .chr,
                     "_group", .i, "_", Sys.Date(), ".png")
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
              paste0(input$sample_name, "_fused_peak_chr", .chr,
                     "_group", .i, "_", Sys.Date(), ".rds")
            },
            content = function(file) {
              reads_plot <- if (inherits(p, "patchwork")) p[[1]] else p
              plot_data  <- as.data.table(reads_plot$data)
              plot_data  <- plot_data[, .(chrom, pos, read_id, IS_REF)]
              seg_dt     <- attr(p, "seg_data")
              saveRDS(
                list(
                  plot_data        = plot_data,
                  seg_data         = seg_dt,
                  haplotype_label  = attr(p, "haplotype_label"),
                  chromosome       = .chr,
                  fused_group      = .i,
                  sample_name      = input$sample_name,
                  app_version      = APP_VERSION
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
    # Read snp_peaks and chr_span as reactive dependencies so the observer
    # re-runs when they are set (chr_span is assigned after snp_coverage in the
    # analysis observer, so an explicit read here ensures .x_max_kb is
    # recomputed from chromosome lengths rather than falling back to pos_kb).
    snp_peaks_all <- if (!is.null(results$snp_peaks) && nrow(results$snp_peaks) > 0)
      copy(results$snp_peaks)[, chrom := as.character(chrom)]
    else NULL
    chr_span_local <- results$chr_span   # explicit reactive dep; may be NULL on first fire
    snp_cov_all[, chrom := as.character(chrom)]
    fits_all[,    chrom := as.character(chrom)]
    if (!is.null(peaks_all) && nrow(peaks_all) > 0)
      peaks_all <- copy(peaks_all)[, chrom := as.character(chrom)]

    # Per-chromosome x-axis: build a lookup of max Kbp per chromosome from
    # chr_span (chromosome length) so each panel is scaled to its own extent
    # rather than the genome-wide maximum.  Falls back to observed SNP positions
    # when chr_span is unavailable.
    chr_x_max_kb <- if (!is.null(chr_span_local) && nrow(chr_span_local) > 0) {
      span_dt <- copy(chr_span_local)
      span_dt[, chrom := as.character(chrom)]
      setNames(
        ceiling(ceiling(span_dt$length / 1000) / 25) * 25,
        span_dt$chrom
      )
    } else {
      # Fallback: max observed SNP position per chromosome
      snp_cov_all[, .(x_max = ceiling(ceiling(max(pos_kb)) / 25) * 25), by = chrom] |>
        (\(dt) setNames(dt$x_max, dt$chrom))()
    }

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
        # Per-chromosome x-axis upper limit (Kbp), rounded to nearest 100 Kb
        .x_max_kb      <- chr_x_max_kb[.chr] %||%
                          ceiling(ceiling(max(snp_cov_all[chrom == .chr, pos_kb])) / 25) * 25
        .snp           <- snp_cov_all[chrom == .chr]
        .fits          <- fits_all[chrom == .chr]
        .peaks         <- if (!is.null(peaks_all) && nrow(peaks_all) > 0)
                            peaks_all[chrom == .chr] else NULL
        .snp_peaks_chr <- if (!is.null(snp_peaks_all) && nrow(snp_peaks_all) > 0)
                            snp_peaks_all[chrom == .chr]
                          else NULL
        # LOH data is read inside renderPlot (not captured here) so that it
        # picks up results$loh_map reactively after analysis completes, and
        # so that copy() prevents in-place pos mutation across renders.

        output[[.plot_id]] <- renderPlot({
          y_ceil <- max(30, max(.snp$n))

          # Compute x-axis limit reactively inside renderPlot so it always
          # reflects the live chr_span value rather than a value captured
          # at observer-construction time (which may be NULL / stale).
          x_max_kb_live <- if (!is.null(results$chr_span) && nrow(results$chr_span) > 0) {
            span_row <- results$chr_span[as.character(chrom) == .chr]
            if (nrow(span_row) > 0)
              ceiling(ceiling(span_row$length / 1000) / 25) * 25
            else
              ceiling(ceiling(max(.snp$pos_kb)) / 25) * 25
          } else {
            ceiling(ceiling(max(.snp$pos_kb)) / 25) * 25
          }

          # Read LOH segments fresh each render — use pre-collapsed table
          # from compute_loh_map()$loh_segments so no inline RLE is needed.
          .loh_chr <- if (!is.null(results$loh_segments) && nrow(results$loh_segments) > 0) {
            segs <- copy(results$loh_segments)
            segs[, chrom := as.character(chrom)]
            segs[chrom == .chr & loh_state %in% c("REF_fixed", "ALT_fixed")]
          } else NULL

          # Read events fresh each render — set after the chain event caller
          # (input$run_chain) completes; NULL beforehand.
          .events_chr <- results$event_table

          p <- ggplot(.snp, aes(x = pos_kb, y = n)) +
            geom_line(
              data  = .fits,
              aes(x = uniform_pos / 1000, y = uniform_fit),
              color = "firebrick", linewidth = 0.8, alpha = 0.7
            ) +
            geom_point(color = "black", alpha = 0.5, size = 0.8, shape = 21) +
            scale_x_continuous(limits = c(0, x_max_kb_live), minor_breaks = seq(0, x_max_kb_live, 100)) +
            xlab("Position (Kbp)") +
            ylab("Number of Reads") +
            ylim(NA, y_ceil) +
            ggtitle(paste("Chromosome", .chr)) +
            theme_bw() +
            theme(
              panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
              panel.grid.major.x = element_line(linewidth = 0.05, color = "red"),
              axis.text          = element_text(size = rel(1.2)),
              axis.title         = element_text(size = rel(1.5))
              
            )

          # ── Peak highlight points ──────────────────────────────────────────
          if (!is.null(.snp_peaks_chr) && nrow(.snp_peaks_chr) > 0) {
            peak_highlight <- merge(
              .snp_peaks_chr[!is.na(snp_pos), .(pos = snp_pos)],
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

          # ── LOH band at the base of the plot ──────────────────────────────
          # Drawn as a thin geom_rect strip (4 % of y ceiling) flush with the
          # x-axis baseline.  REF-fixed = dodgerblue, ALT-fixed = firebrick.
          # Uses pre-collapsed loh_segments — no inline rleid/half_step needed.
          if (!is.null(.loh_chr) && nrow(.loh_chr) > 0) {
            .loh_chr[, xmin := start / 1000]
            .loh_chr[, xmax := end   / 1000]

            loh_band_h  <- y_ceil * -0.10 # negative to put under number line
            loh_colours <- c(REF_fixed = "dodgerblue", ALT_fixed = "firebrick")

            # Resolve strain display names from results (fall back to generic)
            s_ref <- if (!is.null(results$strain_ref) && nzchar(results$strain_ref))
              results$strain_ref else "REF"
            s_alt <- if (!is.null(results$strain_alt) && nzchar(results$strain_alt))
              results$strain_alt else "ALT"

            loh_labels  <- c(
              REF_fixed = s_ref,
              ALT_fixed = s_alt
            )
            loh_caption <- paste0(
              "LOH band (bottom): blue\u202f=\u202f", s_ref,
              ", red\u202f=\u202f", s_alt
            )

            p <- p +
              geom_rect(
                data        = .loh_chr,
                aes(xmin = xmin, xmax = xmax,
                    ymin = 0,    ymax = loh_band_h,
                    fill = loh_state),
                alpha       = 0.85,
                inherit.aes = FALSE
              ) +
              scale_fill_manual(
                values = loh_colours,
                labels = loh_labels,
                name   = "LOH",
                drop   = FALSE
              ) +
              labs(caption = loh_caption) +
              theme(
                plot.caption    = element_blank(), #element_text(size = 7, colour = "grey40"),
                legend.position = "right",
                legend.title    = element_text(size = rel(1.1), face = "bold"),
                legend.text     = element_text(size = rel(1.1))
              )

            p <- add_event_symbols(p, .events_chr, band_ymin = 0, band_ymax = loh_band_h,
                                    chrom_filter = .chr)
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

  # ── Post Fusion Peak Summary table ──────────────────────────────────────────
  # Builds the blended summary: one row per initial peak, with the post-fusion
  # haplotype classification substituted in for fused groups.
  # Singletons retain their original classification from snp_peaks.

  build_post_fusion_summary <- reactive({
    req(results$snp_peaks)

    peaks <- copy(results$snp_peaks)
    peaks[, chrom := as.character(chrom)]

    # Start with the initial (pre-fusion) haplotype label
    has_hap <- "haplotype_label" %in% names(peaks)

    out <- peaks[, {
      list(
        Chromosome            = chrom,
        `Peak Position (Kb)`  = round(peak_pos   / 1000, 2),
        `Qualifying SNP (Kb)` = ifelse(is.na(snp_pos), "None above cutoff",
                                       as.character(round(snp_pos / 1000, 2))),
        `Peak Start (Kb)`     = round(peak_start / 1000, 2),
        `Peak End (Kb)`       = round(peak_end   / 1000, 2),
        `Peak Height`         = round(peak_height, 2),
        `Chimeric Reads`      = ifelse(is.na(chimeric_reads_at_snp), "",
                                       as.character(chimeric_reads_at_snp)),
        `Initial Classification` = if (has_hap)
          ifelse(is.na(haplotype_label), "\u2014", gsub("_", " ", haplotype_label))
        else "\u2014",
        # Placeholders filled below
        `Fusion Status`          = "singleton",
        `Fusion Group`           = NA_integer_,
        `Sub-peaks in Group`     = 1L,
        `Post-Fusion Classification` = if (has_hap)
          ifelse(is.na(haplotype_label), "\u2014", gsub("_", " ", haplotype_label))
        else "\u2014",
        `Classification Changed` = "no"
      )
    }]

    # If fusion has been run, overlay fused group information
    if (!is.null(results$fused_peaks) && !is.null(results$fused_hap_labels)) {

      fp  <- copy(results$fused_peaks)
      fp[, chrom := as.character(chrom)]
      fhl <- copy(results$fused_hap_labels)

      # fused_peaks has one row per original peak_id with fusion_group membership
      # Join the fused haplotype label from fhl via fusion_group
      fhl_keyed <- fhl[, .(fusion_group, fused_classification)]

      # Build a peak_id -> fusion_group + fused_classification lookup
      pk_lookup <- fp[, .(peak_id, chrom, snp_pos, fusion_group, n_sub_peaks,
                          best_edge_type, best_fusion_mode)]
      pk_lookup <- merge(pk_lookup, fhl_keyed, by = "fusion_group", all.x = TRUE)

      # We need to match by chrom + snp_pos (bp) since out doesn't have peak_id.
      # snp_peaks and fused_peaks share the same rows/snp_pos values.
      # Re-attach peak_id to out via snp_pos + chrom
      peaks_id <- peaks[, .(chrom, snp_pos, peak_start, peak_end)]
      peaks_id[, row_seq := .I]

      out[, row_seq := .I]

      # Join pk_lookup onto out row by row using chrom + snp_pos
      pk_lookup[, snp_pos_chr := paste(chrom, snp_pos)]
      out_join_key <- paste(
        peaks[, chrom],
        peaks[, snp_pos]
      )

      for (ri in seq_len(nrow(out))) {
        key <- out_join_key[ri]
        match_rows <- pk_lookup[paste(chrom, snp_pos) == key]
        if (nrow(match_rows) == 0) next

        mr <- match_rows[1]
        fg <- mr$fusion_group
        ns <- mr$n_sub_peaks
        et <- mr$best_edge_type
        fm <- mr$best_fusion_mode
        fc <- mr$fused_classification

        status <- if (!is.na(fm) && fm == "automatic") {
          if (ns > 1L) "fused (auto)" else "singleton"
        } else if (!is.na(fm) && fm == "supervised") {
          if (ns > 1L) "fused (supervised)" else "singleton"
        } else {
          "singleton"
        }

        fused_label <- if (!is.na(fc) && fc != "") {
          gsub("_", " ", fc)
        } else {
          out$`Initial Classification`[ri]
        }

        changed <- if (fused_label != out$`Initial Classification`[ri] &&
                       !is.na(fused_label) && fused_label != "\u2014") "yes" else "no"

        out$`Fusion Status`[ri]              <- status
        out$`Fusion Group`[ri]               <- fg
        out$`Sub-peaks in Group`[ri]         <- ns
        out$`Post-Fusion Classification`[ri] <- fused_label
        out$`Classification Changed`[ri]     <- changed
      }
      out[, row_seq := NULL]
    }

    setorder(out, Chromosome, `Peak Position (Kb)`)
    out
  })

  output$post_fusion_peak_summary_table <- renderTable({
    req(results$snp_peaks)
    df <- build_post_fusion_summary()

    # If fusion not yet run, hide the fusion-specific columns
    if (is.null(results$fused_peaks)) {
      cols_keep <- c("Chromosome", "Peak Position (Kb)", "Qualifying SNP (Kb)",
                     "Peak Start (Kb)", "Peak End (Kb)", "Peak Height",
                     "Chimeric Reads", "Initial Classification")
      df[, ..cols_keep]
    } else {
      df
    }
  }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s",
     na = "\u2014")

  # Legend for Post Fusion Peak Summary
  output$post_fusion_summary_legend <- renderUI({
    if (is.null(results$fused_peaks)) {
      return(helpText("Run 'Run Peak Fusion' to populate fusion status and post-fusion classifications."))
    }
    tags$div(
      style = "margin-bottom: 8px; font-size: 0.85em; color: #444;",
      strong("Classification types: "),
      tags$span("binary", style = "background:#cce5ff; padding:2px 7px; border-radius:3px; margin-right:5px;"),
      tags$span("gene conversion", style = "background:#d4edda; padding:2px 7px; border-radius:3px; margin-right:5px;"),
      tags$span("internal crossover", style = "background:#fff3cd; padding:2px 7px; border-radius:3px; margin-right:5px;"),
      tags$span("undefined", style = "background:#e2e3e5; padding:2px 7px; border-radius:3px; margin-right:5px;"),
      br(), br(),
      strong("Fusion status: "),
      tags$span("singleton", style = "background:#f8f9fa; padding:2px 7px; border-radius:3px; margin-right:5px; border:1px solid #ccc;"),
      tags$span("fused (auto)", style = "background:#d4edda; padding:2px 7px; border-radius:3px; margin-right:5px;"),
      tags$span("fused (supervised)", style = "background:#fff3cd; padding:2px 7px; border-radius:3px; margin-right:5px;"),
      br(), br(),
      tags$em("'Classification Changed' = yes indicates the post-fusion re-classification differs from the initial per-peak label.")
    )
  })

  output$download_post_fusion_summary <- downloadHandler(
    filename = function() paste0(input$sample_name, "_post_fusion_peak_summary_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(results$snp_peaks)
      fwrite(build_post_fusion_summary(), file)
    }
  )

  # ── LOH Regions table ────────────────────────────────────────────────────────
  output$loh_regions_table <- renderTable({
    req(results$loh_segments)

    segs <- copy(results$loh_segments)
    segs <- segs[loh_state %in% c("REF_fixed", "ALT_fixed")]

    if (nrow(segs) == 0)
      return(data.frame(Message = "No fixed LOH regions detected."))

    segs[, chrom := as.character(chrom)]

    # Resolve strain display names
    s_ref <- if (!is.null(results$strain_ref) && nzchar(results$strain_ref))
      results$strain_ref else "REF"
    s_alt <- if (!is.null(results$strain_alt) && nzchar(results$strain_alt))
      results$strain_alt else "ALT"

    segs[, `LOH State` := fcase(
      loh_state == "REF_fixed", paste0("REF fixed (", s_ref, ")"),
      loh_state == "ALT_fixed", paste0("ALT fixed (", s_alt, ")")
    )]

    out <- segs[, .(
      Chromosome          = chrom,
      `Start (bp)`        = format(start,      big.mark = ",", scientific = FALSE),
      `End (bp)`          = format(end,         big.mark = ",", scientific = FALSE),
      `Length (kb)`       = round(length_bp / 1000, 2),
      `SNPs in Region`    = n_snps,
      `LOH State`         = `LOH State`,
      `Mean Allele Balance` = round(balance_mean, 3),
      `SD Allele Balance`   = ifelse(is.na(balance_sd), "\u2014",
                                     as.character(round(balance_sd, 3)))
    )]

    setorder(out, Chromosome, `Start (bp)`)
    out
  }, striped = TRUE, hover = TRUE, bordered = TRUE, spacing = "s", na = "\u2014")

  output$download_loh_regions <- downloadHandler(
    filename = function() paste0(input$sample_name, "_loh_regions_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(results$loh_segments)

      segs <- copy(results$loh_segments)
      segs <- segs[loh_state %in% c("REF_fixed", "ALT_fixed")]
      segs[, chrom := as.character(chrom)]

      s_ref <- if (!is.null(results$strain_ref) && nzchar(results$strain_ref))
        results$strain_ref else "REF"
      s_alt <- if (!is.null(results$strain_alt) && nzchar(results$strain_alt))
        results$strain_alt else "ALT"

      segs[, loh_state_label := fcase(
        loh_state == "REF_fixed", paste0("REF_fixed (", s_ref, ")"),
        loh_state == "ALT_fixed", paste0("ALT_fixed (", s_alt, ")")
      )]

      out <- segs[, .(
        sample      = input$sample_name,
        chrom,
        start,
        end,
        length_bp,
        length_kb   = round(length_bp / 1000, 2),
        n_snps,
        loh_state   = loh_state_label,
        balance_mean = round(balance_mean, 4),
        balance_sd   = round(balance_sd,   4)
      )]

      setorder(out, chrom, start)
      fwrite(out, file)
    }
  )
  
  # --- Event table output
  
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

  # ── Download handlers ────────────────────────────────────────────────────────

  output$download_plot <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chromosome_tracking_", Sys.Date(), ".png"),
    content  = function(file) {
      req(results$snp_coverage, results$chromosome_fits)
      p       <- build_overview_plot(results)
      n_chr   <- length(unique(results$snp_coverage$chrom))
      has_loh <- !is.null(results$loh_segments) &&
                 nrow(results$loh_segments) > 0 &&
                 any(results$loh_segments$loh_state %in% c("REF_fixed", "ALT_fixed"), na.rm = TRUE)
      loh_bonus <- if (has_loh) n_chr * 0.12 else 0
      png_h   <- max(3, min(18, n_chr * 1.2 + loh_bonus))
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
