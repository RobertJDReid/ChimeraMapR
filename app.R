library(shiny)
source("chimera_functions.R")   # loads all packages + whittaker + run_chimera_analysis etc.
# APP_VERSION is now defined inside chimera_functions.R

# ─────────────────────────────────────────────
#                   UI
# ─────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    .well .form-group              { margin-top: 2px; margin-bottom: 2px; }
    .well .shiny-input-container   { margin-top: 2px; margin-bottom: 2px; }
    .well hr                       { margin-top: 2px; margin-bottom: 2px; }
    .well h3                       { margin-top: 2px; margin-bottom: 2px; }
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

#      hr(),

      h3("Analysis Parameters"),
      textInput("sample_name",
                "Sample Name:",
                value = "Sample_01"),

      numericInput("mapq_cutoff",
                   "Minimum MAPQ Value:",
                   value = 20,
                   min = 0,
                   step = 1),
#      helpText("Reads below this score are excluded"),

      numericInput("baseq_cutoff",
                   "Base Quality Minimum:",
                   value = 10,
                   min = 0,
                   step = 1),
#      helpText("Minimum base quality score at SNP"),

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
#      helpText("Set to ~1/2 of median read depth"),

      numericInput("lambda",
                   "Whittaker Lambda (λ):",
                   value = 1,
                   min   = 0.01,
                   step  = 0.5),
      helpText("Smoothness penalty for Whittaker smoother. Lower = tighter fit (preserves sharp peaks); higher = smoother curve, lower = tighter"),

#      hr(),

      actionButton("run_analysis",
                   "Run Analysis",
                   class = "btn-primary btn-lg"),
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
                 helpText("Shows the lambda (λ) value used for each chromosome in the analysis"),
                 br(),
                 h5("Export Whittaker fitted curves"),
                 helpText("Downloads the fitted Whittaker curves for all chromosomes from this run, including run parameters for later comparison across runs."),
                 downloadButton("download_curve_fits", "Download Curve Fits")
        ),
        
        # NEW: dynamic/closeable Selected Region tab placeholder
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
                   tags$li("Applies Whittaker smoothing to identify peaks of per-read haplotype switches")
                 ),
                 h5("Input Files:"),
                 tags$ul(
                   tags$li(strong("Read Data:"), "CSV file with SNP position information from BAM file (columns: chrom, pos, read_id, call, is_del, etc.). For csv files > 200 Mb, compress with ",em("gzip")," or ",em("pigz")," prior to upload."),
                   tags$li(strong("SNP Data:"), "CSV file with SNP positions (columns: CHROM, POS, REF, ALT), or a VCF file (plain or gzipped). Multi-allelic VCF sites are split into one row per ALT allele."),
                   tags$li(strong("Chromosome Size:"), "FAI index file with chromosome lengths")
                 ),
                 h5("Whittaker Lambda (λ):"),
                 tags$ul(
                   tags$li(strong("Low λ (0.01–1):"), "Tight fit — preserves sharp, narrow peaks well. Recommended for impulsive boundary-count signals."),
                   tags$li(strong("High λ (10–1000):"), "Heavy smoothing — useful for broad signal or very noisy data, but will attenuate sharp peaks.")
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

  # Shared trigger for region plot building — incremented by both the overview
  # button and the per-chromosome buttons
  rv_trigger_region <- reactiveVal(0L)

  # Reactive values to store analysis results
  results <- reactiveValues(
    rt_df                 = NULL,
    transition_pos        = NULL,   # boundary positions only (new method)
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
    selected_region_data  = NULL
  )

  # ── RLE helper (operates on plain vectors) ──────────────────────
  rle_helper <- function(x) {
    r  <- rle(x)[[1]]   # run lengths
    rn <- rep(r, r)     # expand to per-position run lengths
    return(rn)
  }

  # ── Main analysis ────────────────────────────────────────────────────────────
  observeEvent(input$run_analysis, {

    req(input$read_data_file, input$snp_data_file, input$chr_size_file)

    # Reset region-selection outputs on re-run
    results$selected_region      <- NULL
    results$selected_region_plot <- NULL
    results$selected_region_data <- NULL

    # Remove Selected Region tab if present
    try(removeTab(inputId = "main_tabs", target = "Selected Region"), silent = TRUE)

    withProgress(message = "Processing data...", value = 0, {
      
      incProgress(0.1, detail = "Running analysis")

      # Shiny strips extensions from uploaded file temp paths, which breaks
      # format detection in load_snp_data (VCF vs CSV) and other loaders.
      # Restore the original suffix (handles double extensions like .vcf.gz).
      restore_ext <- function(datapath, original_name) {
        suffix <- sub("^[^.]+", "", original_name)   # e.g. ".vcf.gz", ".csv", ".fai"
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
      
      # Peak plots are built here — they depend on Shiny output IDs and
      # are not worth moving into chimera_functions.R.
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

          # Reads whose transition boundary falls inside the peak window
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

          # Layer 2: explicit data frame so the RDS download handler can access it
          peak_point_df <- data.frame(
            pos_kb     = snp_p / 1000,
            peak_start = pk_start / 1000,
            peak_end   = pk_end   / 1000,
            peak_height = pk$peak_height
          )

          # ── Haplotype segment data (classify2 method) ──────────────────────
          # Summarise allele balance per position within the peak window,
          # then run-length-encode into contiguous REF/HET/ALT blocks.
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
              xmin    = min(pos) / 1000,
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

          # ── Haplotype segment panel ────────────────────────────────────────
          if (!is.null(seg_data) && nrow(seg_data) > 0) {
            p_seg <- ggplot(seg_data) +
              geom_rect(
                aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
                    fill = SNP_call),
                alpha = 0.45
              ) +
              scale_fill_manual(
                values  = c(ALT = "firebrick", HET = "gray60", REF = "dodgerblue"),
                drop    = FALSE,
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
                panel.grid    = element_blank(),
                axis.text.y   = element_blank(),
                axis.ticks.y  = element_blank(),
                panel.border  = element_blank(),
                legend.position = "none"
              )
            p <- patchwork::wrap_plots(p_reads, p_seg, ncol = 1,
                                       heights = c(8, 1))
          } else {
            p <- p_reads + xlab("Position (Kbp)")
          }

          # Attach seg_data as an attribute so the RDS handler can retrieve it
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

    showNotification("Analysis complete", type = "message", duration = 3)
  })


  # ── Debug outputs ────────────────────────────────────────────────────────────
  # removed

  # ── Outputs ──────────────────────────────────────────────────────────────────

  output$selected_region_text <- renderText({
    reg <- results$selected_region
    if (is.null(reg) || is.null(reg$chrom)) {
      return("No region selected.")
    }

    paste0(
      "Selected region: Chr ", reg$chrom,
      " : ", format(reg$start, big.mark = ","),
      " - ", format(reg$end, big.mark = ","),
      " bp (",
      round((reg$end - reg$start) / 1000, 2),
      " kb; fixed max window)"
    )
  })

  # Overview plot — built fresh inside renderPlot so Shiny correctly registers
  # the coordinate domain and panel mappings needed for brush to work
  output$chr_plot <- renderPlot({
    req(results$snp_coverage, results$chromosome_fits)
    build_overview_plot(results)
  }, height = function() {
    req(results$snp_coverage)
    n_chr <- length(unique(results$snp_coverage$chrom))
    min(1600, max(400, n_chr * 120))
  })

  # Peaks table
  output$peaks_table <- renderTable({
    req(results$peaks_genomic)
    out <- copy(results$snp_peaks)
    out[, `:=`(
      peak_pos_kb   = round(peak_pos   / 1000, 2),
      peak_start_kb = round(peak_start / 1000, 2),
      peak_end_kb   = round(peak_end   / 1000, 2),
      peak_height   = round(peak_height, 2)
    )]
    # Format SNP position and read count as character so NA can display a label
    out[, snp_pos_kb_str := ifelse(
      is.na(snp_pos),
      "None above cutoff",
      as.character(round(snp_pos / 1000, 2))
    )]
    out[, snp_n_str := ifelse(
      is.na(snp_n),
      "",
      as.character(snp_n)
    )]
    out[, chimeric_reads_str := ifelse(
      is.na(chimeric_reads_at_snp),
      "",
      as.character(chimeric_reads_at_snp)
    )]
    out <- out[, .(
      Chromosome                     = chrom,
      `Peak Position (Kb)` = peak_pos_kb,
      `Qualifying SNP (Kb)`          = snp_pos_kb_str,
      `Raw Count at SNP`             = snp_n_str,
      `Peak Start (Kb)`              = peak_start_kb,
      `Peak End (Kb)`                = peak_end_kb,
      `Peak Height`        = peak_height,
      `Chimeric Reads at SNP`        = chimeric_reads_str
    )]
    setorder(out, Chromosome, `Peak Position (Kb)`)
    out
  })

  # Dynamic UI: per-chromosome tabs of individual peak plots
  output$peak_plots_tabs <- renderUI({
    req(results$peak_plots_by_chr)

    if (is.null(results$peak_plots_by_chr) || length(results$peak_plots_by_chr) == 0) {
      return(h4("No peaks detected in the analysis."))
    }

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

  # Render individual peak plots dynamically, with PNG and RDS download handlers
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
              # For patchwork objects, extract read count from the first ggplot
              reads_plot <- if (inherits(p, "patchwork")) p[[1]] else p
              n_reads    <- length(unique(ggplot_build(reads_plot)$data[[1]]$group))
              # Add extra height for the haplotype segment panel
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
              # For patchwork objects the read data lives in the first sub-plot;
              # for plain ggplot objects it lives directly in p$data.
              reads_plot <- if (inherits(p, "patchwork")) p[[1]] else p

              plot_data <- as.data.table(reads_plot$data)

              plot_data <- plot_data[, .(
                chrom,
                pos,
                read_id,
                IS_REF
              )]

              # peak_points layer is always layer 2 of the reads sub-plot
              peak_points_dt <- as.data.table(reads_plot$layers[[2]]$data)

              # Retrieve seg_data stored as an attribute during plot construction
              seg_dt <- attr(p, "seg_data")

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

  # ── Chromosome plots tab: dynamic subtabs ────────────────────────────────────
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
          column(3,
            actionButton(btn_id, "Plot Selected Region", class = "btn-primary")
          ),
          column(9,
            verbatimTextOutput(txt_id)
          )
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
    snp_cov_all[,  chrom := as.character(chrom)]
    fits_all[,     chrom := as.character(chrom)]
    if (!is.null(peaks_all) && nrow(peaks_all) > 0)
      peaks_all <- copy(peaks_all)[, chrom := as.character(chrom)]
    
    lapply(chr_levels, function(chr_name) {
      chr_c    <- as.character(chr_name)
      plot_id  <- paste0("chr_cov_plot_", chr_c)
      brush_id <- paste0("chr_cov_brush_", chr_c)
      btn_id   <- paste0("chr_plot_btn_",  chr_c)
      txt_id   <- paste0("chr_region_text_", chr_c)
      local({
        .chr      <- chr_c
        .plot_id  <- plot_id
        .brush_id <- brush_id
        .btn_id   <- btn_id
        .txt_id   <- txt_id
        .snp      <- snp_cov_all[chrom == .chr]
        .fits     <- fits_all[chrom == .chr]
        .peaks    <- if (!is.null(peaks_all) && nrow(peaks_all) > 0)
                       peaks_all[chrom == .chr] else NULL
        .snp_peaks_chr <- if (!is.null(snp_peaks_all) && nrow(snp_peaks_all) > 0)
          snp_peaks_all[chrom == .chr]
        else NULL
      
        # Render the coverage plot for this chromosome
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

#          if (!is.null(.peaks) && nrow(.peaks) > 0)
#            p <- p + geom_vline(aes(xintercept = peak_pos / 1000),
#                                data = .peaks, color = "lightgreen", alpha = 0.5)
          p
        })
      
        # Capture brush coordinates for this chromosome
        observeEvent(input[[.brush_id]], {
          brush <- input[[.brush_id]]
          if (is.null(brush)) return(NULL)
      
          # Guard against pixel-space coords returned before plot fully registers
          chr_max_kb <- max(.snp$pos_kb, na.rm = TRUE)
          if (brush$xmin > chr_max_kb * 2 || brush$xmax > chr_max_kb * 10) return(NULL)
      
          xmin_kb <- min(brush$xmin, brush$xmax)
          xmax_kb <- max(brush$xmin, brush$xmax)
          xmin_bp <- as.integer(round(xmin_kb * 1000))
          xmax_bp <- as.integer(round(xmax_kb * 1000))
      
          # Store per-chromosome so brushes don't overwrite each other
          results$selected_regions[[.chr]] <- list(
            chrom = .chr,
            start = xmin_bp,
            end   = xmax_bp
          )
        }, ignoreNULL = TRUE)
      
        # Region label for this chromosome
        output[[.txt_id]] <- renderText({
          reg <- results$selected_regions[[.chr]]
          if (is.null(reg)) return("No region selected.")
          paste0(
            "Selected: ", format(reg$start, big.mark = ","),
            " - ", format(reg$end, big.mark = ","), " bp  (",
            round((reg$end - reg$start) / 1000, 2), " kb)"
          )
        })
      
        # Button: push this chromosome's brush region to the shared active region,
        # then increment the trigger to fire build_region_plot()
        observeEvent(input[[.btn_id]], {
          reg <- results$selected_regions[[.chr]]
          req(!is.null(reg))
          results$selected_region <- reg
          rv_trigger_region(isolate(rv_trigger_region()) + 1L)
        }, ignoreNULL = TRUE)
      })        # closes local()
    })          # closes lapply()
  })            # closes observe()

  # Whittaker lambda table
  output$span_table <- renderTable({
    req(results$chr_span)
    out <- copy(results$chr_span)
    out[, length_kb := round(length / 1000, 1)]
    out <- out[, .(
      Chromosome    = chrom,
      `Length (Kb)` = length_kb,
      `SNP Count`   = snp_count,
      `Lambda (λ)`  = lambda
    )]
    setorder(out, Chromosome)
    out
  })

  # Summary statistics
  output$summary_stats <- renderText({
    req(results$rt_df, results$snp_coverage, results$peaks_genomic)

    n_chimeric_reads   <- uniqueN(results$rt_df$read_id)
    n_chromosomes      <- uniqueN(results$snp_coverage$chrom)
    n_peaks            <- nrow(results$peaks_genomic)
    total_snps         <- nrow(results$snp_coverage)
    boundary_snps      <- results$snp_coverage[n > 0, .N]
    total_boundaries   <- sum(results$snp_coverage$n)

    paste0(
      "Sample: ", input$sample_name, "\n\n",
      "Chimeric Reads Detected: ", n_chimeric_reads, "\n",
      "Chromosomes Analyzed: ", n_chromosomes, "\n",
      "Total SNP Positions: ", total_snps, "\n",
      "SNP Positions with Transition Boundaries: ", boundary_snps,
      " (", round(100 * boundary_snps / total_snps, 1), "%)\n",
      "Total Transition Boundary Events: ", total_boundaries, "\n",
      "Peaks Detected: ", n_peaks, "\n\n",
      "Analysis Parameters:\n",
      "  MAPQ Cutoff: ",              input$mapq_cutoff,     "\n",
      "  Base Quality Cutoff: ",      input$baseq_cutoff,    "\n",
      "  Min Run Length: ",           input$min_run,         "\n",
      "  Min Peak Height: ",          input$min_peak_height, "\n",
      "  Whittaker Lambda (λ): ",     input$lambda
    )
  })

  # ── Shared helper: build and display the selected-region read plot ───────────
  build_region_plot <- function() {
    req(results$selected_region, results$rt_df, results$transition_pos)

    reg <- results$selected_region

    # Select reads that contributed a transition boundary inside the region;
    # then plot their full span from rt_df for context.
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
            column(
              3,
              downloadButton(
                "download_selected_region_plot",
                "Download Plot Image (.png)"
              )
            ),
            column(
              4,
              downloadButton(
                "download_selected_region_rds",
                "Download Plot Data (.rds)"
              )
            ),
            column(
              3,
              actionButton(
                "close_selected_region_tab",
                "Close Tab",
                class = "btn-secondary"
              )
            )
          )
        )
      )
    } else {
      updateTabsetPanel(session, "main_tabs", selected = "Selected Region")
    }
  }

  # ── Selected-region plot generation ─────────────────────────────────────────

  # Triggered by any per-chromosome "Plot Selected Region" button
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

  # ── Overview plot: download as R object for replotting ───────────────────────
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
    #content  = function(file) fwrite(results$peaks_genomic, file)
    content  = function(file) fwrite(results$snp_peaks, file)
  )

  output$download_read_ids <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chimeric_read_ids_", Sys.Date(), ".txt"),
    content  = function(file) writeLines(results$chimeric_read_ids, file)
  )

  # selected region plot image
  output$download_selected_region_plot <- downloadHandler(
    filename = function() {
      req(results$selected_region)
      reg <- results$selected_region
      paste0(
        input$sample_name,
        "_selected_region_chr", reg$chrom, "_",
        reg$start, "_", reg$end, "_", Sys.Date(), ".png"
      )
    },
    content = function(file) {
      req(results$selected_region_plot)
      ggsave(file, plot = results$selected_region_plot, width = 10, height = 12, dpi = 300)
    }
  )
  # selected region data
  output$download_selected_region_rds <- downloadHandler(
    filename = function() {
      req(results$selected_region)
      reg <- results$selected_region
      
      paste0(
        input$sample_name,
        "_selected_region_chr", reg$chrom, "_",
        reg$start, "_", reg$end, "_",
        Sys.Date(), ".rds"
      )
    },
    content = function(file) {
      req(results$selected_region_data, results$selected_region)
      
      reg <- results$selected_region
      
      plot_data <- as.data.table(results$selected_region_data)
      
      plot_data <- plot_data[, .(
        chrom,
        pos,
        read_id,
        IS_REF,
        ALLELE
      )]
      
      saveRDS(
        list(
          plot_data = plot_data,
          selected_region = reg,
          sample_name = input$sample_name,
          app_version = APP_VERSION
        ),
        file
      )
    }
  )
  
  output$download_curve_fits <- downloadHandler(
    filename = function() {
      paste0(input$sample_name, "_curve_fits_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(results$chromosome_fits)
      fwrite(results$chromosome_fits, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server)