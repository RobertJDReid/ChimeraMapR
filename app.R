library(shiny)
library(data.table)
library(pracma)
library(ggplot2)

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
                   value = 4,
                   min = 1,
                   step = 1),
      helpText("Minimum consecutive same-allele calls to count as a run; increase for noisier data"),

      numericInput("min_peak_height",
                   "Minimum Peak Height:",
                   value = 5,
                   min = 1,
                   step = 1),
#      helpText("Set to ~1/2 of median read depth"),

      numericInput("points_per_window",
                   "Points Per Window:",
                   value = 25,
                   min = 5,
                   step = 1),
      helpText("Number of SNP points per LOESS window"),

#      hr(),

      actionButton("run_analysis",
                   "Run Analysis",
                   class = "btn-primary btn-lg"),

      width = 3
    ),

    mainPanel(
      tabsetPanel(
        id = "main_tabs",

        tabPanel("Overview Plot",
                 helpText("Genome-wide overview of chimeric read coverage and LOESS-smoothed signal across chromosomes."),
                 plotOutput(
                   "chr_plot",
                   height = "1200px"
                 ),
                 downloadButton("download_plot", "Download Plot")
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
                          Colors indicate REF (yellow) vs ALT (purple) alleles."),
                 br(),
                 uiOutput("peak_plots_tabs")
        ),

        tabPanel("LOESS Spans",
                 h4("Chromosome-Specific LOESS Spans"),
                 tableOutput("span_table"),
                 helpText("Shows the span value used for each chromosome in the analysis"),
                 br(),
                 h5("Export LOESS fitted curves"),
                 helpText("Downloads the fitted LOESS curves for all chromosomes from this run, including run parameters for later comparison across runs."),
                 downloadButton("download_loess_fits", "Download LOESS Fits")
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
                   tags$li("Applies LOESS smoothing to identify peaks (recombination hotspots)")
                 ),
                 h5("Input Files:"),
                 tags$ul(
                   tags$li(strong("Read Data:"), "CSV file with SNP position information from BAM file (columns: chrom, pos, read_id, call, is_del, etc.)"),
                   tags$li(strong("SNP Data:"), "CSV file with SNP positions (columns: CHROM, POS, REF, ALT), or a VCF file (plain or gzipped). Multi-allelic VCF sites are split into one row per ALT allele."),
                   tags$li(strong("Chromosome Size:"), "FAI index file with chromosome lengths")
                 ),
                 h5("LOESS Span Methods:"),
                 tags$ul(
                   tags$li(strong("Fixed Span:"), "Uses the same span value for all chromosomes"),
                   tags$li(strong("Dynamic:"), "Calculates chromosome-specific spans based on chromosome length and SNP density. Formula: span = points_per_window × (1/SNP_density) / chromosome_length")
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

  # Set max upload size to 100 MB
  options(shiny.maxRequestSize = 100 * 1024^2)

  # Shared trigger for region plot building — incremented by both the overview
  # button and the per-chromosome buttons
  rv_trigger_region <- reactiveVal(0L)

  # Reactive values to store analysis results
  results <- reactiveValues(
    rt_df                 = NULL,
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

      # ── 1. Load read data ───────────────────────────────────────────────────
      incProgress(0.1, detail = "Loading read data")
      read_data  <- fread(input$read_data_file$datapath)
      chromosomes <- unique(read_data$chrom)

      # ── 2. Load SNP / allele data ───────────────────────────────────────────
      incProgress(0.2, detail = "Loading SNP data")
      snp_path <- input$snp_data_file$datapath
      snp_name <- tolower(input$snp_data_file$name)
      is_vcf   <- grepl("\\.vcf(\\.gz)?$", snp_name)

      if (is_vcf) {
        allele_data <- fread(
          snp_path,
          sep       = "\t",
          skip      = "#",
          col.names = c("CHROM","POS","ID","REF","ALT","QUAL",
                        "FILTER","INFO","FORMAT","DET")
        )
        allele_data <- allele_data[nchar(REF) == 1 & nchar(ALT) == 1,
                                   .(CHROM, POS = as.integer(POS), REF, ALT, QUAL)]
      } else {
        allele_data <- fread(snp_path)
        allele_data <- allele_data[, .(CHROM, POS, REF, ALT)]
      }

      snp_number <- nrow(allele_data)

      # ── 3. Load chromosome size data (FAI) ──────────────────────────────────
      incProgress(0.3, detail = "Loading chromosome sizes")
      chr_size <- fread(
        input$chr_size_file$datapath,
        col.names = c("CHROM", "length", "offset", "col1", "col2")
      )
      chr_size <- chr_size[, .(CHROM, length)]

      # Preserve FASTA chromosome order as factor levels
      fasta_chr_order <- chr_size$CHROM
      genome_size     <- sum(chr_size$length)
      snp_density     <- snp_number / genome_size

      # ── 4. Filter reads and classify alleles ────────────────────────────────
      incProgress(0.4, detail = "Classifying alleles")

      full_read <- read_data[
        mapq      >= input$mapq_cutoff &
        base_qual >= input$baseq_cutoff &
        is_del    == 0
      ]

      full_read <- merge(
        full_read, allele_data,
        by.x  = c("chrom", "pos"),
        by.y  = c("CHROM", "POS"),
        all.x = TRUE
      )

      full_read[, IS_REF := call == REF]
      full_read[, ALLELE := fcase(
        call == REF, "REF",
        call == ALT, "ALT",
        default = "OTHER"
      )]

      full_read <- full_read[ALLELE != "OTHER"]
      setorder(full_read, read_id, pos)

      # ── 5. Detect chimeric reads via RLE ────────────────────────────────────
      incProgress(0.5, detail = "Detecting chimeric reads")
      min_run <- input$min_run

      full_read[, runs := rle_helper(ALLELE), by = read_id]
      full_read <- full_read[runs >= min_run]
      full_read[, new_runs := rle_helper(ALLELE), by = read_id]

      rt_df <- full_read[
        full_read[, .I[.N > min_run & new_runs[1] != .N], by = read_id]$V1
      ]

      results$rt_df             <- rt_df
      results$chimeric_read_ids <- unique(rt_df$read_id)

      # ── 6. Count chimeric reads per SNP position ────────────────────────────
      incProgress(0.6, detail = "Counting SNP coverage")

      pos_count <- rt_df[, .(n = .N), by = .(chrom, pos)]

      snp_coverage <- allele_data[CHROM %in% chromosomes, .(chrom = CHROM, pos = POS)]

      snp_coverage <- merge(
        snp_coverage, pos_count,
        by    = c("chrom", "pos"),
        all.x = TRUE
      )
      snp_coverage[is.na(n), n := 0L]
      snp_coverage[, chrom := factor(chrom, levels = fasta_chr_order)]
      snp_coverage[, pos_kb := pos / 1000]

      results$snp_coverage <- snp_coverage
      
      # ── 7. Calculate LOESS spans per chromosome ─────────────────────────────
      incProgress(0.65, detail = "Calculating chromosome-specific spans")
      
      chr_span <- chr_size[CHROM %in% chromosomes]
      chr_span[, chrom := factor(CHROM, levels = fasta_chr_order)]
      
      chr_span[, lspan := input$points_per_window * (1 / snp_density) / length]
      
      chr_span <- chr_span[, .(chrom, lspan, length)]
      
      results$chr_span <- chr_span
      
      # ── 8. Fit LOESS models and find peaks ──────────────────────────────────
      incProgress(0.7, detail = "Fitting models and finding peaks")

      snp_coverage[, chrom := droplevels(chrom)]
      snp_by_chr <- split(snp_coverage, by = "chrom", keep.by = TRUE)

      span_lookup <- setNames(chr_span$lspan, as.character(chr_span$chrom))
      
      model_results <- lapply(snp_by_chr, function(snps_dt) {
        chr_name <- as.character(snps_dt$chrom[1])
        lspan    <- span_lookup[chr_name]
        
        mdl         <- loess(n ~ pos, data = snps_dt, span = lspan)
        uniform_pos <- seq(min(snps_dt$pos), max(snps_dt$pos), by = 200)
        uniform_fit <- predict(mdl, newdata = data.frame(pos = uniform_pos))
        
        # ── (A) Original findpeaks ──────────────────────────────────────────────
        raw_peaks <- pracma::findpeaks(
          uniform_fit,
          minpeakheight = input$min_peak_height,
          threshold     = 2
        )
        
        findpeaks_positions <- if (!is.null(raw_peaks) && nrow(raw_peaks) > 0) {
          uniform_pos[raw_peaks[, 2]]   # column 2 = index of peak apex
        } else {
          numeric(0)
        }
        
        # ── (B) Region-based fallback ───────────────────────────────────────────
        # Parameters — could be exposed as UI inputs later
        max_region_width_bp <- 50000   # ignore regions wider than 50 Kb (likely baseline noise)
        min_region_width_bp <- 1000        # optional: skip tiny blips
        
        above_thresh <- !is.na(uniform_fit) & uniform_fit >= input$min_peak_height
        
        # Run-length encode the threshold-crossing vector
        rle_above   <- rle(above_thresh)
        region_ends <- cumsum(rle_above$lengths)
        region_starts <- c(1L, head(region_ends, -1) + 1L)
        
        region_peaks <- lapply(seq_along(rle_above$lengths), function(k) {
          if (!rle_above$values[k]) return(NULL)   # below threshold, skip
          
          idx_range <- region_starts[k]:region_ends[k]
          pos_range <- uniform_pos[idx_range]
          region_width <- max(pos_range) - min(pos_range)
          
          # Apply width filter
          if (region_width > max_region_width_bp) return(NULL)
          if (region_width < min_region_width_bp) return(NULL)
          
          # Peak position = argmax within region (or use median(pos_range) for centroid)
          local_max_idx <- idx_range[which.max(uniform_fit[idx_range])]
          
          data.frame(
            peak_height = uniform_fit[local_max_idx],
            peak_index  = local_max_idx,
            peak_start  = region_starts[k],
            peak_end    = region_ends[k]
          )
        })
        
        region_peaks_df <- do.call(rbind, Filter(Negate(is.null), region_peaks))
        
        # ── (C) Deduplicate: keep region peaks not already covered by findpeaks ─
        # A region peak is "already found" if findpeaks placed a peak within its bounds
        if (!is.null(region_peaks_df) && nrow(region_peaks_df) > 0) {
          region_peaks_df$is_new <- vapply(seq_len(nrow(region_peaks_df)), function(i) {
            rstart <- uniform_pos[region_peaks_df$peak_start[i]]
            rend   <- uniform_pos[region_peaks_df$peak_end[i]]
            !any(findpeaks_positions >= rstart & findpeaks_positions <= rend)
          }, logical(1))
          
          novel_region_peaks <- region_peaks_df[region_peaks_df$is_new, , drop = FALSE]
        } else {
          novel_region_peaks <- NULL
        }
        
        # ── (D) Combine findpeaks + novel region peaks ──────────────────────────
        combined_peaks <- raw_peaks   # may be NULL
        
        if (!is.null(novel_region_peaks) && nrow(novel_region_peaks) > 0) {
          novel_mat <- as.matrix(novel_region_peaks[, c("peak_height","peak_index",
                                                        "peak_start","peak_end")])
          combined_peaks <- if (is.null(combined_peaks)) {
            novel_mat
          } else {
            rbind(combined_peaks, novel_mat)
          }
        }
        
        list(
          chrom       = chr_name,
          lspan       = lspan,
          uniform_pos = uniform_pos,
          uniform_fit = uniform_fit,
          peaks       = combined_peaks
        )
      })

      # ── 9. Extract peak positions into a flat data.table ───────────────────
      peaks_list <- lapply(model_results, function(res) {
        if (is.null(res$peaks) || nrow(res$peaks) == 0) return(NULL)

        pk  <- as.data.table(res$peaks)
        setnames(pk, c("peak_height", "peak_index", "peak_start", "peak_end"))

        pk[, chrom      := res$chrom]
        pk[, peak_pos   := res$uniform_pos[peak_index]]
        pk[, peak_start := res$uniform_pos[peak_start]]
        pk[, peak_end   := res$uniform_pos[peak_end]]
        pk[, peak_index := NULL]
        pk
      })

      peaks_genomic <- rbindlist(peaks_list, fill = TRUE)
      if (nrow(peaks_genomic) > 0) {
        peaks_genomic[, chrom := factor(chrom, levels = fasta_chr_order)]
      }

      # ── 10. Map LOESS peaks to nearest actual SNP position ──────────────────
      if (nrow(peaks_genomic) > 0) {
        snp_peaks <- copy(peaks_genomic)
        snp_peaks[, snp_pos := {
          chr_snps <- snp_coverage[chrom == .BY$chrom, pos]
          vapply(peak_pos, function(pp) {
            chr_snps[which.min(abs(chr_snps - pp))]
          }, numeric(1))
        }, by = chrom]
      } else {
        snp_peaks <- peaks_genomic[, snp_pos := numeric(0)]
      }

      results$peaks_genomic <- peaks_genomic
      results$snp_peaks     <- snp_peaks

      # ── 11. Build chromosome_fits for overview plot ─────────────────────────
      incProgress(0.8, detail = "Creating overview plot")

      chromosome_fits <- rbindlist(lapply(model_results, function(res) {
        data.table(
          sample_name       = input$sample_name,
          chrom             = factor(res$chrom, levels = fasta_chr_order),
          uniform_pos       = res$uniform_pos,
          uniform_fit       = res$uniform_fit,
          pos_kb            = res$uniform_pos / 1000,
          span_method       = "dynamic",
          lspan             = res$lspan,
          mapq_cutoff       = input$mapq_cutoff,
          baseq_cutoff      = input$baseq_cutoff,
          min_run           = input$min_run,
          min_peak_height   = input$min_peak_height,
          loess_span        = NA_real_,
          points_per_window = input$points_per_window,
          run_date          = as.character(Sys.Date())
        )
      }))

      results$chromosome_fits <- chromosome_fits

      # ── 12. Overview plot ───────────────────────────────────────────────────
      chr_plot <- ggplot(snp_coverage, aes(x = pos_kb, y = n)) +
        geom_line(
          data  = chromosome_fits,
          aes(x = uniform_pos / 1000, y = uniform_fit),
          color = "firebrick", linewidth = 0.6, alpha = 0.5
        ) +
        geom_point(color = "black", alpha = 0.5, size = 0.5, shape = 21) +
        scale_x_continuous(minor_breaks = seq(0, 1600, 100)) +
        xlab("Position (Kbp)") +
        ylab("Number of Reads") +
        ylim(0, max(30, max(snp_coverage$n))) +
        facet_grid(chrom ~ ., switch = "y") +
        theme_bw() +
        theme(
          panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
          panel.grid.major.x = element_line(linewidth = 0.05, color = "red"),
          strip.background   = element_blank(),
          strip.placement    = "outside"
        )

      if (nrow(peaks_genomic) > 0) {
        chr_plot <- chr_plot +
          geom_vline(
            aes(xintercept = peak_pos / 1000),
            data  = peaks_genomic,
            color = "blue", alpha = 0.5
          )
      }

      results$plot <- chr_plot

      # ── 13. Individual peak plots ───────────────────────────────────────────
      incProgress(0.9, detail = "Creating individual peak plots")

      if (nrow(snp_peaks) > 0) {

        peak_chrs <- levels(droplevels(snp_peaks$chrom))

        peak_plots_list <- lapply(peak_chrs, function(chr_name) {
          chr_peaks  <- snp_peaks[chrom == chr_name]
          chr_rt     <- rt_df[chrom == chr_name]

          peak_plots <- lapply(seq_len(nrow(chr_peaks)), function(i) {
            peak_data <- chr_peaks[i]

            read_ids <- chr_rt[pos == peak_data$snp_pos, unique(read_id)]

            if (length(read_ids) == 0) return(NULL)

            pltdf <- chr_rt[read_id %in% read_ids]

            ggplot(pltdf, aes(x = pos / 1000, y = 1, colour = IS_REF)) +
              geom_vline(
                xintercept = peak_data$snp_pos / 1000,
                color = "grey80", linewidth = 3, alpha = 0.5
              ) +
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
                "Chr ", chr_name, " - Peak ", i,
                " (Position: ", round(peak_data$snp_pos / 1000, 2), " Kb)"
              ))
          })

          peak_plots <- Filter(Negate(is.null), peak_plots)
          list(chromosome = chr_name, plots = peak_plots)
        })

        results$peak_plots_by_chr <- peak_plots_list

      } else {
        results$peak_plots_by_chr <- NULL
      }

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

    snp_coverage    <- copy(results$snp_coverage)
    chromosome_fits <- copy(results$chromosome_fits)
    peaks_genomic   <- results$peaks_genomic

    # Convert chrom to character — Shiny's brush panel detection fails with factors
    snp_coverage[,    chrom := as.character(chrom)]
    chromosome_fits[, chrom := as.character(chrom)]
    if (!is.null(peaks_genomic) && nrow(peaks_genomic) > 0)
      peaks_genomic <- copy(peaks_genomic)[, chrom := as.character(chrom)]

    p <- ggplot(snp_coverage, aes(x = pos_kb, y = n)) +
      geom_line(
        data  = chromosome_fits,
        aes(x = uniform_pos / 1000, y = uniform_fit),
        color = "firebrick", linewidth = 0.6, alpha = 0.5
      ) +
      geom_point(color = "black", alpha = 0.5, size = 0.5, shape = 21) +
      scale_x_continuous(minor_breaks = seq(0, 1600, 100)) +
      xlab("Position (Kbp)") +
      ylab("Number of Reads") +
      ylim(0, max(30, max(snp_coverage$n))) +
      facet_grid(chrom ~ ., switch = "y") +
      theme_bw() +
      theme(
        panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
        panel.grid.major.x = element_line(linewidth = 0.05, color = "red"),
        strip.background   = element_blank(),
        strip.placement    = "outside",
        axis.text          = element_text(size = 12),
        axis.title         = element_text(size = 15),
        strip.text.y       = element_text(size = 9, angle = 0, hjust = 0, face = "bold")
      )

    if (!is.null(peaks_genomic) && nrow(peaks_genomic) > 0) {
      p <- p +
        geom_vline(
          aes(xintercept = peak_pos / 1000),
          data  = peaks_genomic,
          color = "blue", alpha = 0.5
        )
    }

    p
  })

  # Peaks table
  output$peaks_table <- renderTable({
    req(results$peaks_genomic)
    out <- copy(results$peaks_genomic)
    out[, `:=`(
      peak_pos_kb   = round(peak_pos   / 1000, 2),
      peak_start_kb = round(peak_start / 1000, 2),
      peak_end_kb   = round(peak_end   / 1000, 2),
      peak_height   = round(peak_height, 2)
    )]
    out <- out[, .(
      Chromosome           = chrom,
      `Peak Position (Kb)` = peak_pos_kb,
      `Peak Start (Kb)`    = peak_start_kb,
      `Peak End (Kb)`      = peak_end_kb,
      `Peak Height`        = peak_height
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
        plotOutput(paste0("peak_plot_", chr_name, "_", i), height = "600px")
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

  # Render individual peak plots dynamically
  observe({
    req(results$peak_plots_by_chr)

    lapply(results$peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots    <- chr_data$plots

      lapply(seq_along(plots), function(i) {
        output_name <- paste0("peak_plot_", chr_name, "_", i)
        local({
          p <- plots[[i]]
          output[[output_name]] <- renderPlot({ p })
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
          if (!is.null(.peaks) && nrow(.peaks) > 0)
            p <- p + geom_vline(aes(xintercept = peak_pos / 1000),
                                data = .peaks, color = "blue", alpha = 0.6)
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

  # LOESS span table
  output$span_table <- renderTable({
    req(results$chr_span)
    out <- copy(results$chr_span)
    out[, `:=`(
      length_kb = round(length / 1000, 1),
      lspan     = round(lspan, 5)
    )]
    out <- out[, .(
      Chromosome    = chrom,
      `Length (Kb)` = length_kb,
      `LOESS Span`  = lspan
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
    covered_snps     <- results$snp_coverage[n > 0, .N]

    span_method_text <- paste0("Dynamic (Points Per Window: ", input$points_per_window, ")")

    paste0(
      "Sample: ", input$sample_name, "\n\n",
      "Chimeric Reads Detected: ", n_chimeric_reads, "\n",
      "Chromosomes Analyzed: ", n_chromosomes, "\n",
      "Total SNP Positions: ", total_snps, "\n",
      "SNP Positions with Chimeric Reads: ", covered_snps,
      " (", round(100 * covered_snps / total_snps, 1), "%)\n",
      "Peaks Detected: ", n_peaks, "\n\n",
      "Analysis Parameters:\n",
      "  MAPQ Cutoff: ",         input$mapq_cutoff,     "\n",
      "  Base Quality Cutoff: ", input$baseq_cutoff,    "\n",
      "  Min Run Length: ",      input$min_run,         "\n",
      "  Min Peak Height: ",     input$min_peak_height, "\n",
      "  LOESS Span Method: ",   span_method_text
    )
  })

  # ── Shared helper: build and display the selected-region read plot ───────────
  build_region_plot <- function() {
    req(results$selected_region, results$rt_df)

    reg <- results$selected_region

    touching_ids <- results$rt_df[
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
          downloadButton("download_selected_region_plot", "Download Plot"),
          actionButton("close_selected_region_tab", "Close Tab", class = "btn-secondary")
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
      snp_coverage    <- results$snp_coverage
      chromosome_fits <- results$chromosome_fits
      peaks_genomic   <- results$peaks_genomic
      p <- ggplot(snp_coverage, aes(x = pos_kb, y = n)) +
        geom_line(data = chromosome_fits, aes(x = uniform_pos / 1000, y = uniform_fit),
                  color = "firebrick", linewidth = 0.6, alpha = 0.5) +
        geom_point(color = "black", alpha = 0.5, size = 0.5, shape = 21) +
        scale_x_continuous(minor_breaks = seq(0, 1600, 100)) +
        xlab("Position (Kbp)") + ylab("Number of Reads") +
        ylim(0, max(30, max(snp_coverage$n))) +
        facet_grid(chrom ~ ., switch = "y") +
        theme_bw() +
        theme(panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
              panel.grid.major.x = element_line(linewidth = 0.05, color = "red"),
              strip.background = element_blank(), strip.placement = "outside",
              axis.text = element_text(size = 12), axis.title = element_text(size = 15),
              strip.text.y = element_text(size = 9, angle = 0, hjust = 0, face = "bold"))
      if (!is.null(peaks_genomic) && nrow(peaks_genomic) > 0)
        p <- p + geom_vline(aes(xintercept = peak_pos / 1000),
                            data = peaks_genomic, color = "blue", alpha = 0.5)
      ggsave(file, plot = p, width = 12, height = 16, dpi = 300)
    }
  )

  output$download_peaks <- downloadHandler(
    filename = function() paste0(input$sample_name, "_peaks_", Sys.Date(), ".csv"),
    content  = function(file) fwrite(results$peaks_genomic, file)
  )

  output$download_read_ids <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chimeric_read_ids_", Sys.Date(), ".txt"),
    content  = function(file) writeLines(results$chimeric_read_ids, file)
  )

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
  
  output$download_loess_fits <- downloadHandler(
    filename = function() {
      paste0(input$sample_name, "_loess_fits_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(results$chromosome_fits)
      fwrite(results$chromosome_fits, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server)