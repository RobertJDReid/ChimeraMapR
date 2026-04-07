library(shiny)
library(data.table)
library(pracma)
library(ggplot2)

# ─────────────────────────────────────────────
#                   UI
# ─────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("ChimeraMapR: Chimeric SNP Detection in Long-Read Sequencing Data"),
  
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
      
      hr(),
      
      h3("Analysis Parameters"),
      textInput("sample_name", 
                "Sample Name:", 
                value = "Sample_01"),
      
      numericInput("mapq_cutoff",
                   "MAPQ Cutoff:",
                   value = 20,
                   min = 0,
                   step = 1),
      helpText("Minimum mapping quality score; reads below this are excluded"),
      
      numericInput("baseq_cutoff",
                   "Base Quality Cutoff:",
                   value = 10,
                   min = 0,
                   step = 1),
      helpText("Minimum base quality score at SNP position"),
      
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
      helpText("Set to ~1/2 of median read depth"),
      
      radioButtons("span_method",
                   "LOESS Span Method:",
                   choices = c("Fixed Span" = "fixed",
                              "Dynamic (Per Chromosome)" = "dynamic"),
                   selected = "dynamic"),
      
      conditionalPanel(
        condition = "input.span_method == 'fixed'",
        numericInput("loess_span", 
                     "LOESS Span:", 
                     value = 0.015, 
                     min = 0.001, 
                     max = 1, 
                     step = 0.005),
        helpText("Lower values for low coverage datasets (typically 0.015)")
      ),
      
      conditionalPanel(
        condition = "input.span_method == 'dynamic'",
        numericInput("points_per_window", 
                     "Points Per Window:", 
                     value = 25, 
                     min = 5, 
                     step = 1),
        helpText("Number of SNP points per LOESS window. Span calculated as: ppw × (1/SNP_density) / chr_length")
      ),
      
      hr(),
      
      actionButton("run_analysis", 
                   "Run Analysis", 
                   class = "btn-primary btn-lg"),
      
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Overview Plot",
                 plotOutput("chr_plot", height = "1200px"),
                 downloadButton("download_plot", "Download Plot")
        ),
        tabPanel("Peak Summary",
                 h4("Detected Peaks"),
                 tableOutput("peaks_table"),
                 downloadButton("download_peaks", "Download Peak Data")
        ),
        # NEW: Individual Peak Plots tab with nested tabs by chromosome
        tabPanel("Individual Peak Plots",
                 h4("Detailed Peak Visualizations by Chromosome"),
                 helpText("Each peak shows all chimeric reads that intersect that position. 
                          Colors indicate REF (yellow) vs ALT (purple) alleles."),
                 br(),
                 uiOutput("peak_plots_tabs")
        ),
        tabPanel("Read Statistics",
                 h4("Analysis Summary"),
                 verbatimTextOutput("summary_stats"),
                 br(),
                 downloadButton("download_read_ids", "Download Chimeric Read IDs")
        ),
        tabPanel("LOESS Spans",
                 h4("Chromosome-Specific LOESS Spans"),
                 tableOutput("span_table"),
                 helpText("Shows the span value used for each chromosome in the analysis")
        ),
        tabPanel("About",
                 h4("About This Analysis"),
                 p("This application identifies chimeric reads in sequencing data by tracking 
                   allele changes across chromosomes. Chimeric reads contain sequence from 
                   multiple parental chromosomes and can indicate recombination events."),
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

  # Reactive values to store analysis results
  results <- reactiveValues(
    rt_df            = NULL,
    snp_coverage     = NULL,
    peaks_genomic    = NULL,
    snp_peaks        = NULL,
    chromosome_fits  = NULL,
    chr_span         = NULL,
    plot             = NULL,
    chimeric_read_ids = NULL,
    peak_plots_by_chr = NULL
  )

  # ── RLE helper (unchanged – operates on plain vectors) ──────────────────────
  rle_helper <- function(x) {
    r  <- rle(x)[[1]]   # run lengths
    rn <- rep(r, r)     # expand to per-position run lengths
    return(rn)
  }

  # ── Main analysis ────────────────────────────────────────────────────────────
  observeEvent(input$run_analysis, {

    req(input$read_data_file, input$snp_data_file, input$chr_size_file)

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
        # fread skips comment lines that start with "#" via skip="#"
        allele_data <- fread(
          snp_path,
          sep       = "\t",
          skip      = "#",
          col.names = c("CHROM","POS","ID","REF","ALT","QUAL",
                        "FILTER","INFO","FORMAT","DET")
        )
        # Keep only single-nucleotide variants and required columns
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

      # Apply quality filters and remove deletions
      full_read <- read_data[
        mapq     >= input$mapq_cutoff  &
        base_qual >= input$baseq_cutoff &
        is_del   == 0
      ]

      # Join allele info (left join: keep all read rows, add REF/ALT/QUAL)
      full_read <- merge(
        full_read, allele_data,
        by.x  = c("chrom", "pos"),
        by.y  = c("CHROM", "POS"),
        all.x = TRUE
      )

      # Classify each base call
      full_read[, IS_REF := call == REF]
      full_read[, ALLELE := fcase(
        call == REF, "REF",
        call == ALT, "ALT",
        default = "OTHER"
      )]

      # Drop rows that are neither REF nor ALT
      full_read <- full_read[ALLELE != "OTHER"]

      # Sort by read then position
      setorder(full_read, read_id, pos)

      # ── 5. Detect chimeric reads via RLE ────────────────────────────────────
      incProgress(0.5, detail = "Detecting chimeric reads")
      min_run <- input$min_run

      # First RLE pass: compute run lengths per read
      full_read[, runs := rle_helper(ALLELE), by = read_id]

      # Remove short (noisy) runs
      full_read <- full_read[runs > min_run]

      # Second RLE pass on denoised data
      full_read[, new_runs := rle_helper(ALLELE), by = read_id]

      # Keep reads that:
      #   (a) still have more than min_run rows after denoising, AND
      #   (b) do not consist of a single uninterrupted run
      #       (new_runs[1] != .N means the first run doesn't span everything)
      rt_df <- full_read[
        full_read[, .I[.N > min_run & new_runs[1] != .N], by = read_id]$V1
      ]

      results$rt_df             <- rt_df
      results$chimeric_read_ids <- unique(rt_df$read_id)

      # ── 6. Count chimeric reads per SNP position ────────────────────────────
      incProgress(0.6, detail = "Counting SNP coverage")

      pos_count <- rt_df[, .(n = .N), by = .(chrom, pos)]

      # Start from all SNP positions on chromosomes present in read data
      snp_coverage <- allele_data[CHROM %in% chromosomes, .(chrom = CHROM, pos = POS)]

      snp_coverage <- merge(
        snp_coverage, pos_count,
        by    = c("chrom", "pos"),
        all.x = TRUE
      )
      snp_coverage[is.na(n), n := 0L]
      snp_coverage[, chrom := factor(chrom, levels = fasta_chr_order)]

      results$snp_coverage <- snp_coverage

      # ── 7. Calculate LOESS spans per chromosome ─────────────────────────────
      incProgress(0.65, detail = "Calculating chromosome-specific spans")

      chr_span <- chr_size[CHROM %in% chromosomes]
      chr_span[, chrom := factor(CHROM, levels = fasta_chr_order)]

      if (input$span_method == "dynamic") {
        chr_span[, lspan := input$points_per_window * (1 / snp_density) / length]
      } else {
        chr_span[, lspan := input$loess_span]
      }
      chr_span <- chr_span[, .(chrom, lspan, length)]

      results$chr_span <- chr_span

      # ── 8. Fit LOESS models and find peaks ──────────────────────────────────
      # data.table has no list-column / nest equivalent so we use split + lapply
      incProgress(0.7, detail = "Fitting models and finding peaks")

      # Split coverage into a list of data.tables, one per chromosome.
      # Drop unused factor levels first — split() on a factor produces an empty
      # data.table for every level, including chromosomes not in the data, which
      # causes loess() to crash with "invalid 'x'" on the empty slices.
      snp_coverage[, chrom := droplevels(chrom)] # added to squash bug
      snp_by_chr <- split(snp_coverage, by = "chrom", keep.by = TRUE)

      # Build a named lookup for spans
      span_lookup <- setNames(chr_span$lspan, as.character(chr_span$chrom))

      model_results <- lapply(snp_by_chr, function(snps_dt) {
        chr_name <- as.character(snps_dt$chrom[1])
        lspan    <- span_lookup[chr_name]

        # Fit LOESS
        mdl         <- loess(n ~ pos, data = snps_dt, span = lspan)
        uniform_pos <- seq(min(snps_dt$pos), max(snps_dt$pos), by = 200)
        uniform_fit <- predict(mdl, newdata = data.frame(pos = uniform_pos))

        # Find peaks
        raw_peaks <- pracma::findpeaks(
          uniform_fit,
          minpeakheight = input$min_peak_height,
          threshold     = 2
        )

        list(
          chrom       = chr_name,
          uniform_pos = uniform_pos,
          uniform_fit = uniform_fit,
          peaks       = raw_peaks    # matrix or NULL
        )
      })

      # ── 9. Extract peak positions into a flat data.table ───────────────────
      peaks_list <- lapply(model_results, function(res) {
        if (is.null(res$peaks) || nrow(res$peaks) == 0) return(NULL)

        pk  <- as.data.table(res$peaks)
        setnames(pk, c("peak_height", "peak_index", "peak_start", "peak_end"))

        pk[, chrom       := res$chrom]
        pk[, peak_pos    := res$uniform_pos[peak_index]]
        pk[, peak_start  := res$uniform_pos[peak_start]]
        pk[, peak_end    := res$uniform_pos[peak_end]]
        pk[, peak_index  := NULL]
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
          chrom       = factor(res$chrom, levels = fasta_chr_order),
          uniform_pos = res$uniform_pos,
          uniform_fit = res$uniform_fit
        )
      }))

      results$chromosome_fits <- chromosome_fits

      # ── 12. Overview plot ───────────────────────────────────────────────────
      chr_plot <- ggplot(snp_coverage, aes(x = pos / 1000, y = n)) +
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
        facet_grid(chrom ~ ., switch = "y", scales = "free_x") +
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

        # Unique chromosomes that have peaks, in factor order
        peak_chrs <- levels(droplevels(snp_peaks$chrom))

        peak_plots_list <- lapply(peak_chrs, function(chr_name) {
          chr_peaks  <- snp_peaks[chrom == chr_name]
          chr_rt     <- rt_df[chrom == chr_name]

          peak_plots <- lapply(seq_len(nrow(chr_peaks)), function(i) {
            peak_data <- chr_peaks[i]

            # Reads that cover this exact peak SNP position
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
                axis.text.y    = element_blank(),
                axis.ticks.y   = element_blank(),
                axis.title.y   = element_blank(),
                legend.position = "none",
                plot.background = element_blank(),
                strip.background = element_blank(),
                panel.border   = element_rect(linewidth = 0.1, linetype = 3),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text.y   = element_blank(),
                panel.spacing  = unit(0, "mm")
              ) +
              xlab("Position (Kbp)") +
              ggtitle(paste0(
                "Chr ", chr_name, " - Peak ", i,
                " (Position: ", round(peak_data$snp_pos / 1000, 2), " Kb)"
              ))
          })

          # Drop NULLs
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

  # ── Outputs ──────────────────────────────────────────────────────────────────

  # Overview plot
  output$chr_plot <- renderPlot({
    req(results$plot)
    results$plot +
      theme(
        axis.text    = element_text(size = 12),
        axis.title   = element_text(size = 15),
        strip.text.y = element_text(size = 9, angle = 0, hjust = 0, face = "bold")
      )
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

    span_method_text <- if (input$span_method == "dynamic") {
      paste0("Dynamic (Points Per Window: ", input$points_per_window, ")")
    } else {
      paste0("Fixed (Span: ", input$loess_span, ")")
    }

    paste0(
      "Sample: ", input$sample_name, "\n\n",
      "Chimeric Reads Detected: ", n_chimeric_reads, "\n",
      "Chromosomes Analyzed: ", n_chromosomes, "\n",
      "Total SNP Positions: ", total_snps, "\n",
      "SNP Positions with Chimeric Reads: ", covered_snps,
      " (", round(100 * covered_snps / total_snps, 1), "%)\n",
      "Peaks Detected: ", n_peaks, "\n\n",
      "Analysis Parameters:\n",
      "  MAPQ Cutoff: ",        input$mapq_cutoff,  "\n",
      "  Base Quality Cutoff: ", input$baseq_cutoff, "\n",
      "  Min Run Length: ",     input$min_run,       "\n",
      "  Min Peak Height: ",    input$min_peak_height, "\n",
      "  LOESS Span Method: ",  span_method_text
    )
  })

  # ── Download handlers ────────────────────────────────────────────────────────

  output$download_plot <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chromosome_tracking_", Sys.Date(), ".png"),
    content  = function(file) ggsave(file, plot = results$plot, width = 12, height = 16, dpi = 300)
  )

  output$download_peaks <- downloadHandler(
    filename = function() paste0(input$sample_name, "_peaks_", Sys.Date(), ".csv"),
    content  = function(file) fwrite(results$peaks_genomic, file)
  )

  output$download_read_ids <- downloadHandler(
    filename = function() paste0(input$sample_name, "_chimeric_read_ids_", Sys.Date(), ".txt"),
    content  = function(file) writeLines(results$chimeric_read_ids, file)
  )
}

# Run the application
shinyApp(ui = ui, server = server)
