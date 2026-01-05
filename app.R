library(shiny)
library(tidyverse)
library(pracma)

# Define UI
ui <- fluidPage(
  titlePanel("Chromosome Tracking: Chimeric Read Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Files"),
      fileInput("read_data_file", 
                "Read Data File (CSV/GZ):",
                accept = c(".csv", ".gz", ".csv.gz")),
      
      fileInput("snp_data_file", 
                "SNP Data File (CSV):",
                accept = c(".csv")),
      
      fileInput("chr_size_file", 
                "Chromosome Size File (FAI):",
                accept = c(".fai", ".txt")),
      
      hr(),
      
      h3("Analysis Parameters"),
      textInput("sample_name", 
                "Sample Name:", 
                value = "Sample_01"),
      
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
                     value = 20, 
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
                   tags$li(strong("SNP Data:"), "CSV file with SNP positions (columns: CHROM, POS, REF, ALT)"),
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

# Define server logic
server <- function(input, output, session) {

  # Set max upload size to 100 MB
  options(shiny.maxRequestSize = 100*1024^2)
  
  # Reactive values to store analysis results
  results <- reactiveValues(
    filt_read = NULL,
    snp_coverage = NULL,
    peaks_genomic = NULL,
    snp_peaks = NULL,  # NEW: Store SNP peaks for individual plots
    chromosome_fits = NULL,
    chr_span = NULL,
    plot = NULL,
    chimeric_read_ids = NULL,
    peak_plots_by_chr = NULL  # NEW: Store individual peak plots
  )
  
  # Main analysis function
  observeEvent(input$run_analysis, {
    
    # Validate inputs
    req(input$read_data_file, input$snp_data_file, input$chr_size_file)
    
    withProgress(message = 'Processing data...', value = 0, {
      
      # Load read data
      incProgress(0.1, detail = "Loading read data")
      read_data <- read_csv(input$read_data_file$datapath, show_col_types = FALSE)
      chromosomes <- unique(read_data$chrom)
      
      # Load SNP data
      incProgress(0.2, detail = "Loading SNP data")
      allele_data <- read_csv(input$snp_data_file$datapath, show_col_types = FALSE) %>%
        select(CHROM, POS, REF, ALT)
      
      snp_number <- nrow(allele_data)
      
      # Load chromosome size data
      incProgress(0.3, detail = "Loading chromosome sizes")
      chr_size <- read_table(input$chr_size_file$datapath,
                            col_names = c("CHROM", "length", "offset", "col1", "col2"),
                            show_col_types = FALSE) %>%
        select(CHROM, length)
      
      genome_size <- sum(chr_size$length)
      snp_density <- snp_number / genome_size
      
      # Process read data
      incProgress(0.4, detail = "Classifying alleles")
      full_read <- read_data %>%
        filter(is_del == 0) %>%
        left_join(allele_data, by = join_by(chrom == CHROM, pos == POS)) %>%
        mutate(IS_REF = call == REF,
               ALLELE = case_when(
                 call == REF ~ "REF",
                 call == ALT ~ "ALT",
                 .default = "OTHER"
               ))
      
      # Split data by read_id
      incProgress(0.5, detail = "Detecting chimeric reads")
      read_list <- split(full_read, full_read$read_id)
      
      # Apply run-length encoding
      read_alleles <- map(read_list, `[[`, "ALLELE")
      read_rle <- map(read_alleles, function(x) {
        r <- rle(x)[[1]]
        rep(r, r)
      })
      
      reads_trim <- map2(read_alleles, read_rle, function(l, r) {
        l[which(r != 1)]
      })
      
      reads_f <- map(reads_trim, function(x) {
        length(rle(x)[[1]]) > 1
      })
      
      # Filter for reads with changes
      new_change_list <- map(read_list, function(x) {
        sum(x$IS_REF) > 1 & sum(!x$IS_REF) > 1
      })
      
      ch_ch_changes <- reads_f[which(reads_f == TRUE)]
      
      # Filter read data
      filt_read <- full_read %>%
        filter(read_id %in% names(ch_ch_changes)) %>%
        select(-is_del, -is_refskip) %>%
        select(chrom, pos, read_id, call, IS_REF, everything()) %>%
        arrange(read_id, pos)
      
      results$filt_read <- filt_read
      results$chimeric_read_ids <- names(ch_ch_changes)
      
      # Count SNP positions
      incProgress(0.6, detail = "Counting SNP coverage")
      pos_count <- filt_read %>%
        count(chrom, pos)
      
      snp_coverage <- allele_data %>%
        filter(CHROM %in% chromosomes) %>%
        select(chrom = CHROM, pos = POS) %>%
        left_join(pos_count, by = c("chrom", "pos")) %>%
        replace_na(list(n = 0)) %>%
        mutate(chrom = str_remove(chrom, "S288C_chr")) %>%
        mutate(chrom = factor(chrom, levels = as.character(as.roman(1:16))))
      
      results$snp_coverage <- snp_coverage
      
      # Calculate chromosome-specific spans
      if (input$span_method == "dynamic") {
        incProgress(0.65, detail = "Calculating chromosome-specific spans")
        chr_span <- chr_size %>%
          filter(CHROM %in% chromosomes) %>%
          mutate(chrom = str_remove(CHROM, "S288C_chr")) %>%
          mutate(chrom = factor(chrom, levels = as.character(as.roman(1:16)))) %>%
          mutate(lspan = input$points_per_window * (1/snp_density) / length) %>%
          select(chrom, lspan, length)
        
        results$chr_span <- chr_span
      } else {
        chr_span <- chr_size %>%
          filter(CHROM %in% chromosomes) %>%
          mutate(chrom = str_remove(CHROM, "S288C_chr")) %>%
          mutate(chrom = factor(chrom, levels = as.character(as.roman(1:16)))) %>%
          mutate(lspan = input$loess_span) %>%
          select(chrom, lspan, length)
        
        results$chr_span <- chr_span
      }
      
      # Fit LOESS models by chromosome
      incProgress(0.7, detail = "Fitting models and finding peaks")
      snp_by_chromosome <- snp_coverage %>%
        group_by(chrom) %>%
        nest(.key = "snps") %>%
        left_join(chr_span, by = "chrom")
      
      model_by_chromosome <- snp_by_chromosome %>%
        mutate(
          model = map2(snps, lspan, ~loess(n ~ pos, data = .x, span = .y)),
          uniform_pos = map(snps, ~seq(min(.x$pos), max(.x$pos), by = 200)),
          uniform_fit = map2(model, uniform_pos, predict),
          peaks = map(uniform_fit, pracma::findpeaks, 
                     minpeakheight = input$min_peak_height, 
                     threshold = 2),
          peaks = map(peaks, function(x) {
            if(is.null(x)) return(tibble())
            as_tibble(x, .name_repair = "unique")
          })
        )
      
      # Extract peak positions
      peaks_genomic <- model_by_chromosome %>%
        select(chrom, uniform_pos, peaks) %>%
        unnest(peaks) %>%
        rename("peak_height" = `...1`, "peak_index" = `...2`, 
               "peak_start" = `...3`, "peak_end" = `...4`) %>%
        mutate(
          peak_pos = map2_dbl(uniform_pos, peak_index, ~ .x[.y]),
          peak_start = map2_dbl(uniform_pos, peak_start, ~ .x[.y]),
          peak_end = map2_dbl(uniform_pos, peak_end, ~ .x[.y])
        ) %>%
        select(-uniform_pos)
      
      # Map peaks to actual SNP positions
      snp_peaks <- peaks_genomic %>%
        group_by(chrom) %>%
        nest(.key = "peaks") %>%
        left_join(select(snp_by_chromosome, chrom, snps), by = "chrom") %>%
        unnest(peaks) %>%
        mutate(
          snp_pos = map2_dbl(peak_pos, snps, ~{
            i <- which.min(abs(.y$pos - .x))
            .y$pos[i]
          })
        ) %>%
        select(-snps)
      
      results$peaks_genomic <- peaks_genomic
      results$snp_peaks <- snp_peaks  # NEW: Store for individual plots
      
      # Prepare data for plotting
      incProgress(0.8, detail = "Creating overview plot")
      chromosome_fits <- model_by_chromosome %>%
        select(chrom, uniform_pos, uniform_fit) %>%
        unnest(cols = c(uniform_pos, uniform_fit))
      
      results$chromosome_fits <- chromosome_fits
      
      # Create overview plot
      chr_plot <- snp_coverage %>%
        ggplot(aes(x = pos/1000, y = n)) +
        geom_line(data = chromosome_fits, 
                 aes(x = uniform_pos/1000, y = uniform_fit),
                 color = "firebrick", linewidth = 0.6, alpha = 0.5) +
        geom_point(color = "black", alpha = 0.5, size = 0.5, shape = 21) +
        scale_x_continuous(minor_breaks = seq(0, 1600, 100)) +
        xlab("Position (Kbp)") +
        ylab("Number of Reads") +
        ylim(0, max(30, max(snp_coverage$n))) +
        facet_grid(chrom ~ .) +
        theme_bw() +
        theme(
          panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
          panel.grid.major.x = element_line(linewidth = 0.05, color = "red"),
          strip.background = element_blank()
        )
      
      # Add peak lines if any peaks detected
      if(nrow(peaks_genomic) > 0) {
        chr_plot <- chr_plot + 
          geom_vline(aes(xintercept = peak_pos/1000), 
                    data = peaks_genomic, 
                    color = "blue", 
                    alpha = 0.5)
      }
      
      results$plot <- chr_plot
      
      # NEW: Generate individual peak plots
      incProgress(0.9, detail = "Creating individual peak plots")
      if(nrow(snp_peaks) > 0) {
        # First, prepare a mapping of clean chromosome names to original names
        chr_mapping <- snp_coverage %>%
          select(chrom) %>%
          distinct() %>%
          mutate(original_chrom = chromosomes[match(as.character(chrom), 
                                                     str_remove(chromosomes, "S288C_chr"))])
        
        # Create plots grouped by chromosome
        peak_plots_list <- snp_peaks %>%
          mutate(peak_num = row_number()) %>%
          left_join(chr_mapping, by = "chrom") %>%
          group_by(chrom) %>%
          group_map(~ {
            chr_name <- .y$chrom
            chr_peaks <- .x
            
            # Get the original chromosome name from the mapping
            original_chr <- chr_peaks$original_chrom[1]
            clean_chr <- str_remove(original_chr, "S288C_chr")
            
            # Create a plot for each peak in this chromosome
            peak_plots <- map(1:nrow(chr_peaks), function(i) {
              peak_data <- chr_peaks[i, ]
              
              # Find all reads touching this peak
              read_ids <- filt_read %>%
                filter(chrom == original_chr,
                       pos == peak_data$snp_pos) %>%
                pull(read_id)
              
              # Filter for reads touching this peak
              if(length(read_ids) > 0) {
                pltdf <- filt_read %>%
                  filter(read_id %in% read_ids,
                         chrom == original_chr)
                
                p <- pltdf %>%
                  group_by(read_id) %>%
                  ggplot(aes(x = pos/1000, y = 1, colour = IS_REF)) +
                  geom_vline(xintercept = peak_data$snp_pos/1000, 
                            color = "green", linewidth = 3, alpha = 0.5) +
                  geom_point() +
                  facet_grid(read_id ~ .) +
                  scale_color_viridis_d(option = "turbo", begin = 0.87, end = 0.2) +
                  theme_bw() +
                  theme(
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = "none",
                    plot.background = element_blank(),
                    strip.background = element_blank(),
                    panel.border = element_rect(linewidth = 0.1, linetype = 3),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    strip.text.y = element_blank(),
                    panel.spacing = unit(0, "mm")
                  ) +
                  xlab("Position (Kbp)") +
                  ggtitle(paste0("Chr ", chr_name, " - Peak ", i, 
                                " (Position: ", round(peak_data$snp_pos/1000, 2), " Kb)"))
                
                return(p)
              } else {
                return(NULL)
              }
            })
            
            # Remove NULL plots
            peak_plots <- peak_plots[!sapply(peak_plots, is.null)]
            
            list(chromosome = chr_name, plots = peak_plots)
          }, .keep = TRUE)
        
        results$peak_plots_by_chr <- peak_plots_list
      } else {
        results$peak_plots_by_chr <- NULL
      }
      
      incProgress(1, detail = "Complete")
    })
    
    showNotification("Analysis complete", type = "message", duration = 3)
  })
  
  # Output: Main plot
  output$chr_plot <- renderPlot({
    req(results$plot)
    results$plot +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
      )
  })
  
  # Output: Peaks table
  output$peaks_table <- renderTable({
    req(results$peaks_genomic)
    results$peaks_genomic %>%
      mutate(
        peak_pos_kb = round(peak_pos/1000, 2),
        peak_start_kb = round(peak_start/1000, 2),
        peak_end_kb = round(peak_end/1000, 2),
        peak_height = round(peak_height, 2)
      ) %>%
      select(Chromosome = chrom, 
             `Peak Position (Kb)` = peak_pos_kb,
             `Peak Start (Kb)` = peak_start_kb,
             `Peak End (Kb)` = peak_end_kb,
             `Peak Height` = peak_height) %>%
      arrange(Chromosome, `Peak Position (Kb)`)
  })
  
  # NEW: Dynamic UI for peak plots organized by chromosome
  output$peak_plots_tabs <- renderUI({
    req(results$peak_plots_by_chr)
    
    if(is.null(results$peak_plots_by_chr) || length(results$peak_plots_by_chr) == 0) {
      return(h4("No peaks detected in the analysis."))
    }
    
    # Create tabs for each chromosome
    chr_tabs <- map(results$peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots <- chr_data$plots
      
      if(length(plots) == 0) {
        return(NULL)
      }
      
      # Create plot outputs for this chromosome
      plot_outputs <- map(1:length(plots), function(i) {
        output_name <- paste0("peak_plot_", chr_name, "_", i)
        plotOutput(output_name, height = "600px")
      })
      
      tabPanel(
        paste("Chr", chr_name),
        h5(paste0("Chromosome ", chr_name, " - ", length(plots), " peak(s) detected")),
        hr(),
        do.call(tagList, plot_outputs)
      )
    })
    
    # Remove NULL tabs
    chr_tabs <- chr_tabs[!sapply(chr_tabs, is.null)]
    
    if(length(chr_tabs) == 0) {
      return(h4("No peaks detected in the analysis."))
    }
    
    do.call(tabsetPanel, chr_tabs)
  })
  
  # NEW: Render individual peak plots dynamically
  observe({
    req(results$peak_plots_by_chr)
    
    walk(results$peak_plots_by_chr, function(chr_data) {
      chr_name <- chr_data$chromosome
      plots <- chr_data$plots
      
      walk(1:length(plots), function(i) {
        output_name <- paste0("peak_plot_", chr_name, "_", i)
        
        output[[output_name]] <- renderPlot({
          plots[[i]]
        })
      })
    })
  })
  
  # Output: Chromosome span table
  output$span_table <- renderTable({
    req(results$chr_span)
    results$chr_span %>%
      mutate(
        length_kb = round(length/1000, 1),
        lspan = round(lspan, 5)
      ) %>%
      select(Chromosome = chrom,
             `Length (Kb)` = length_kb,
             `LOESS Span` = lspan) %>%
      arrange(Chromosome)
  })
  
  # Output: Summary statistics
  output$summary_stats <- renderText({
    req(results$filt_read, results$snp_coverage, results$peaks_genomic)
    
    n_chimeric_reads <- length(unique(results$filt_read$read_id))
    n_chromosomes <- length(unique(results$snp_coverage$chrom))
    n_peaks <- nrow(results$peaks_genomic)
    total_snps <- nrow(results$snp_coverage)
    covered_snps <- sum(results$snp_coverage$n > 0)
    
    span_method_text <- if(input$span_method == "dynamic") {
      paste0("Dynamic (Points Per Window: ", input$points_per_window, ")")
    } else {
      paste0("Fixed (Span: ", input$loess_span, ")")
    }
    
    paste0(
      "Sample: ", input$sample_name, "\n\n",
      "Chimeric Reads Detected: ", n_chimeric_reads, "\n",
      "Chromosomes Analyzed: ", n_chromosomes, "\n",
      "Total SNP Positions: ", total_snps, "\n",
      "SNP Positions with Chimeric Reads: ", covered_snps, " (", 
      round(100 * covered_snps/total_snps, 1), "%)\n",
      "Peaks Detected: ", n_peaks, "\n\n",
      "Analysis Parameters:\n",
      "  Min Peak Height: ", input$min_peak_height, "\n",
      "  LOESS Span Method: ", span_method_text
    )
  })
  
  # Download handlers
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0(input$sample_name, "_chromosome_tracking_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = results$plot, width = 12, height = 16, dpi = 300)
    }
  )
  
  output$download_peaks <- downloadHandler(
    filename = function() {
      paste0(input$sample_name, "_peaks_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(results$peaks_genomic, file)
    }
  )

  output$download_read_ids <- downloadHandler(
    filename = function() {
      paste0(input$sample_name, "_chimeric_read_ids_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write_lines(results$chimeric_read_ids, file)
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)
