# =============================================================================
#  chimera_functions.R
#  Pure analysis functions for ChimeraMapR.
#
#  No Shiny dependency — safe to source from app.R or any CLI/batch script.
#
#  Main entry point:
#    results <- run_chimera_analysis(read_path, snp_path, fai_path, ...)
#
#  Plot builder (used by both the Shiny download handler and the CLI):
#    p <- build_overview_plot(results)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(pracma)
  library(ggplot2)
  library(Matrix)
})

APP_VERSION <- "0.5.0"


# -----------------------------------------------------------------------------
#  Whittaker smoother
#  y      : numeric vector (the signal to smooth)
#  lambda : smoothness penalty — smaller = tighter fit, larger = smoother
#  d      : penalty derivative order (2 = standard, penalises curvature)
#
#  Based on Eilers P, Bloemberg T, Wehrens R (2026).
#  _ptw: Parametric Time Warping_. doi:10.32614/CRAN.package.ptw
# -----------------------------------------------------------------------------
whittaker <- function(y, lambda = 1, d = 2) {
  n <- length(y)
  if (n < d + 2L) return(y)
  E <- Matrix::Diagonal(n)
  D <- Matrix::diff(E, differences = d)
  as.numeric(Matrix::solve(E + lambda * Matrix::t(D) %*% D, y))
}


# -----------------------------------------------------------------------------
#  RLE helper — returns per-position run lengths for a vector x
# -----------------------------------------------------------------------------
rle_helper <- function(x) {
  r  <- rle(x)[[1]]
  rn <- rep(r, r)
  rn
}


# -----------------------------------------------------------------------------
#  File loaders
# -----------------------------------------------------------------------------

#' Load read data CSV (plain or gzipped)
load_read_data <- function(path) {
  fread(path)
}

#' Load SNP data — auto-detects VCF (.vcf, .vcf.gz) vs plain CSV.
#' Returns a data.table with columns: CHROM, POS, REF, ALT  (+ QUAL for VCF)
load_snp_data <- function(path) {
  is_vcf <- grepl("\\.vcf(\\.gz)?$", tolower(path))

  if (is_vcf) {
    dt <- fread(
      path,
      sep       = "\t",
      skip      = "#",
      col.names = c("CHROM", "POS", "ID", "REF", "ALT",
                    "QUAL", "FILTER", "INFO", "FORMAT", "DET")
    )
    # Keep only biallelic SNPs (single-character REF and ALT)
    dt <- dt[nchar(REF) == 1 & nchar(ALT) == 1,
             .(CHROM, POS = as.integer(POS), REF, ALT, QUAL)]
  } else {
    dt <- fread(path)
    dt <- dt[, .(CHROM, POS, REF, ALT)]
  }
  dt
}

#' Load chromosome size from a FASTA index (.fai).
#' Returns a data.table with columns: CHROM, length
load_chr_size <- function(path) {
  dt <- fread(path,
              col.names = c("CHROM", "length", "offset", "col1", "col2"))
  dt[, .(CHROM, length)]
}


# -----------------------------------------------------------------------------
#  Core analysis pipeline
#
#  Returns a named list that mirrors the `results` reactiveValues in app.R:
#    $rt_df              — chimeric reads (full allele table)
#    $chimeric_read_ids  — unique read IDs
#    $transition_pos     — boundary positions
#    $snp_coverage       — per-SNP transition counts
#    $peaks_genomic      — raw peak table from findpeaks / fallback
#    $snp_peaks          — peaks mapped to best qualifying SNP
#    $chromosome_fits    — Whittaker fitted curves (exportable)
#    $chr_span           — per-chromosome lambda / SNP count table
#    $params             — the parameter list used for this run
# -----------------------------------------------------------------------------
run_chimera_analysis <- function(
    read_data_path,
    snp_data_path,
    chr_size_path,
    sample_name     = "Sample_01",
    mapq_cutoff     = 20L,
    baseq_cutoff    = 10L,
    min_run         = 2L,
    min_peak_height = 10L,
    lambda          = 1,
    warn_fn         = function(msg) message("WARNING: ", msg)
) {

  # ── 1. Load data ─────────────────────────────────────────────────────────────
  message("  Loading read data ...")
  read_data   <- load_read_data(read_data_path)
  chromosomes <- unique(read_data$chrom)

  message("  Loading SNP data ...")
  allele_data <- load_snp_data(snp_data_path)

  message("  Loading chromosome sizes ...")
  chr_size    <- load_chr_size(chr_size_path)

  # Preserve FASTA chromosome order; restrict to chromosomes in read data
  fasta_chr_order  <- chr_size$CHROM
  allele_data_used <- allele_data[CHROM %in% chromosomes]
  chr_size_used    <- chr_size[CHROM %in% chromosomes]

  snp_number   <- nrow(allele_data_used)
  genome_size  <- sum(chr_size_used$length)
  snp_density  <- snp_number / genome_size

  # ── 2. Filter reads and classify alleles ─────────────────────────────────────
  message("  Filtering reads and classifying alleles ...")
  full_read <- read_data[
    mapq      >= mapq_cutoff  &
    base_qual >= baseq_cutoff &
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
  full_read <- full_read[, .(chrom, pos, read_id, IS_REF, ALLELE)]
  setorder(full_read, read_id, pos)

  # ── 3. Detect chimeric reads via RLE ─────────────────────────────────────────
  message("  Detecting chimeric reads ...")
  full_read[, runs := rle_helper(ALLELE), by = read_id]
  full_read <- full_read[runs >= min_run]
  full_read[, new_runs := rle_helper(ALLELE), by = read_id]

  rt_df <- full_read[
    full_read[, .I[.N > min_run & new_runs[1] != .N], by = read_id]$V1
  ]
  rt_df[, c("runs", "new_runs") := NULL]

  chimeric_read_ids <- unique(rt_df$read_id)

  # ── 4. Extract transition boundaries ─────────────────────────────────────────
  message("  Extracting transition boundaries ...")
  transition_pos <- rt_df[, {
    is_last_of_run  <- c(ALLELE[-1] != ALLELE[-.N], FALSE)
    is_first_of_run <- c(FALSE, ALLELE[-1] != ALLELE[-.N])
    .SD[is_last_of_run | is_first_of_run, .(chrom, pos)]
  }, by = read_id]

  pos_count <- transition_pos[, .(n = .N), by = .(chrom, pos)]

  snp_coverage <- allele_data[CHROM %in% chromosomes, .(chrom = CHROM, pos = POS)]
  snp_coverage <- merge(snp_coverage, pos_count,
                        by = c("chrom", "pos"), all.x = TRUE)
  snp_coverage[is.na(n), n := 0L]
  snp_coverage[, chrom  := factor(chrom, levels = fasta_chr_order)]
  snp_coverage[, pos_kb := pos / 1000]

  # ── 5. Build per-chromosome lambda table ─────────────────────────────────────
  snp_counts_by_chr <- allele_data[CHROM %in% chromosomes, .N, by = CHROM]
  chr_span <- chr_size[CHROM %in% chromosomes]
  chr_span[, chrom := factor(CHROM, levels = fasta_chr_order)]
  chr_span <- merge(chr_span, snp_counts_by_chr, by = "CHROM", all.x = TRUE)
  chr_span[is.na(N), N := 0L]
  chr_span[, lambda := lambda]
  chr_span <- chr_span[, .(chrom, lambda, length, snp_count = N)]

  # ── 6. Whittaker smoothing and peak finding (per chromosome) ─────────────────
  message("  Fitting Whittaker smoother and finding peaks ...")
  snp_coverage[, chrom := droplevels(chrom)]
  snp_by_chr <- split(snp_coverage, by = "chrom", keep.by = TRUE)

  model_results <- lapply(snp_by_chr, function(snps_dt) {
    snps_dt     <- snps_dt[order(pos)]
    chr_name    <- as.character(snps_dt$chrom[1])
    uniform_pos <- seq(min(snps_dt$pos), max(snps_dt$pos), by = 200)

    uniform_fit <- tryCatch({
      med_step    <- if (nrow(snps_dt) > 1L) median(diff(snps_dt$pos)) else 200
      scaled_lam  <- lambda * (200 / med_step)^2
      smoothed    <- whittaker(snps_dt$n, lambda = scaled_lam, d = 2)
      approx(snps_dt$pos, smoothed, xout = uniform_pos, rule = 2)$y
    }, error = function(e) {
      warn_fn(paste0(
        "Whittaker smoother failed for ", chr_name,
        " (too few SNPs) — raw data shown, no curve."
      ))
      rep(NA_real_, length(uniform_pos))
    })

    # (A) findpeaks on the Whittaker fit
    raw_peaks <- if (all(is.na(uniform_fit))) NULL else pracma::findpeaks(
      uniform_fit,
      minpeakheight = 1,
      threshold     = 0
    )

    findpeaks_positions <- if (!is.null(raw_peaks) && nrow(raw_peaks) > 0) {
      uniform_pos[raw_peaks[, 2]]
    } else {
      numeric(0)
    }

    # (B) Raw-signal fallback: cluster high-count SNP positions
    cluster_gap_bp <- 10000L
    raw_hi         <- snps_dt[n >= min_peak_height]

    if (nrow(raw_hi) > 0) {
      raw_hi  <- raw_hi[order(raw_hi$pos), ]
      gaps    <- c(cluster_gap_bp + 1L, diff(raw_hi$pos))
      raw_hi$cl <- cumsum(gaps > cluster_gap_bp)

      region_peaks_df <- do.call(rbind, lapply(split(raw_hi, raw_hi$cl), function(cl_rows) {
        best <- cl_rows[which.max(cl_rows$n), ]
        ui   <- which.min(abs(uniform_pos - best$pos))
        data.frame(
          peak_height = best$n,
          peak_index  = ui,
          peak_start  = max(1L, ui - 5L),
          peak_end    = min(length(uniform_pos), ui + 5L)
        )
      }))
    } else {
      region_peaks_df <- NULL
    }

    # (C) Deduplicate — keep fallback peaks not already found by findpeaks
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

    # (D) Combine findpeaks + novel fallback peaks
    combined_peaks <- raw_peaks
    if (!is.null(novel_region_peaks) && nrow(novel_region_peaks) > 0) {
      novel_mat <- as.matrix(novel_region_peaks[,
        c("peak_height", "peak_index", "peak_start", "peak_end")])
      combined_peaks <- if (is.null(combined_peaks)) novel_mat else
                         rbind(combined_peaks, novel_mat)
    }

    list(
      chrom       = chr_name,
      lambda      = lambda,
      uniform_pos = uniform_pos,
      uniform_fit = uniform_fit,
      peaks       = combined_peaks
    )
  })

  # ── 7. Extract peak positions into a flat data.table ─────────────────────────
  peaks_list <- lapply(model_results, function(res) {
    if (is.null(res$peaks) || nrow(res$peaks) == 0) return(NULL)
    pk <- as.data.table(res$peaks)
    setnames(pk, c("peak_height", "peak_index", "peak_start", "peak_end"))
    pk[, chrom      := res$chrom]
    pk[, peak_pos   := res$uniform_pos[peak_index]]
    pk[, peak_start := res$uniform_pos[peak_start]]
    pk[, peak_end   := res$uniform_pos[peak_end]]
    pk[, peak_index := NULL]
    pk
  })
  peaks_genomic <- rbindlist(peaks_list, fill = TRUE)
  if (nrow(peaks_genomic) > 0)
    peaks_genomic[, chrom := factor(chrom, levels = fasta_chr_order)]

  # ── 8. Map each peak interval to its best qualifying SNP ─────────────────────
  message("  Mapping peaks to SNP positions ...")
  if (nrow(peaks_genomic) > 0) {

    snp_peaks_list <- lapply(seq_len(nrow(peaks_genomic)), function(i) {
      row      <- peaks_genomic[i]
      chr_name <- as.character(row$chrom)
      chr_cov  <- snp_coverage[as.character(chrom) == chr_name]

      in_interval <- chr_cov[
        pos >= row$peak_start &
        pos <= row$peak_end   &
        n   >= min_peak_height
      ]

      if (nrow(in_interval) > 0) {
        max_n      <- max(in_interval$n)
        candidates <- in_interval[n == max_n]

        if (nrow(candidates) > 1L) {
          candidates[, dist_to_peak := abs(pos - row$peak_pos)]
          min_dist   <- min(candidates$dist_to_peak)
          candidates <- candidates[dist_to_peak == min_dist]
          candidates[, dist_to_peak := NULL]
        }

        best <- if (nrow(candidates) == 1L) candidates else candidates[sample(.N, 1L)]

        data.table(
          chrom       = row$chrom,
          peak_pos    = row$peak_pos,
          peak_height = row$peak_height,
          peak_start  = row$peak_start,
          peak_end    = row$peak_end,
          snp_pos     = best$pos,
          snp_n       = best$n
        )
      } else {
        data.table(
          chrom       = row$chrom,
          peak_pos    = row$peak_pos,
          peak_height = row$peak_height,
          peak_start  = row$peak_start,
          peak_end    = row$peak_end,
          snp_pos     = NA_real_,
          snp_n       = NA_integer_
        )
      }
    })

    snp_peaks <- rbindlist(snp_peaks_list, fill = TRUE)
    snp_peaks[, chrom := factor(as.character(chrom),
                                levels = levels(peaks_genomic$chrom))]

    snp_peaks[, chimeric_reads_at_snp := {
      chr_name <- as.character(.BY$chrom)
      vapply(snp_pos, function(sp) {
        if (is.na(sp)) return(NA_integer_)
        uniqueN(rt_df[as.character(chrom) == chr_name & pos == sp, read_id])
      }, integer(1))
    }, by = chrom]

  } else {
    snp_peaks <- data.table(
      chrom                 = factor(character(0), levels = fasta_chr_order),
      peak_pos              = numeric(0),
      peak_height           = numeric(0),
      peak_start            = numeric(0),
      peak_end              = numeric(0),
      snp_pos               = numeric(0),
      snp_n                 = integer(0),
      chimeric_reads_at_snp = integer(0)
    )
  }

  # ── 9. Build chromosome_fits table (exportable curve data) ───────────────────
  chromosome_fits <- rbindlist(lapply(model_results, function(res) {
    data.table(
      sample_name     = sample_name,
      chrom           = factor(res$chrom, levels = fasta_chr_order),
      uniform_pos     = res$uniform_pos,
      uniform_fit     = res$uniform_fit,
      pos_kb          = res$uniform_pos / 1000,
      smoother        = "whittaker",
      lambda          = res$lambda,
      mapq_cutoff     = mapq_cutoff,
      baseq_cutoff    = baseq_cutoff,
      min_run         = min_run,
      min_peak_height = min_peak_height,
      run_date        = as.character(Sys.Date())
    )
  }))

  # ── Return all results ────────────────────────────────────────────────────────
  list(
    rt_df             = rt_df,
    chimeric_read_ids = chimeric_read_ids,
    transition_pos    = transition_pos,
    snp_coverage      = snp_coverage,
    peaks_genomic     = peaks_genomic,
    snp_peaks         = snp_peaks,
    chromosome_fits   = chromosome_fits,
    chr_span          = chr_span,
    params = list(
      sample_name     = sample_name,
      mapq_cutoff     = mapq_cutoff,
      baseq_cutoff    = baseq_cutoff,
      min_run         = min_run,
      min_peak_height = min_peak_height,
      lambda          = lambda,
      app_version     = APP_VERSION
    )
  )
}


# -----------------------------------------------------------------------------
#  Overview plot builder
#
#  Identical logic to the Shiny renderPlot / download_plot handler, but
#  expressed as a plain function that returns a ggplot object.
#
#  results : the list returned by run_chimera_analysis()
# -----------------------------------------------------------------------------
build_overview_plot <- function(results) {
  snp_coverage    <- copy(results$snp_coverage)
  chromosome_fits <- copy(results$chromosome_fits)
  peaks_genomic   <- results$peaks_genomic
  snp_peaks       <- results$snp_peaks

  # Convert chrom to character so factors don't cause facet alignment issues
  snp_coverage[,    chrom := as.character(chrom)]
  chromosome_fits[, chrom := as.character(chrom)]

  p <- ggplot(snp_coverage, aes(x = pos_kb, y = n)) +
    geom_line(
      data      = chromosome_fits,
      aes(x = uniform_pos / 1000, y = uniform_fit),
      color     = "firebrick", linewidth = 0.6, alpha = 0.5
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

  # Highlight peak SNP points in blue
  if (!is.null(peaks_genomic) && nrow(peaks_genomic) > 0 &&
      !is.null(snp_peaks)     && nrow(snp_peaks) > 0) {

    snp_peaks_ch <- copy(snp_peaks)
    snp_peaks_ch[, chrom := as.character(chrom)]

    peak_highlight <- merge(
      snp_peaks_ch[, .(chrom, pos = snp_pos)],
      snp_coverage[,  .(chrom, pos, pos_kb, n)],
      by = c("chrom", "pos")
    )
    if (nrow(peak_highlight) > 0) {
      p <- p + geom_point(
        data        = peak_highlight,
        aes(x = pos_kb, y = n),
        color       = "black",
        fill        = "dodgerblue",
        size        = 2.5,
        shape       = 21,
        alpha       = 0.9,
        inherit.aes = FALSE
      )
    }
  }

  p
}
