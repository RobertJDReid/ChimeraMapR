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
# mclust is intentionally NOT loaded here at the top level.
# It exports em() which conflicts with data.table's internal gzip detection
# and causes a fatal error at startup.  compute_loh_map() checks for it via
# requireNamespace() and calls it with mclust:: qualification so it is loaded
# lazily only when the LOH analysis step runs.

APP_VERSION <- "0.6.0"


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

  # Snapshot of all filtered reads BEFORE chimeric-only subsetting.
  # Returned as $full_read so the caller can build a population-level LOH map
  # that covers fixed-haplotype regions even where no chimeric reads exist.
  full_read_all <- copy(full_read)

  # ── LOH-specific read set: MAPQ + is_del only, NO base-quality filter ────────
  # Allele balance for LOH calling benefits from maximum SNP depth; the base
  # quality filter applied to the chimeric-read pipeline is intentionally
  # omitted here so that bulk allele frequencies reflect the full pileup.
  full_read_loh <- read_data[
    mapq   >= mapq_cutoff &
    is_del == 0
  ]
  full_read_loh <- merge(
    full_read_loh, allele_data,
    by.x  = c("chrom", "pos"),
    by.y  = c("CHROM", "POS"),
    all.x = TRUE
  )
  full_read_loh[, IS_REF := call == REF]
  full_read_loh[, ALLELE := fcase(
    call == REF, "REF",
    call == ALT, "ALT",
    default = "OTHER"
  )]
  full_read_loh <- full_read_loh[ALLELE != "OTHER"]
  full_read_loh <- full_read_loh[, .(chrom, pos, read_id, IS_REF, ALLELE)]

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
    full_read         = full_read_all,   # complete filtered reads (pre-chimeric filter); used for LOH
    full_read_loh     = full_read_loh,   # MAPQ-only filtered reads (no base-qual); used for LOH allele balance
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
#  compute_loh_map()
#
#  Replaces the simple threshold-based compute_loh_map_full().  Uses a
#  three-component Gaussian mixture model (Mclust) to classify per-position
#  allele balance into AC = 0 (ALT-fixed / REF homozygous), AC = 1 (het), or
#  AC = 2 (REF-fixed / ALT homozygous), then applies run-length smoothing to
#  remove isolated short runs that are likely noise.
#
#  Arguments
#    full_read_loh : data.table from run_chimera_analysis()$full_read_loh
#                   (MAPQ-only filtered; NO base-quality filter so that
#                   allele-balance estimates reflect the full pileup depth)
#    min_depth     : minimum reads at a position to attempt a call (default 5)
#    min_run_snps  : short runs <= this length are subjected to strict balance
#                   re-checking; those that fail are removed before final
#                   run re-labelling (default 2, mirrors Rmd logic)
#
#  Returns a data.table with columns:
#    chrom, pos, n_ref, n_total, balance, AC, loh_state
#  where loh_state is one of "ALT_fixed" | "REF_fixed" | "HET" | NA
#  (NA = fewer than min_depth reads at that position or Mclust could not fit).
# -----------------------------------------------------------------------------
compute_loh_map <- function(full_read_loh,
                            min_depth    = 5L,
                            min_run_snps = 2L,
                            warn_fn      = function(msg) message("WARNING: ", msg)) {

  if (!requireNamespace("mclust", quietly = TRUE))
    stop("Package 'mclust' is required for LOH calling. Install with: install.packages('mclust')")
  
  # Attach mclust temporarily so Mclust()'s internal eval(mc, parent.frame())
  # can resolve MclustBIC. We detach on exit to avoid the em()/data.table conflict.
  # library() would permanently attach it; this scopes it to this function call only.
  already_attached <- "package:mclust" %in% search()
  if (!already_attached) {
    suppressPackageStartupMessages(
      attachNamespace("mclust")
    )
    on.exit({
      try(detach("package:mclust", unload = FALSE, character.only = TRUE), silent = TRUE)
    }, add = TRUE)
  }
  
  # ── 1. Aggregate to per-position allele balance ────────────────────────────
  summary_dt <- full_read_loh[, .(
    n_ref   = sum(IS_REF, na.rm = TRUE),
    n_total = .N,
    balance = sum(IS_REF, na.rm = TRUE) / .N
  ), by = .(chrom, pos)]

  # Drop positions with insufficient depth — Mclust can't use them
  summary_dt <- summary_dt[n_total >= min_depth]

  if (nrow(summary_dt) < 10L) {
    warn_fn("Too few SNP positions with sufficient depth for Mclust LOH calling; returning NA.")
    summary_dt[, `:=`(AC = NA_integer_, loh_state = NA_character_)]
    return(summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)])
  }
  
  # ── 2. Fit Gaussian mixture model on allele balance ────────────────────────
  # Allow Mclust to select the best number of components (G = 1:3).
  # Forcing G = 3 can fail silently on samples that are predominantly HET
  # (only 1–2 components present), returning NULL and causing no LOH to be
  # called.  Letting Mclust pick means it may return G = 1 or G = 2;
  # the cluster-ranking step below handles any G gracefully.
  fit <- withCallingHandlers(
    tryCatch(
      mclust::Mclust(summary_dt$balance, G = 1:3, verbose = FALSE),
      error = function(e) { warn_fn(paste0("Mclust error: ", e$message)); NULL }
    ),
    warning = function(w) {
      warn_fn(paste0("Mclust warning: ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  
  if (is.null(fit)) {
    warn_fn("Mclust returned NULL; LOH map will contain no called regions.")
    summary_dt[, `:=`(AC = NA_integer_, loh_state = NA_character_)]
    return(summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)])
  }
  
  summary_dt[, cluster := fit$classification]

  # Rank clusters by mean balance: lowest mean → AC=0, middle → AC=1, highest → AC=2.
  # When G < 3 is selected, only 1 or 2 AC values will be present; positions in
  # the single cluster of a G=1 fit are all called HET (AC=1).
  n_clusters <- length(fit$parameters$mean)
  cluster_means <- data.table(
    cluster      = seq_len(n_clusters),
    mean_balance = fit$parameters$mean
  )
  setorder(cluster_means, mean_balance)

  if (n_clusters == 1L) {
    # Only one component — all positions are HET; no LOH regions
    cluster_means[, AC := 1L]
  } else if (n_clusters == 2L) {
    # Two components: assign AC=0 (low balance) and AC=2 (high balance)
    # This handles a sample with only REF-LOH and HET, or ALT-LOH and HET,
    # or (rarely) pure REF-LOH vs ALT-LOH.
    # Use the midpoint: if both means are near 0/1 extremes → AC 0 and 2;
    # if one is near 0.5 → AC 0/2 and 1 respectively.
    mid <- mean(cluster_means$mean_balance)
    if (mid < 0.25) {
      cluster_means[, AC := c(0L, 1L)]   # ALT-fixed + HET
    } else if (mid > 0.75) {
      cluster_means[, AC := c(1L, 2L)]   # HET + REF-fixed
    } else {
      cluster_means[, AC := c(0L, 2L)]   # ALT-fixed + REF-fixed (no HET)
    }
  } else {
    # Three components: standard AC=0/1/2 ranking by mean balance
    cluster_means[, AC := 0L:2L]
  }

  summary_dt[cluster_means, AC := i.AC, on = "cluster"]
  summary_dt[, cluster := NULL]

  # ── 3. Run-length smoothing — remove noisy short runs ─────────────────────
  # First pass: label runs per chromosome
  setorder(summary_dt, chrom, pos)
  summary_dt[, LOH_factor := rleid(AC), by = chrom]
  summary_dt[, AC_runs    := .N,        by = .(chrom, LOH_factor)]

  # For runs of <= min_run_snps, keep only rows that pass strict balance gates.
  # Positions that fail are dropped so that formerly-split runs can merge.
  strict_pass <- function(ac, bal) {
    (ac == 0L & bal < 0.02) |
    (ac == 1L & bal > 0.40 & bal < 0.60) |
    (ac == 2L & bal > 0.98)
  }
  summary_dt <- summary_dt[
    !(AC_runs <= min_run_snps & !strict_pass(AC, balance))
  ]

  # Re-label after deletions so formerly-interrupted runs merge
  summary_dt[, LOH_factor := rleid(AC), by = chrom]
  summary_dt[, AC_runs    := .N,        by = .(chrom, LOH_factor)]

  # ── 4. Map AC to loh_state ──────────────────────────────────────────────────
  summary_dt[, loh_state := fcase(
    AC == 0L, "ALT_fixed",
    AC == 2L, "REF_fixed",
    AC == 1L, "HET",
    default  = NA_character_
  )]

  snp_table <- summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)]

  # ── 5. Collapse per-SNP rows into contiguous LOH segments ──────────────────
  # Mirrors the LOH_table aggregation from the Rmd workflow.  Adjacent
  # same-state positions are merged into one row per run, giving start/end
  # coordinates that downstream plotting and gap_has_loh() can use directly
  # without re-running RLE compression on every render.
  setorder(summary_dt, chrom, pos)
  # LOH_factor was set in step 3; re-derive on the final filtered set so
  # indices are contiguous after any row deletions.
  summary_dt[, LOH_factor := rleid(AC), by = chrom]

  loh_segments <- summary_dt[, .(
    start        = min(pos),
    end          = max(pos),
    length_bp    = max(pos) - min(pos) + 1L,
    n_snps       = .N,
    AC           = unique(AC),
    loh_state    = unique(loh_state),
    balance_mean = mean(balance, na.rm = TRUE),
    balance_sd   = sd(balance,   na.rm = TRUE)
  ), by = .(chrom, LOH_factor)]
  loh_segments[, LOH_factor := NULL]

  list(
    snp_table    = snp_table,
    loh_segments = loh_segments
  )
}


# -----------------------------------------------------------------------------
#  Overview plot builder
#
#  Expressed as a plain function that returns a single ggplot object.
#
#  results : the list returned by run_chimera_analysis()
#
#  Strip labels use facet_grid(chrom ~ ., switch = "y") with angle = 0 so
#  that chromosome names read horizontally on the left side of each panel.
#
#  The X axis upper limit is data-driven: the maximum chromosome length from
#  results$chr_span (in bp, converted to Kb) so the scale is always sized to
#  the longest chromosome in the dataset.
#
#  When results$loh_map is present and contains REF_fixed / ALT_fixed regions,
#  a thin LOH band (4% of the y-axis height) is drawn at the bottom of each
#  chromosome facet panel via geom_rect(), ensuring perfect per-chromosome
#  alignment regardless of the number of chromosomes.  When no loh_map is
#  available the function returns the plain overview plot unchanged.
# -----------------------------------------------------------------------------
build_overview_plot <- function(results) {

  snp_cov   <- copy(results$snp_coverage)
  fits      <- copy(results$chromosome_fits)
  peaks     <- results$peaks_genomic
  snp_peaks <- results$snp_peaks
  loh_segs_all <- results$loh_segments   # pre-built segment table; NULL if not yet run

  # Strain display names (set by UI inputs; fall back to generic labels)
  strain_ref <- if (!is.null(results$strain_ref) && nzchar(results$strain_ref))
    results$strain_ref else "REF"
  strain_alt <- if (!is.null(results$strain_alt) && nzchar(results$strain_alt))
    results$strain_alt else "ALT"

  # ── Data-driven X axis upper limit ────────────────────────────────────────
  # Use the longest chromosome length from chr_span (bp → Kb).
  # Falls back to the observed maximum SNP position if chr_span is absent.
  if (!is.null(results$chr_span) && nrow(results$chr_span) > 0) {
    x_max_kb <- ceiling(max(results$chr_span$length, na.rm = TRUE) / 1000)
  } else {
    x_max_kb <- ceiling(max(snp_cov$pos_kb, na.rm = TRUE))
  }
  # Round up to the nearest 100 Kb for clean minor-break grid lines
  x_max_kb <- ceiling(x_max_kb / 25) * 25

  # Convert chrom to character so factors don't cause facet alignment issues
  snp_cov[, chrom := as.character(chrom)]
  fits[,    chrom := as.character(chrom)]

  # ── Main coverage panel ────────────────────────────────────────────────────
  p_main <- ggplot(snp_cov, aes(x = pos_kb, y = n)) +
    geom_line(
      data  = fits,
      aes(x = uniform_pos / 1000, y = uniform_fit),
      color = "firebrick", linewidth = 0.6, alpha = 0.7
    ) +
    geom_point(color = "black", alpha = 0.5, size = 0.5, shape = 21) +
    scale_x_continuous(
      limits       = c(0, x_max_kb),
      minor_breaks = seq(0, x_max_kb, 100)
    ) +
    xlab("Position (Kbp)") +
    ylab("Number of Reads") +
    ylim(NA, max(30, max(snp_cov$n))) +
    facet_grid(chrom ~ ., switch = "y") +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_line(linewidth = 0.05, color = "black"),
      panel.grid.major.x = element_line(linewidth = 0.05, color = "red"),
      strip.background   = element_blank(),
      strip.placement    = "outside",
      axis.text          = element_text(size = rel(1.2)),
      axis.title         = element_text(size = rel(1.5)),
      strip.text.y       = element_text(size = rel(1.1), angle = 0, hjust = 0, face = "bold")
    )

  # Highlight peak SNP points in blue
  if (!is.null(snp_peaks) && nrow(snp_peaks) > 0) {
    snp_cov_char   <- copy(snp_cov)   # already character from above
    snp_peaks_char <- copy(snp_peaks)
    snp_peaks_char[, chrom := as.character(chrom)]

    peak_highlight <- merge(
      snp_peaks_char[!is.na(snp_pos), .(chrom, pos = snp_pos)],
      snp_cov_char[, .(chrom, pos, pos_kb, n)],
      by = c("chrom", "pos")
    )
    if (nrow(peak_highlight) > 0) {
      p_main <- p_main +
        geom_point(
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

  # ── LOH overlay ──────────────────────────────────────────────────────────────────────
  # Rendered as a geom_rect band at the base of each facet panel so that
  # LOH regions are always aligned with their chromosome row without relying
  # on patchwork to synchronize two separate faceted plots.
  #
  # loh_segments is a pre-collapsed segment table from compute_loh_map():
  # one row per contiguous run of the same loh_state.  No inline RLE needed.
  has_loh <- !is.null(loh_segs_all) && nrow(loh_segs_all) > 0 &&
             any(loh_segs_all$loh_state %in% c("REF_fixed", "ALT_fixed"), na.rm = TRUE)

  if (!has_loh) {
    return(p_main)   # nothing to add — return the plain overview
  }

  loh_fixed <- copy(loh_segs_all[loh_state %in% c("REF_fixed", "ALT_fixed")])
  if (nrow(loh_fixed) == 0) return(p_main)

  # Convert bp coordinates to Kb and ensure chrom matches the coverage plot
  loh_fixed[, xmin  := start / 1000]
  loh_fixed[, xmax  := end   / 1000]
  loh_fixed[, chrom := as.character(chrom)]
  loh_fixed <- loh_fixed[chrom %in% unique(snp_cov$chrom)]
  loh_colours <- c(REF_fixed = "dodgerblue", ALT_fixed = "firebrick")

  # Height of the LOH band in data (read-count) units: 4% of the y ceiling.
  # Placed at ymin = 0 so it sits flush with the x-axis baseline in every facet.
  y_ceiling  <- max(30, max(snp_cov$n))
  loh_band_h <- y_ceiling * -0.10 # negative to put it below number line

  loh_labels <- c(
    REF_fixed = paste0(strain_ref, " (blue)"),
    ALT_fixed = paste0(strain_alt, " (red)")
  )
  loh_caption <- paste0(
    "LOH band (bottom of each panel): blue\u202f=\u202f", strain_ref,
    ", red\u202f=\u202f", strain_alt
  )

  p_main <- p_main +
    geom_rect(
      data        = loh_fixed,
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
      plot.caption     = element_text(size = rel(1), colour = "grey40"),
      legend.position  = "bottom",
      legend.title     = element_text(size = rel(1), face = "bold"),
      legend.text      = element_text(size = rel(1))
    )

  p_main
}
