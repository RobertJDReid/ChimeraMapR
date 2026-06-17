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
  if (requireNamespace("igraph", quietly = TRUE)) library(igraph)
})

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
  dt <- fread(path, header = FALSE)
  setnames(dt, 1:2, c("CHROM", "length"))
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
#  Classifies per-position allele balance into three haplotype states using a
#  custom EM algorithm (beta-binomial mixture) followed by Viterbi HMM
#  segmentation.  No external mixture-model packages are required.
#
#  States (keyed on ALT allele count k = n_alt):
#    HOM_REF     — binomial near 0 ALT  → REF_fixed (AC = 2)
#    DIP_HET_0.5 — beta-binomial at 0.5 → HET       (AC = 1)
#    HOM_ALT     — binomial near 1 ALT  → ALT_fixed (AC = 0)
#
#  Arguments
#    full_read_loh : data.table from run_chimera_analysis()$full_read_loh
#                   (MAPQ-only filtered; NO base-quality filter so that
#                   allele-balance estimates reflect the full pileup depth)
#    min_depth     : hard minimum reads per position; positions below this are
#                   excluded before fitting (default 5)
#    min_run_snps  : retained for API compatibility; not used by the HMM
#    trans_stay    : HMM self-transition probability — higher values produce
#                   fewer state changes / longer segments (default 0.999)
#
#  Returns a named list:
#    $snp_table    : per-SNP data.table
#                   (chrom, pos, n_ref, n_total, balance, AC, loh_state)
#    $loh_segments : collapsed contiguous-run table
#                   (chrom, start, end, length_bp, n_snps, AC, loh_state,
#                    balance_mean, balance_sd)
# -----------------------------------------------------------------------------
compute_loh_map <- function(full_read_loh,
                            min_depth    = 5L,
                            min_run_snps = 2L,
                            trans_stay   = 0.999,
                            warn_fn      = function(msg) message("WARNING: ", msg)) {

  # ── Internal model functions ───────────────────────────────────────────────

  dbinom_err <- function(k, N, p)
    dbinom(k, N, pmax(pmin(p, 1 - 1e-9), 1e-9))

  dbetabinom <- function(k, N, mu, theta) {
    alpha <- mu * theta
    beta_ <- (1 - mu) * theta
    out   <- exp(lchoose(N, k) + lbeta(k + alpha, N - k + beta_) - lbeta(alpha, beta_))
    out[!is.finite(out)] <- 0
    out
  }

  wll_theta <- function(theta, k, N, mu, w) {
    ll <- dbetabinom(k, N, mu, theta)
    ll[ll <= 0] <- 1e-300
    sum(w * log(ll))
  }

  fit_ab_mixture <- function(k, N, theta_init = 80, eps_init = 0.02,
                             max_iter = 200, tol = 1e-5) {
    labels <- c("HOM_REF", "DIP_HET_0.5", "HOM_ALT")
    pi     <- rep(1/3, 3)
    eps    <- eps_init
    theta  <- theta_init
    r      <- matrix(0, length(k), 3)
    for (iter in seq_len(max_iter)) {
      old    <- c(pi, eps, theta)
      r[, 1] <- pi[1] * dbinom_err(k, N, eps)
      r[, 2] <- pi[2] * dbetabinom(k, N, 0.5, theta)
      r[, 3] <- pi[3] * dbinom_err(k, N, 1 - eps)
      rs     <- rowSums(r); rs[rs == 0] <- 1e-300
      r      <- r / rs
      pi     <- colMeans(r)
      w_ref  <- r[, 1]; w_alt <- r[, 3]
      num    <- sum(w_ref * k / N) + sum(w_alt * (1 - k / N))
      den    <- sum(w_ref) + sum(w_alt)
      if (den > 0) eps <- pmax(pmin(num / den, 0.15), 1e-4)
      wj <- r[, 2]
      if (sum(wj) >= 1)
        theta <- optimize(wll_theta, c(5, 5000),
                          k = k, N = N, mu = 0.5, w = wj,
                          maximum = TRUE)$maximum
      if (max(abs(c(pi, eps, theta) - old) / (abs(old) + 1)) < tol) break
    }
    list(
      mixing_proportions = setNames(pi, labels),
      error_rate         = eps,
      theta              = theta,
      posterior          = `colnames<-`(r, labels),
      assignments        = labels[max.col(r)]
    )
  }

  viterbi_segment <- function(k, N, chrom, fit, trans_stay) {
    n_states  <- 3L
    trans_off <- (1 - trans_stay) / (n_states - 1)
    A_log     <- matrix(log(trans_off), n_states, n_states)
    diag(A_log) <- log(trans_stay)
    eps           <- fit$error_rate
    theta         <- fit$theta
    labels        <- names(fit$mixing_proportions)
    log_emit      <- matrix(0, length(k), n_states)
    log_emit[, 1] <- log(pmax(dbinom_err(k, N, eps),          1e-300))
    log_emit[, 2] <- log(pmax(dbetabinom(k, N, 0.5, theta),   1e-300))
    log_emit[, 3] <- log(pmax(dbinom_err(k, N, 1 - eps),      1e-300))
    log_pi        <- log(fit$mixing_proportions)
    assignments   <- character(length(k))
    for (chr in unique(chrom)) {
      idx <- which(chrom == chr)
      T_  <- length(idx)
      E   <- log_emit[idx, , drop = FALSE]
      delta      <- matrix(-Inf, T_, n_states)
      psi        <- matrix(0L,   T_, n_states)
      delta[1, ] <- log_pi + E[1, ]
      if (T_ > 1L) {
        for (t in 2L:T_) {
          scores     <- matrix(delta[t - 1L, ], n_states, n_states) + A_log
          psi[t, ]   <- apply(scores, 2, which.max)
          delta[t, ] <- apply(scores, 2, max) + E[t, ]
        }
        path     <- integer(T_)
        path[T_] <- which.max(delta[T_, ])
        for (t in (T_ - 1L):1L)
          path[t] <- psi[t + 1L, path[t + 1L]]
        assignments[idx] <- labels[path]
      } else {
        assignments[idx] <- labels[which.max(delta[1L, ])]
      }
    }
    assignments
  }

  # ── 1. Aggregate to per-position allele balance ────────────────────────────
  summary_dt <- full_read_loh[, .(
    n_ref   = sum(IS_REF, na.rm = TRUE),
    n_total = .N,
    n_alt   = .N - sum(IS_REF, na.rm = TRUE),
    balance = sum(IS_REF, na.rm = TRUE) / .N
  ), by = .(chrom, pos)]

  # Adaptive depth floor: median - 3*IQR (data-driven), clipped below at min_depth
  depth_cutoff <- max(min_depth,
                      floor(median(summary_dt$n_total) -
                              3 * IQR(summary_dt$n_total)),
                      na.rm = TRUE)
  summary_dt   <- summary_dt[n_total >= depth_cutoff]

  empty_segs <- data.table(
    chrom = character(), start = integer(), end = integer(),
    length_bp = integer(), n_snps = integer(), AC = integer(),
    loh_state = character(), balance_mean = numeric(), balance_sd = numeric()
  )

  if (nrow(summary_dt) < 10L) {
    warn_fn("Too few SNP positions with sufficient depth for LOH calling; returning NA.")
    summary_dt[, `:=`(AC = NA_integer_, loh_state = NA_character_)]
    return(list(
      snp_table    = summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)],
      loh_segments = empty_segs
    ))
  }

  # ── 2. Fit 3-state EM mixture ──────────────────────────────────────────────
  fit <- tryCatch(
    fit_ab_mixture(k = summary_dt$n_alt, N = summary_dt$n_total),
    error = function(e) {
      warn_fn(paste0("EM fit error: ", e$message))
      NULL
    }
  )
  if (is.null(fit)) {
    warn_fn("EM fit failed; LOH map will contain no called regions.")
    summary_dt[, `:=`(AC = NA_integer_, loh_state = NA_character_)]
    return(list(
      snp_table    = summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)],
      loh_segments = empty_segs
    ))
  }

  # ── 3. Viterbi HMM segmentation ───────────────────────────────────────────
  setorder(summary_dt, chrom, pos)
  hmm_states <- tryCatch(
    viterbi_segment(
      k          = summary_dt$n_alt,
      N          = summary_dt$n_total,
      chrom      = as.character(summary_dt$chrom),
      fit        = fit,
      trans_stay = trans_stay
    ),
    error = function(e) {
      warn_fn(paste0("Viterbi segmentation error: ", e$message))
      NULL
    }
  )
  if (is.null(hmm_states)) {
    summary_dt[, `:=`(AC = NA_integer_, loh_state = NA_character_)]
    return(list(
      snp_table    = summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)],
      loh_segments = empty_segs
    ))
  }

  # ── 4. Map HMM states → AC and loh_state ──────────────────────────────────
  summary_dt[, hmm_state := hmm_states]
  summary_dt[, AC := fcase(
    hmm_state == "HOM_REF",     2L,
    hmm_state == "DIP_HET_0.5", 1L,
    hmm_state == "HOM_ALT",     0L,
    default = NA_integer_
  )]
  summary_dt[, loh_state := fcase(
    AC == 0L, "ALT_fixed",
    AC == 1L, "HET",
    AC == 2L, "REF_fixed",
    default  = NA_character_
  )]
  summary_dt[, hmm_state := NULL]

  snp_table <- summary_dt[, .(chrom, pos, n_ref, n_total, balance, AC, loh_state)]

  # ── 5. Collapse into contiguous LOH segments ──────────────────────────────
  setorder(summary_dt, chrom, pos)
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

# ---------------------------------------------------------------------------
# classify_peak_haplotype()
#   Classifies the haplotype structure under a single detected peak by
#   examining the run-length pattern of SNP_call (REF / ALT / HET) within
#   the peak window.  The window is expanded outward (by SNP-position
#   counting, up to the median read length of reads in the peak) when
#   fewer than zone_min_snps + 1 positions exist on either side of snp_pos.
#
#   Arguments:
#     pk          : single-row data.table from snp_peaks
#     chr_name    : chromosome name (character)
#     rt_df       : full RLE-smoothed read data (data.table with chrom, pos,
#                   read_id, IS_REF columns)
#     touching_ids: read IDs that overlap the peak window (character vector)
#     zone_min_snps: minimum SNPs per side required (integer, default 2L)
#
#   Returns a named list:
#     label      : character label for the peak
#     seg_data   : data.table used for the segment bar (or NULL)
#     win_start  : actual window start used (bp, possibly expanded)
#     win_end    : actual window end used (bp, possibly expanded)
#     expanded   : logical — was the window expanded beyond peak_start/end?
# ---------------------------------------------------------------------------
classify_peak_haplotype <- function(pk, chr_name, rt_df, touching_ids,
                                    zone_min_snps = 2L) {

  snp_p    <- pk$snp_pos
  pk_start <- pk$peak_start
  pk_end   <- pk$peak_end
  min_each <- zone_min_snps + 1L        # SNPs required on each side

  # All positions available from reads touching this peak on this chromosome
  avail_pos <- sort(unique(
    rt_df[read_id %in% touching_ids & as.character(chrom) == chr_name, pos]
  ))

  if (length(avail_pos) == 0)
    return(list(label = "undefined", seg_data = NULL,
                win_start = pk_start, win_end = pk_end, expanded = FALSE))

  # Median read length for the expansion upper limit
  read_lengths <- rt_df[read_id %in% touching_ids & as.character(chrom) == chr_name,
                        .(len = max(pos) - min(pos)), by = read_id]$len
  max_expand   <- median(read_lengths, na.rm = TRUE)
  if (is.na(max_expand) || max_expand <= 0) max_expand <- 0

  # ── Helper: count positions strictly on each side of snp_pos ────────────
  count_left  <- function(s) sum(avail_pos >= s & avail_pos <  snp_p)
  count_right <- function(e) sum(avail_pos >  snp_p & avail_pos <= e)

  win_start <- pk_start
  win_end   <- pk_end
  expanded  <- FALSE

  # ── Expand left if needed ────────────────────────────────────────────────
  if (count_left(win_start) < min_each) {
    cands <- avail_pos[avail_pos < snp_p]
    limit <- snp_p - max_expand          # hard lower bound
    cands <- cands[cands >= limit]
    if (length(cands) >= min_each) {
      new_start <- cands[length(cands) - min_each + 1L]
      if (new_start < win_start) {
        win_start <- new_start
        expanded  <- TRUE
      }
    }
    # If still insufficient after expansion, we proceed; label will be "undefined"
  }

  # ── Expand right if needed ───────────────────────────────────────────────
  if (count_right(win_end) < min_each) {
    cands <- avail_pos[avail_pos > snp_p]
    limit <- snp_p + max_expand
    cands <- cands[cands <= limit]
    if (length(cands) >= min_each) {
      new_end <- cands[min_each]
      if (new_end > win_end) {
        win_end  <- new_end
        expanded <- TRUE
      }
    }
  }

  # ── Build seg_data from the (possibly expanded) window ──────────────────
  win_pos <- avail_pos[avail_pos >= win_start & avail_pos <= win_end]

  # Check we have at least zone_min_snps + 1 on each side now
  n_left  <- sum(win_pos <  snp_p)
  n_right <- sum(win_pos >  snp_p)
  if (n_left < min_each || n_right < min_each) {
    return(list(label = "undefined", seg_data = NULL,
                win_start = win_start, win_end = win_end, expanded = expanded))
  }

  # Aggregate IS_REF by position within the window to get SNP_call
  read_win_df <- rt_df[
    read_id %in% touching_ids &
      as.character(chrom) == chr_name &
      pos >= win_start & pos <= win_end
  ]
  if (nrow(read_win_df) == 0)
    return(list(label = "undefined", seg_data = NULL,
                win_start = win_start, win_end = win_end, expanded = expanded))

  peak_summary <- read_win_df[, .(
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
    , xmax := data.table::shift(xmin, type = "lead", fill = win_end / 1000)
  ][, SNP_call := factor(SNP_call, levels = c("ALT", "HET", "REF"))]

  # ── Classify by run pattern ──────────────────────────────────────────────
  runs <- as.character(seg_data$SNP_call)   # ordered sequence of run states
  n_runs <- length(runs)

  label <- if (n_runs == 1L) {
    if      (runs[1] == "HET")                           "internal_crossover"
    else                                                  "undefined"          # single REF or ALT

  } else if (n_runs == 2L) {
    if (runs[1] %in% c("REF","ALT") && runs[2] %in% c("REF","ALT") &&
        runs[1] != runs[2])                               "binary"
    else                                                  "undefined"

  } else if (n_runs == 3L) {
    pat <- paste(runs, collapse = "-")
    if      (pat %in% c("ALT-REF-ALT", "REF-ALT-REF")) "gene_conversion"
    else if (pat %in% c("HET-ALT-HET", "HET-REF-HET")) "internal_crossover"
    else                                                  "undefined"

  } else {
    "undefined"
  }

  list(
    label     = label,
    seg_data  = seg_data,
    win_start = win_start,
    win_end   = win_end,
    expanded  = expanded
  )
}

# ---------------------------------------------------------------------------
# classify_fused_peak_haplotype()
#   Variant of classify_peak_haplotype() for a fused peak group.
#   The "peak" is represented by fused_pos_bp (anchor), fused_start_bp, and
#   fused_end_bp (window).  touching_ids is the union of all reads from all
#   constituent sub-peaks.  Classification and window-expansion logic is
#   identical to the original; the only differences are which columns are read
#   from the fused group row and that the snp_pos anchor is the group mean.
#
#   Arguments:
#     fg           : single-row data.table from fused_peaks (one row per group,
#                    as produced by the group_rep summary — must contain
#                    fused_pos_bp, fused_start_bp, fused_end_bp)
#     chr_name     : chromosome name (character)
#     rt_df        : full RLE-smoothed read data (data.table)
#     touching_ids : union of reads from all constituent sub-peaks
#     zone_min_snps: minimum SNPs per side required (integer, default 2L)
#
#   Returns the same named list as classify_peak_haplotype():
#     label, seg_data, win_start, win_end, expanded
# ---------------------------------------------------------------------------
classify_fused_peak_haplotype <- function(fg, chr_name, rt_df, touching_ids,
                                          zone_min_snps = 2L) {

  snp_p    <- fg$fused_pos_bp
  pk_start <- fg$fused_start_bp
  pk_end   <- fg$fused_end_bp
  min_each <- zone_min_snps + 1L

  avail_pos <- sort(unique(
    rt_df[read_id %in% touching_ids & as.character(chrom) == chr_name, pos]
  ))

  if (length(avail_pos) == 0)
    return(list(label = "undefined", seg_data = NULL,
                win_start = pk_start, win_end = pk_end, expanded = FALSE))

  read_lengths <- rt_df[read_id %in% touching_ids & as.character(chrom) == chr_name,
                        .(len = max(pos) - min(pos)), by = read_id]$len
  max_expand   <- median(read_lengths, na.rm = TRUE)
  if (is.na(max_expand) || max_expand <= 0) max_expand <- 0

  count_left  <- function(s) sum(avail_pos >= s & avail_pos <  snp_p)
  count_right <- function(e) sum(avail_pos >  snp_p & avail_pos <= e)

  win_start <- pk_start
  win_end   <- pk_end
  expanded  <- FALSE

  if (count_left(win_start) < min_each) {
    cands <- avail_pos[avail_pos < snp_p]
    limit <- snp_p - max_expand
    cands <- cands[cands >= limit]
    if (length(cands) >= min_each) {
      new_start <- cands[length(cands) - min_each + 1L]
      if (new_start < win_start) { win_start <- new_start; expanded <- TRUE }
    }
  }

  if (count_right(win_end) < min_each) {
    cands <- avail_pos[avail_pos > snp_p]
    limit <- snp_p + max_expand
    cands <- cands[cands <= limit]
    if (length(cands) >= min_each) {
      new_end <- cands[min_each]
      if (new_end > win_end) { win_end <- new_end; expanded <- TRUE }
    }
  }

  win_pos <- avail_pos[avail_pos >= win_start & avail_pos <= win_end]
  n_left  <- sum(win_pos <  snp_p)
  n_right <- sum(win_pos >  snp_p)
  if (n_left < min_each || n_right < min_each)
    return(list(label = "undefined", seg_data = NULL,
                win_start = win_start, win_end = win_end, expanded = expanded))

  read_win_df <- rt_df[
    read_id %in% touching_ids &
      as.character(chrom) == chr_name &
      pos >= win_start & pos <= win_end
  ]
  if (nrow(read_win_df) == 0)
    return(list(label = "undefined", seg_data = NULL,
                win_start = win_start, win_end = win_end, expanded = expanded))

  peak_summary <- read_win_df[, .(
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
    , xmax := data.table::shift(xmin, type = "lead", fill = win_end / 1000)
  ][, SNP_call := factor(SNP_call, levels = c("ALT", "HET", "REF"))]

  runs   <- as.character(seg_data$SNP_call)
  n_runs <- length(runs)

  label <- if (n_runs == 1L) {
    if (runs[1] == "HET") "internal_crossover" else "undefined"
  } else if (n_runs == 2L) {
    if (runs[1] %in% c("REF","ALT") && runs[2] %in% c("REF","ALT") &&
        runs[1] != runs[2]) "binary" else "undefined"
  } else if (n_runs == 3L) {
    pat <- paste(runs, collapse = "-")
    if      (pat %in% c("ALT-REF-ALT", "REF-ALT-REF")) "gene_conversion"
    else if (pat %in% c("HET-ALT-HET", "HET-REF-HET")) "internal_crossover"
    else                                                  "undefined"
  } else {
    "undefined"
  }

  list(label = label, seg_data = seg_data,
       win_start = win_start, win_end = win_end, expanded = expanded)
}

# =============================================================================
#  PEAK FUSION FUNCTIONS
#  Moved from app.R so chimera_cli.R can call compute_peak_pairs() without
#  loading Shiny.
# =============================================================================

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

  classifiable <- state_df[!is.na(state_df$state_L) &
                           !is.na(state_df$state_R), ]

  no_middle <- all(is.na(state_df$state_M))

  if (nrow(classifiable) == 0) return("unresolvable")

  if (no_middle) {
    lr_patterns <- paste(classifiable$state_L, classifiable$state_R, sep = "-")
    unique_pats <- unique(lr_patterns)
    if (length(unique_pats) == 2) {
      sorted <- sort(unique_pats)
      if (sorted[1] == "ALT-REF" && sorted[2] == "REF-ALT")
        return("crossover")
    }
    return("ambiguous")
  }

  classifiable_3 <- classifiable[!is.na(classifiable$state_M), ]

  if (nrow(classifiable_3) == 0) {
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

  is_return <- function(p) {
    parts <- strsplit(p, "-")[[1]]
    parts[1] == parts[3] && parts[1] != parts[2]
  }
  all_return <- all(sapply(unique_pats, is_return))
  if (all_return && length(unique_pats) <= 2) return("gene_conversion")

  if (length(unique_pats) == 2) {
    sorted <- sort(unique_pats)
    if ((sorted[1] == "ALT-REF-REF" && sorted[2] == "REF-ALT-ALT"))
      return("crossover")
  }

  no_return <- all(!sapply(unique_pats, is_return))
  if (no_return) return("independent_events")

  return("ambiguous")
}

# ─────────────────────────────────────────────────────────────────────────────
#   FUSION HEURISTICS
# ─────────────────────────────────────────────────────────────────────────────

FUSION_HEURISTICS <- list(

  excluded_peak_classes  = c("gene_conversion", "internal_crossover"),
  eligible_peak_classes  = c("binary", "undefined"),
  unfusable_edge_types   = c("independent_events", "unresolvable"),
  auto_edge_types        = c("gene_conversion", "crossover"),
  supervised_edge_types  = c("gene_conversion", "crossover", "ambiguous"),

  edge_priority = c(gene_conversion    = 1,
                    crossover          = 2,
                    ambiguous          = 3,
                    independent_events = 4,
                    unresolvable       = 5),

  loh_min_depth    = 5L,
  loh_alt_max      = 0.15,
  loh_ref_min      = 0.85,
  loh_gap_min_frac = 0.5
)

# ---------------------------------------------------------------------------
# gap_has_loh()
#   Returns TRUE when any fixed-haplotype segment in loh_segments overlaps
#   the open interval (pos_start, pos_end) on chr_name.
# ---------------------------------------------------------------------------
gap_has_loh <- function(loh_segments, chr_name, pos_start, pos_end) {
  if (is.null(loh_segments) || nrow(loh_segments) == 0) return(FALSE)
  any(
    as.character(loh_segments$chrom) == chr_name &
    loh_segments$loh_state %in% c("REF_fixed", "ALT_fixed") &
    loh_segments$start < pos_end &
    loh_segments$end   > pos_start
  )
}

# gap_loh_state()
#   Returns the loh_state ("REF_fixed" or "ALT_fixed") of the first fixed
#   segment that overlaps (pos_start, pos_end), or NA if none.
# ---------------------------------------------------------------------------
gap_loh_state <- function(loh_segments, chr_name, pos_start, pos_end) {
  if (is.null(loh_segments) || nrow(loh_segments) == 0) return(NA_character_)
  idx <- which(
    as.character(loh_segments$chrom) == chr_name &
    loh_segments$loh_state %in% c("REF_fixed", "ALT_fixed") &
    loh_segments$start < pos_end &
    loh_segments$end   > pos_start
  )
  if (length(idx) == 0) return(NA_character_)
  loh_segments$loh_state[idx[1]]
}

# ---------------------------------------------------------------------------
# classify_loh_crossover_edge()
#   LOH-aware two-zone edge classifier for the LOH-crossover probing path.
#
#   state_list: data.frame with columns state_L and state_R (state_M is NA).
#   loh_state:  "REF_fixed" or "ALT_fixed" — the fixed state of the LOH in gap.
#   homog_frac: fraction of classifiable reads that must agree on a crossing
#               or same-state pattern to call CO / GC (default 0.80, matching
#               default_chain_params()$homog_frac). The "left_crossers" /
#               "right_crossers" read sets only test physical reach past the
#               far peak, not allele state along the way — a minority of
#               reads can have already reverted (e.g. a short, independent
#               conversion tract within the gap) before physically reaching
#               the far side. Requiring unanimous agreement among ALL
#               spanning reads let that minority noise flip real CO/GC
#               evidence to "ambiguous"; a majority-vote threshold tolerates it.
#
#   Key difference from classify_edge_type(): a single consistent L-R direction
#   (e.g., all ALT-REF with a REF_fixed LOH) is accepted as crossover evidence
#   because the LOH state itself is the complementary observation.  When
#   classify_edge_type() requires both ALT-REF and REF-ALT patterns, it misses
#   cases where only one crossing direction is observed due to read-length limits.
# ---------------------------------------------------------------------------
classify_loh_crossover_edge <- function(state_list, loh_state, homog_frac = 0.80) {
  loh_allele <- if (isTRUE(loh_state == "REF_fixed")) "REF" else "ALT"
  het_allele <- if (isTRUE(loh_state == "REF_fixed")) "ALT" else "REF"

  classifiable <- state_list[
    !is.na(state_list$state_L) & !is.na(state_list$state_R), ]
  n <- nrow(classifiable)
  if (n == 0) return("unresolvable")

  lr_patterns <- paste(classifiable$state_L, classifiable$state_R, sep = "-")

  co_l_pat   <- paste(het_allele, loh_allele, sep = "-")  # e.g. ALT-REF
  co_r_pat   <- paste(loh_allele, het_allele, sep = "-")  # e.g. REF-ALT
  same_pats  <- c(paste(loh_allele, loh_allele, sep = "-"),
                  paste(het_allele, het_allele, sep = "-"))

  n_crossing <- sum(lr_patterns %in% c(co_l_pat, co_r_pat))
  n_same     <- sum(lr_patterns %in% same_pats)

  if (n_crossing / n >= homog_frac) return("crossover")
  if (n_same     / n >= homog_frac) return("gene_conversion")

  "ambiguous"
}

# ---------------------------------------------------------------------------
# peak_is_fusion_eligible()
# ---------------------------------------------------------------------------
peak_is_fusion_eligible <- function(lbl) {
  is.na(lbl) ||
    lbl %in% FUSION_HEURISTICS$eligible_peak_classes
}

# ---------------------------------------------------------------------------
# decide_fusion_mode()
# ---------------------------------------------------------------------------
decide_fusion_mode <- function(edge_type, jaccard, jaccard_threshold,
                               loh_in_gap = FALSE) {

  if (edge_type %in% FUSION_HEURISTICS$unfusable_edge_types)
    return("none")

  if (jaccard >= jaccard_threshold &&
      edge_type %in% FUSION_HEURISTICS$auto_edge_types)
    return("automatic")

  if (loh_in_gap && edge_type %in% FUSION_HEURISTICS$auto_edge_types)
    return("automatic")

  if (jaccard > 0 &&
      edge_type %in% FUSION_HEURISTICS$supervised_edge_types)
    return("supervised")

  if (loh_in_gap && edge_type %in% FUSION_HEURISTICS$supervised_edge_types)
    return("supervised")

  "none"
}

# ---------------------------------------------------------------------------
# compute_peak_pairs()
#   Takes snp_peaks + rt_df + transition_pos and returns:
#     list(peak_pairs = data.table, fused_peaks = data.table)
#   Pass loh_segments (from compute_loh_map()$loh_segments) to enable LOH
#   promotion of supervised pairs; pass NULL to skip.
# ---------------------------------------------------------------------------
compute_peak_pairs <- function(snp_peaks,
                               rt_df,
                               transition_pos,
                               loh_segments       = NULL,
                               jaccard_threshold  = 0.20,
                               zone_min_snps      = 2L,
                               supervised_override = NULL,
                               homog_frac          = 0.80) {

  if (is.null(snp_peaks) || nrow(snp_peaks) == 0)
    return(list(peak_pairs = NULL, fused_peaks = NULL))

  peaks_dt <- copy(snp_peaks)
  peaks_dt <- peaks_dt[!is.na(snp_pos)]
  if (nrow(peaks_dt) == 0)
    return(list(peak_pairs = NULL, fused_peaks = NULL))
  peaks_dt[, peak_id := .I]
  peaks_dt[, chrom   := as.character(chrom)]

  all_pairs <- list()

  for (chr_name in unique(peaks_dt$chrom)) {

    chr_peaks <- peaks_dt[chrom == chr_name][order(snp_pos)]
    if (nrow(chr_peaks) < 2) next

    get_peak_reads <- function(pk) {
      transition_pos[
        as.character(chrom) == chr_name &
          pos >= pk$peak_start & pos <= pk$peak_end,
        unique(read_id)
      ]
    }

    read_sets <- lapply(seq_len(nrow(chr_peaks)), function(i) get_peak_reads(chr_peaks[i]))
    names(read_sets) <- chr_peaks$peak_id

    for (j in seq_len(nrow(chr_peaks) - 1L)) {

      pk_a <- chr_peaks[j]
      pk_b <- chr_peaks[j + 1L]

      # Peaks immediately flanking this pair on either side, if any. Reads can
      # be longer than the HET-FIX-HET (or HET-FIX) span being evaluated and
      # carry on into the NEXT recombination domain (e.g. chrI's adjacent
      # 91.3-96.5kb and 101.6-103.7kb LOH tracts, or an independent crossover
      # further out) — without a bound, classify_zone_state's majority vote
      # over an unbounded zone mixes that unrelated downstream signal into
      # this pair's haplotype call.
      prev_pk <- if (j >= 2L) chr_peaks[j - 1L] else NULL
      next_pk <- if (j + 2L <= nrow(chr_peaks)) chr_peaks[j + 2L] else NULL

      label_a <- if ("haplotype_label" %in% names(pk_a)) pk_a$haplotype_label else NA_character_
      label_b <- if ("haplotype_label" %in% names(pk_b)) pk_b$haplotype_label else NA_character_

      if (!peak_is_fusion_eligible(label_a) || !peak_is_fusion_eligible(label_b)) next

      reads_a <- read_sets[[as.character(pk_a$peak_id)]]
      reads_b <- read_sets[[as.character(pk_b$peak_id)]]

      all_pair_reads <- union(reads_a, reads_b)
      if (length(all_pair_reads) == 0) next

      read_lengths <- rt_df[read_id %in% all_pair_reads,
                            .(len = max(pos) - min(pos)), by = read_id]$len
      median_read_len <- median(read_lengths, na.rm = TRUE)
      if (is.na(median_read_len) || median_read_len <= 0) next

      if (is.na(pk_a$snp_pos) || is.na(pk_b$snp_pos)) next

      gap_bp <- pk_b$snp_pos - pk_a$snp_pos
      # Evaluate LOH in the gap before applying the read-length filter; the
      # LOH-crossover path relaxes the filter and uses a different spanning strategy.
      loh_in_gap <- if (!is.null(loh_segments)) {
        gap_has_loh(loh_segments, chr_name, pk_a$snp_pos, pk_b$snp_pos)
      } else {
        FALSE
      }

      if (!loh_in_gap && (is.na(gap_bp) || gap_bp > median_read_len)) next

      shared_reads <- intersect(reads_a, reads_b)

      # Activates whenever the LOH bridges the gap and no read is chimeric at
      # both flanking peaks at once. A read that crosses the LOH boundary only
      # once (the expected pattern for a true crossover) can never land in
      # both reads_a and reads_b, so intersect() is empty by construction for
      # that read population. Gating on gap_bp > median_read_len was
      # self-defeating: median_read_len is measured from the very reads that
      # might cross the gap, so good long crossing reads inflate the
      # threshold and prevent this branch from ever firing.
      is_loh_crossover_mode <- loh_in_gap && length(shared_reads) == 0

      # Clip the outer zones to the next peak on either side (exclusive of
      # that peak's own anchor position, which belongs to its own pair's
      # evaluation) so a long read's allele calls beyond this HET-FIX-HET
      # span never leak into this pair's classification.
      zone_l_start <- if (!is.null(prev_pk) && !is.na(prev_pk$snp_pos))
        prev_pk$snp_pos + 1 else -Inf
      zone_l_end   <- pk_a$snp_pos
      zone_m_start <- pk_a$snp_pos
      zone_m_end   <- pk_b$snp_pos
      zone_r_start <- pk_b$snp_pos
      zone_r_end   <- if (!is.null(next_pk) && !is.na(next_pk$snp_pos))
        next_pk$snp_pos - 1 else Inf

      n_spanning <- 0L
      edge_type  <- "unresolvable"
      n_shared   <- length(shared_reads)
      jaccard    <- 0

      if (is_loh_crossover_mode) {
        # Reads are chimeric on ONE side of the LOH only (they cross at a
        # boundary) so shared_reads (intersect(reads_a, reads_b)) is empty —
        # that's the condition that put us in this branch.
        # Find reads from each set that physically reach the opposite zone.
        chr_rt <- rt_df[as.character(chrom) == chr_name]

        left_crossers <- if (length(reads_a) > 0)
          chr_rt[read_id %in% reads_a,
                 .(reaches_right = any(pos >= zone_r_start)), by = read_id
          ][reaches_right == TRUE, read_id]
        else character(0)

        right_crossers <- if (length(reads_b) > 0)
          chr_rt[read_id %in% reads_b,
                 .(reaches_left = any(pos <= zone_l_end)), by = read_id
          ][reaches_left == TRUE, read_id]
        else character(0)

        spanning_ids <- unique(c(left_crossers, right_crossers))
        if (length(spanning_ids) == 0) next

        span_df <- chr_rt[read_id %in% spanning_ids]
        if (nrow(span_df) > 0) {
          state_list <- span_df[, {
            sL <- classify_zone_state(pos, IS_REF, zone_l_start, zone_l_end, zone_min_snps)
            sR <- classify_zone_state(pos, IS_REF, zone_r_start, zone_r_end, zone_min_snps)
            .(state_L = sL, state_M = NA_character_, state_R = sR)
          }, by = read_id]
          n_spanning <- nrow(state_list)
          # Use the LOH-aware classifier: a single consistent crossing direction
          # is accepted as crossover evidence when the LOH allele is the
          # complementary observation (classify_edge_type requires both).
          loh_state_in_gap <- gap_loh_state(
            loh_segments, chr_name, pk_a$snp_pos, pk_b$snp_pos)
          edge_type <- if (!is.na(loh_state_in_gap)) {
            classify_loh_crossover_edge(as.data.frame(state_list), loh_state_in_gap,
                                        homog_frac = homog_frac)
          } else {
            classify_edge_type(as.data.frame(state_list))
          }
        }

      } else {
        # Standard path: reads chimeric at both peaks span the gap naturally.
        n_union  <- length(union(reads_a, reads_b))
        jaccard  <- if (n_union > 0) n_shared / n_union else 0

        spanning_ids <- shared_reads

        if (length(spanning_ids) > 0) {
          span_df <- rt_df[read_id %in% spanning_ids &
                             as.character(chrom) == chr_name]

          if (nrow(span_df) > 0) {
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
      }

      pair_key <- paste0(chr_name, "_", pk_a$peak_id, "_", pk_b$peak_id)

      # LOH-crossover boundary peaks mark the opposite ends of a CO event and
      # must not be fused into a single peak record.
      fusion_mode <- if (is_loh_crossover_mode) {
        "none"
      } else {
        decide_fusion_mode(
          edge_type         = edge_type,
          jaccard           = jaccard,
          jaccard_threshold = jaccard_threshold,
          loh_in_gap        = loh_in_gap
        )
      }

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
        loh_in_gap   = loh_in_gap,
        edge_type    = edge_type,
        fusion_mode  = fusion_mode
      )
    }
  }

  pairs_dt <- if (length(all_pairs) > 0) {
    rbindlist(all_pairs)
  } else {
    data.table(
      pair_key           = character(), chrom        = character(),
      peak_id_a          = integer(),   peak_id_b    = integer(),
      snp_pos_a          = numeric(),   snp_pos_b    = numeric(),
      gap_bp             = numeric(),   median_read_len_bp = numeric(),
      n_reads_a          = integer(),   n_reads_b    = integer(),
      n_shared           = integer(),   n_spanning   = integer(),
      jaccard            = numeric(),   loh_in_gap   = logical(),
      edge_type          = character(), fusion_mode  = character()
    )
  }

  auto_edges <- pairs_dt[fusion_mode == "automatic"]

  peaks_dt[, fusion_group := peak_id]

  if (nrow(auto_edges) > 0) {
    if (!requireNamespace("igraph", quietly = TRUE))
      stop("Package 'igraph' is required for peak fusion. Install with: install.packages('igraph')")
    g <- igraph::graph_from_data_frame(
      d        = auto_edges[, .(from = peak_id_a, to = peak_id_b)],
      directed = FALSE,
      vertices = data.frame(name = peaks_dt$peak_id)
    )
    comps <- igraph::components(g)
    mem   <- comps$membership
    peaks_dt[, fusion_group := mem[as.character(peak_id)]]
  }

  fused_coords <- peaks_dt[, .(
    fused_pos_bp    = round(mean(snp_pos,    na.rm = TRUE)),
    fused_start_bp  = min(peak_start, na.rm = TRUE),
    fused_end_bp    = max(peak_end,   na.rm = TRUE),
    n_sub_peaks     = .N,
    constituent_ids = paste(peak_id, collapse = ",")
  ), by = fusion_group]

  peaks_dt <- merge(peaks_dt, fused_coords, by = "fusion_group", all.x = TRUE)

  get_peak_edge_info <- function(pid) {
    rel_pairs <- pairs_dt[peak_id_a == pid | peak_id_b == pid]
    if (nrow(rel_pairs) == 0)
      return(data.table(peak_id = pid, best_edge_type = NA_character_,
                        best_fusion_mode = NA_character_,
                        adjacent_pair_keys = NA_character_))
    rel_pairs[, priority := FUSION_HEURISTICS$edge_priority[edge_type]]
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
