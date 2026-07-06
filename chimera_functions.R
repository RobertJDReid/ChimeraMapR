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
  library(Rcpp)
  if (requireNamespace("igraph", quietly = TRUE)) library(igraph)
})

APP_VERSION <- "0.6.3"

# -----------------------------------------------------------------------------
#  Compile the beta-binomial EM + Viterbi HMM (src/loh_hmm.cpp), used by
#  compute_loh_map() below. Located relative to this file — not getwd() —
#  so this works whether sourced from app.R (cwd = app dir), chimera_cli.R
#  (may run via a PATH symlink from any directory), or an Rmd.
# -----------------------------------------------------------------------------
.chimera_script_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
  }
  getwd()
})
sourceCpp(file.path(.chimera_script_dir, "src", "loh_hmm.cpp"),
          cacheDir = file.path(.chimera_script_dir, "src", ".rcpp_cache"))

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
    del_rate_cutoff = 0.10,
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

  # ── SNP-level QC: drop SNPs with an excessive local deletion rate ──────────
  # A position where a large fraction of confidently-mapped (MAPQ-passing)
  # reads register a deletion rather than a base call usually sits in/near a
  # short repeat or homopolymer that destabilizes alignment for a subset of
  # reads; the reads that do get a base call there can carry a spurious
  # substitution, which masquerades as a single-SNP flip in an otherwise
  # fixed haplotype block. Excluding such SNPs here removes them from every
  # downstream read set (chimeric-read RLE, LOH allele balance, peak/
  # breakpoint detection) with no special-casing needed further down.
  del_stats <- read_data[
    mapq >= mapq_cutoff,
    .(n_total = .N, n_del = sum(is_del == 1L)),
    by = .(chrom, pos)
  ]
  del_stats[, del_frac := n_del / n_total]

  allele_data_used <- merge(
    allele_data_used, del_stats,
    by.x = c("CHROM", "POS"), by.y = c("chrom", "pos"), all.x = TRUE
  )
  n_high_del <- sum(!is.na(allele_data_used$del_frac) &
                       allele_data_used$del_frac > del_rate_cutoff)
  if (n_high_del > 0L) {
    warn_fn(paste0(
      n_high_del, " SNP position(s) excluded: local deletion rate exceeds ",
      round(del_rate_cutoff * 100), "% of confidently-mapped reads."
    ))
  }
  allele_data_used <- allele_data_used[is.na(del_frac) | del_frac <= del_rate_cutoff]
  allele_data_used[, c("n_total", "n_del", "del_frac") := NULL]

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
    full_read, allele_data_used,
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
    full_read_loh, allele_data_used,
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

  snp_coverage <- allele_data_used[, .(chrom = CHROM, pos = POS)]
  snp_coverage <- merge(snp_coverage, pos_count,
                        by = c("chrom", "pos"), all.x = TRUE)
  snp_coverage[is.na(n), n := 0L]
  snp_coverage[, chrom  := factor(chrom, levels = fasta_chr_order)]
  snp_coverage[, pos_kb := pos / 1000]

  # ── 5. Build per-chromosome lambda table ─────────────────────────────────────
  snp_counts_by_chr <- allele_data_used[, .N, by = CHROM]
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

    # (C) Deduplicate — keep fallback peaks not already found by findpeaks.
    # Use the findpeaks valley-to-valley extent (raw_peaks columns 3–4) to test
    # coverage: a fallback peak is redundant if its SNP position falls inside
    # any findpeaks peak's own extent.  The old approach used the fallback's
    # narrow ±5-step window to test for findpeaks *positions*, which failed when
    # the smoother placed its maximum >1 000 bp away from the raw-SNP cluster.
    if (!is.null(region_peaks_df) && nrow(region_peaks_df) > 0) {
      if (!is.null(raw_peaks) && nrow(raw_peaks) > 0) {
        fp_left  <- uniform_pos[raw_peaks[, 3]]
        fp_right <- uniform_pos[raw_peaks[, 4]]
      } else {
        fp_left  <- numeric(0)
        fp_right <- numeric(0)
      }
      region_peaks_df$is_new <- vapply(seq_len(nrow(region_peaks_df)), function(i) {
        fb_pos <- uniform_pos[region_peaks_df$peak_index[i]]
        !any(fp_left <= fb_pos & fb_pos <= fp_right)
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
      del_rate_cutoff = del_rate_cutoff,
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
  # EM fit + Viterbi decoding are implemented in src/loh_hmm.cpp (compiled at
  # source-time, see top of this file) — these are thin wrappers that keep
  # the same call signatures used below.

  loh_labels <- c("HOM_REF", "DIP_HET_0.5", "HOM_ALT")

  fit_ab_mixture <- function(k, N, theta_init = 80, eps_init = 0.02,
                             max_iter = 200, tol = 1e-5) {
    fit <- fit_ab_mixture_cpp(as.double(k), as.double(N),
                              theta_init, eps_init, max_iter, tol)
    list(
      mixing_proportions = setNames(fit$mixing_proportions, loh_labels),
      error_rate         = fit$error_rate,
      theta              = fit$theta
    )
  }

  viterbi_segment <- function(k, N, chrom, fit, trans_stay) {
    viterbi_segment_cpp(as.double(k), as.double(N), as.character(chrom),
                        as.double(fit$mixing_proportions), fit$error_rate, fit$theta,
                        trans_stay)
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
  # A lone single-SNP run sandwiched between two runs of the identical
  # flanking state (on both sides) is treated as a flicker, not a real
  # transition: a true recombination/gene-conversion tract is expected to be
  # supported by >= 2 consecutive SNPs. Such singletons are absorbed into
  # their flanking state before segments are built, mirroring the run-collapse
  # already applied to chimeric-read patterns in classify_peak_haplotype() /
  # classify_fused_peak_haplotype(). snp_table above is built from the
  # pre-collapse AC/loh_state, so per-SNP displays still show the raw,
  # uncollapsed observed pattern (including the flicker position itself).
  seg_dt <- copy(summary_dt)
  setorder(seg_dt, chrom, pos)
  repeat {
    seg_dt[, run_id := rleid(AC), by = chrom]
    run_ac <- seg_dt[, .(n = .N, ac = AC[1]), by = .(chrom, run_id)]
    setorder(run_ac, chrom, run_id)
    run_ac[, ac_prev := shift(ac, 1L),  by = chrom]
    run_ac[, ac_next := shift(ac, -1L), by = chrom]
    to_absorb <- run_ac[
      n == 1L & !is.na(ac_prev) & !is.na(ac_next) &
        ac_prev == ac_next & ac_prev != ac
    ]
    if (nrow(to_absorb) == 0L) break

    merge_map <- to_absorb[, .(chrom, run_id, new_ac = ac_prev)]
    seg_dt <- merge(seg_dt, merge_map, by = c("chrom", "run_id"), all.x = TRUE)
    seg_dt[!is.na(new_ac), AC := new_ac]
    seg_dt[, new_ac := NULL]
    seg_dt[, loh_state := fcase(
      AC == 0L, "ALT_fixed",
      AC == 1L, "HET",
      AC == 2L, "REF_fixed",
      default  = NA_character_
    )]
  }
  seg_dt[, run_id := rleid(AC), by = chrom]

  loh_segments <- seg_dt[, .(
    start        = min(pos),
    end          = max(pos),
    length_bp    = max(pos) - min(pos) + 1L,
    n_snps       = .N,
    AC           = unique(AC),
    loh_state    = unique(loh_state),
    balance_mean = mean(balance, na.rm = TRUE),
    balance_sd   = sd(balance,   na.rm = TRUE)
  ), by = .(chrom, run_id)]
  loh_segments[, run_id := NULL]

  list(
    snp_table    = snp_table,
    loh_segments = loh_segments
  )
}


# -----------------------------------------------------------------------------
#  compute_coverage_map()
#
#  Models per-position total read depth (NOT allele balance) to tell a
#  genuine sequencing-depth drop (e.g. a hemizygous terminal/arm deletion,
#  where one parental homolog's reads are simply absent) apart from depth
#  that is merely consistent with the rest of the chromosome.
#
#  This is deliberately a separate model from compute_loh_map(): n_total at
#  a SNP position carries no REF/ALT signal, so there's nothing to gain from
#  a balance-aware mixture here.  Architecture mirrors compute_loh_map() —
#  EM mixture fit genome-wide, then per-chromosome Viterbi segmentation —
#  but with a 2-state Gaussian mixture on log(depth) in place of the
#  3-state beta-binomial mixture on allele counts.
#
#  Arguments
#    full_read_loh : data.table from run_chimera_analysis()$full_read_loh
#                    (MAPQ-only filtered; NO base-quality filter — the same
#                    read population compute_loh_map() uses, so the two
#                    depth signals come from the same pileup)
#    min_depth     : hard minimum reads per position; positions below this
#                    are excluded before fitting (default 5)
#    trans_stay    : HMM self-transition probability (default 0.995; a true
#                    deletion boundary is a single sharp step, not a
#                    recombination tract, so this is less sticky than
#                    compute_loh_map()'s 0.999)
#
#  Returns a named list:
#    $coverage_table    : per-position data.table
#                         (chrom, pos, n_total, depth_state, depth_ratio)
#    $coverage_segments : collapsed contiguous-run table
#                         (chrom, start, end, length_bp, n_snps, depth_state,
#                          depth_mean, depth_ratio)
#    $baseline_depth    : median depth among NORMAL_DEPTH-called positions,
#                         genome-wide — the denominator for depth_ratio
# -----------------------------------------------------------------------------
compute_coverage_map <- function(full_read_loh,
                                 min_depth  = 5L,
                                 trans_stay = 0.995,
                                 warn_fn    = function(msg) message("WARNING: ", msg)) {

  # ── Internal model functions ───────────────────────────────────────────────

  fit_depth_mixture <- function(x, max_iter = 200, tol = 1e-6) {
    qs    <- quantile(x, c(0.10, 0.90), na.rm = TRUE)
    mu    <- c(qs[1], qs[2])
    sigma <- rep(max(sd(x) / 2, 0.05), 2)
    pi    <- c(0.5, 0.5)
    r     <- matrix(0, length(x), 2)
    for (iter in seq_len(max_iter)) {
      old    <- c(mu, sigma, pi)
      r[, 1] <- pi[1] * dnorm(x, mu[1], sigma[1])
      r[, 2] <- pi[2] * dnorm(x, mu[2], sigma[2])
      rs     <- rowSums(r); rs[rs == 0] <- 1e-300
      r      <- r / rs
      pi     <- colMeans(r)
      for (j in 1:2) {
        w  <- r[, j]
        sw <- sum(w)
        if (sw > 0) {
          mu[j]    <- sum(w * x) / sw
          sigma[j] <- sqrt(pmax(sum(w * (x - mu[j])^2) / sw, 1e-6))
        }
      }
      if (mu[1] > mu[2]) {     # keep state 1 = LOW, state 2 = NORMAL
        mu <- rev(mu); sigma <- rev(sigma); pi <- rev(pi); r <- r[, 2:1]
      }
      if (max(abs(c(mu, sigma, pi) - old) / (abs(old) + 1)) < tol) break
    }
    list(mu = mu, sigma = sigma, pi = pi,
        labels = c("LOW_DEPTH", "NORMAL_DEPTH"))
  }

  viterbi_segment_depth <- function(x, chrom, fit, trans_stay) {
    n_states  <- 2L
    trans_off <- (1 - trans_stay) / (n_states - 1)
    A_log     <- matrix(log(trans_off), n_states, n_states)
    diag(A_log) <- log(trans_stay)
    log_emit  <- cbind(
      dnorm(x, fit$mu[1], fit$sigma[1], log = TRUE),
      dnorm(x, fit$mu[2], fit$sigma[2], log = TRUE)
    )
    log_pi      <- log(fit$pi)
    labels      <- fit$labels
    assignments <- character(length(x))
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

  # ── 1. Aggregate to per-position total read depth ──────────────────────────
  cov_dt <- full_read_loh[, .(n_total = .N), by = .(chrom, pos)]
  setorder(cov_dt, chrom, pos)

  # Adaptive depth floor, same approach as compute_loh_map()'s depth_cutoff.
  depth_cutoff <- max(min_depth,
                      floor(median(cov_dt$n_total) - 3 * IQR(cov_dt$n_total)),
                      na.rm = TRUE)
  cov_dt <- cov_dt[n_total >= depth_cutoff]

  empty_segs <- data.table(
    chrom = character(), start = integer(), end = integer(),
    length_bp = integer(), n_snps = integer(),
    depth_state = character(), depth_mean = numeric(), depth_ratio = numeric()
  )

  if (nrow(cov_dt) < 10L) {
    warn_fn("Too few positions with sufficient depth for coverage modeling; returning NA.")
    cov_dt[, `:=`(depth_state = NA_character_, depth_ratio = NA_real_)]
    return(list(coverage_table = cov_dt, coverage_segments = empty_segs,
               baseline_depth = NA_real_))
  }

  # ── 2. Fit 2-state EM mixture on log(depth) ─────────────────────────────────
  fit <- tryCatch(
    fit_depth_mixture(log(cov_dt$n_total)),
    error = function(e) {
      warn_fn(paste0("Depth EM fit error: ", e$message))
      NULL
    }
  )
  if (is.null(fit)) {
    warn_fn("Depth EM fit failed; coverage map will contain no called regions.")
    cov_dt[, `:=`(depth_state = NA_character_, depth_ratio = NA_real_)]
    return(list(coverage_table = cov_dt, coverage_segments = empty_segs,
               baseline_depth = NA_real_))
  }

  # ── 3. Viterbi HMM segmentation (per chromosome) ────────────────────────────
  hmm_states <- tryCatch(
    viterbi_segment_depth(log(cov_dt$n_total), as.character(cov_dt$chrom),
                          fit, trans_stay),
    error = function(e) {
      warn_fn(paste0("Depth Viterbi segmentation error: ", e$message))
      NULL
    }
  )
  if (is.null(hmm_states)) {
    cov_dt[, `:=`(depth_state = NA_character_, depth_ratio = NA_real_)]
    return(list(coverage_table = cov_dt, coverage_segments = empty_segs,
               baseline_depth = NA_real_))
  }
  cov_dt[, depth_state := hmm_states]

  # Baseline = median depth among NORMAL_DEPTH-called positions, genome-wide.
  baseline <- cov_dt[depth_state == "NORMAL_DEPTH", median(n_total)]
  if (length(baseline) == 0 || is.na(baseline)) baseline <- median(cov_dt$n_total)
  cov_dt[, depth_ratio := round(n_total / baseline, 3)]

  # ── 4. Collapse into contiguous coverage segments ───────────────────────────
  setorder(cov_dt, chrom, pos)
  cov_dt[, COV_factor := rleid(depth_state), by = chrom]

  coverage_segments <- cov_dt[, .(
    start       = min(pos),
    end         = max(pos),
    length_bp   = max(pos) - min(pos) + 1L,
    n_snps      = .N,
    depth_state = unique(depth_state),
    depth_mean  = round(mean(n_total), 1),
    depth_ratio = round(mean(n_total) / baseline, 3)
  ), by = .(chrom, COV_factor)]
  coverage_segments[, COV_factor := NULL]

  list(
    coverage_table     = cov_dt[, .(chrom, pos, n_total, depth_state, depth_ratio)],
    coverage_segments  = coverage_segments,
    baseline_depth     = baseline
  )
}


# -----------------------------------------------------------------------------
#  get_chromosome_ploidy()
#
#  Estimates copy number per chromosome from per-position read depth.  Uses
#  the same full_read_loh input as compute_coverage_map() — MAPQ-only filtered
#  reads, no base-quality cutoff — so both functions sample the same pileup.
#
#  Algorithm: aggregate to per-position depth, take the genome-wide median as
#  the reference baseline for reference_ploidy copies, divide each
#  chromosome's median depth by that baseline, and round to the nearest
#  integer.  Median is used at both levels for robustness to local depth
#  variation and positional outliers.
#
#  Arguments
#    full_read_loh    : data.table from run_chimera_analysis()$full_read_loh
#    reference_ploidy : expected copy number of a normal chromosome (default 2)
#    min_depth        : hard minimum reads per position; positions below this
#                       are excluded before computing the baseline (default 5)
#
#  Returns a data.table with one row per chromosome:
#    chrom, n_positions, median_depth, depth_ratio, estimated_ploidy
# -----------------------------------------------------------------------------
get_chromosome_ploidy <- function(full_read_loh,
                                   reference_ploidy = 2L,
                                   min_depth        = 5L) {
  cov_dt <- full_read_loh[, .(n_total = .N), by = .(chrom, pos)]
  setorder(cov_dt, chrom, pos)

  depth_cutoff <- max(min_depth,
                      floor(median(cov_dt$n_total) - 3 * IQR(cov_dt$n_total)))
  cov_dt <- cov_dt[n_total >= depth_cutoff]

  baseline <- median(cov_dt$n_total)

  chr_ploidy <- cov_dt[, .(
    n_positions  = .N,
    median_depth = round(median(n_total), 1),
    depth_ratio  = round(median(n_total) / baseline, 3)
  ), by = chrom]

  chr_ploidy[, estimated_ploidy := as.integer(round(depth_ratio * reference_ploidy))]
  chr_ploidy[]
}


# -----------------------------------------------------------------------------
#  Recombination event symbols, overlaid on the LOH band at the bp midpoint
#  of each event's recorded start/end span. One symbol per event_class
#  recognised in EVENT_SYMBOL_MAP; event classes not listed are skipped.
# -----------------------------------------------------------------------------
EVENT_SYMBOL_MAP <- c(
  CO_GC                    = "✖", # HEAVY MULTIPLICATION X
  CO_GC_subres             = "✖", # HEAVY MULTIPLICATION X
  NCO_GC                   = "𝝤", # CAPTIAL OMICRON
  NCO_GC_subres            = "𝝤", # CAPTIAL OMICRON (no flanking LOH tract)
  NCO_GC_in_terminal       = "𝝤", # CAPTIAL OMICRON
  GC_UNRESOLVED            = "𝝤?", # confirmed GC-type tract (LOH flanked by two binary peaks); NCO/CO undetermined
  GC_ONE_SIDED             = "𝝤*", # one confirmed binary junction peak; other boundary has no peak at all
  CROSSOVER_NO_TRACT       = "✖", # HEAVY MULTIPLICATION X (crossover, tract below LOH resolution)
  DOUBLE_GC                = "𝝤𝝤", # two NCO gene conversions in one token
  CO_TERM                  = "TCO",
  CO_TERM_PROBABLE         = "TCO?",
  TCO_CAPTURED_TCO         = "2 X TCO",
  TERMINAL_DELETION        = "Δ",  # GREEK CAPITAL LETTER DELTA
  `AMBIGUOUS(low_coverage)` = "?"  # low spanning reads; shown as review marker
)

#' add_event_symbols()
#'   Adds a bold, centered symbol layer for recombination events to an
#'   existing chromosome-coverage plot. Symbols are centered horizontally on
#'   the event's bp midpoint (start/end from event_tbl, converted to Kb) and
#'   vertically in the middle of the LOH band, so they sit on the same line
#'   as the LOH region rectangles.
#'
#' @param p            ggplot to add the layer to
#' @param event_tbl    data.table with columns event_class, chrom, start, end
#'                      (e.g. results$event_table); NULL/empty is a no-op
#' @param band_ymin    numeric ymin of the LOH band (data units)
#' @param band_ymax    numeric ymax of the LOH band (data units)
#' @param chrom_filter optional chromosome to restrict to (character); NULL
#'                      keeps all rows, relying on facetting (overview plot)
#' @param size         text size passed to geom_text
add_event_symbols <- function(p, event_tbl, band_ymin, band_ymax,
                               chrom_filter = NULL, size = 5) {
  if (is.null(event_tbl) || nrow(event_tbl) == 0) return(p)

  ev <- copy(event_tbl)
  ev[, chrom := as.character(chrom)]
  if (!is.null(chrom_filter)) ev <- ev[chrom == chrom_filter]
  ev <- ev[event_class %in% names(EVENT_SYMBOL_MAP)]
  if (nrow(ev) == 0) return(p)

  ev[, x     := (start + end) / 2 / 1000]
  ev[, y     := (band_ymin + band_ymax) / 2]
  ev[, label := EVENT_SYMBOL_MAP[event_class]]

  p + geom_text(
    data        = ev,
    aes(x = x, y = y, label = label),
    fontface    = "bold",
    size        = size,
    inherit.aes = FALSE
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
#
#  When results$event_table is present (set after the chain event caller has
#  run), recombination events are overlaid as bold symbols centered on the
#  LOH band — see add_event_symbols() / EVENT_SYMBOL_MAP above.
# -----------------------------------------------------------------------------
build_overview_plot <- function(results) {

  snp_cov   <- copy(results$snp_coverage)
  fits      <- copy(results$chromosome_fits)
  peaks     <- results$peaks_genomic
  snp_peaks <- results$snp_peaks
  loh_segs_all <- results$loh_segments   # pre-built segment table; NULL if not yet run
  event_tbl    <- results$event_table    # from run_chain_analysis(); NULL if not yet run

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
  
  # ── Per-chromosome ploidy panel backgrounds ───────────────────────────────
  # 1N → #FAFBAC, 3N → #FEDBFF, 2N → white (no layer added).
  # fill is mapped through aes() (keyed on the display label) so these
  # backgrounds join the same fill legend as the LOH band colours below.
  PLOIDY_BG     <- c("1" = "#FCFCDA", "3" = "#F6E3FC")
  PLOIDY_LABELS <- c("1" = "<2N", "3" = ">2N")
  ploidy_dt <- if (!is.null(results$ploidy_map)) {
    #warning("non-null ploidy_map")
    results$ploidy_map
  } else if (!is.null(results$full_read_loh) && nrow(results$full_read_loh) > 0) {
    #warning("There is a full_read_loh value so use this to call get_chromosome_ploidy")
    get_chromosome_ploidy(results$full_read_loh)
  } else {
    #warning("complete fallback")
    data.table()
  }
  aneuploid <- if (nrow(ploidy_dt) > 0)
    ploidy_dt[estimated_ploidy %in% c(1L, 3L)] else data.table()
  has_ploidy_bg <- nrow(aneuploid) > 0
  # alpha < 1 so the panel grid lines (drawn beneath all geom layers by
  # ggplot2 regardless of code order) remain visible through the background
  ploidy_bg_layers <- lapply(seq_len(nrow(aneuploid)), function(i) {
    ploidy_key <- as.character(aneuploid$estimated_ploidy[i])
    geom_rect(
      data        = data.table(chrom       = as.character(aneuploid$chrom[i]),
                                ploidy_call = PLOIDY_LABELS[[ploidy_key]]),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = ploidy_call),
      color       = "black",
      alpha       = 0.5,
      inherit.aes = FALSE
    )
  })

  #browser()
  
  # ── Main coverage panel ────────────────────────────────────────────────────
  p_main <- ggplot(snp_cov, aes(x = pos_kb, y = n)) +
    ploidy_bg_layers +
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

  loh_fixed <- if (has_loh) copy(loh_segs_all[loh_state %in% c("REF_fixed", "ALT_fixed")]) else data.table()
  if (has_loh && nrow(loh_fixed) > 0) {
    # Convert bp coordinates to Kb and ensure chrom matches the coverage plot
    loh_fixed[, xmin  := start / 1000]
    loh_fixed[, xmax  := end   / 1000]
    loh_fixed[, chrom := as.character(chrom)]
    loh_fixed <- loh_fixed[chrom %in% unique(snp_cov$chrom)]
  } else {
    has_loh <- FALSE
  }
  loh_colours <- c(REF_fixed = "dodgerblue", ALT_fixed = "firebrick")

  # Height of the LOH band in data (read-count) units: 4% of the y ceiling.
  # Placed at ymin = 0 so it sits flush with the x-axis baseline in every facet.
  y_ceiling  <- max(30, max(snp_cov$n))
  loh_band_h <- y_ceiling * -0.15 # negative to put it below number line

  loh_labels <- c(
    REF_fixed = paste0(strain_ref), # will fix later
    ALT_fixed = paste0(strain_alt)
  )
  loh_caption <- paste0(
    "LOH band (bottom of each panel): blue\u202f=\u202f", strain_ref,
    ", red\u202f=\u202f", strain_alt
  )

  # -- Combined fill legend: ploidy backgrounds + LOH band colours --
  # Both sets of colours are mapped onto the same "fill" aesthetic (ploidy
  # backgrounds via aes(fill = ploidy_call) above) so they share one legend.
  fill_values <- character(0)
  fill_labels <- character(0)

  if (has_ploidy_bg) {
    ploidy_keys  <- unique(as.character(aneuploid$estimated_ploidy))
    ploidy_names <- PLOIDY_LABELS[ploidy_keys]
    fill_values  <- c(fill_values, setNames(PLOIDY_BG[ploidy_keys], ploidy_names))
    fill_labels  <- c(fill_labels, setNames(ploidy_names, ploidy_names))
  }

  if (has_loh) {
    fill_values <- c(fill_values, loh_colours)
    fill_labels <- c(fill_labels, loh_labels)

    p_main <- p_main +
      geom_rect(
        data        = loh_fixed,
        aes(xmin = xmin, xmax = xmax,
            ymin = 0,    ymax = loh_band_h,
            fill = loh_state),
        alpha       = 0.85,
        inherit.aes = FALSE
      )
  }

  if (length(fill_values) == 0) {
    return(p_main)   # nothing to add - return the plain overview
  }

  p_main <- p_main +
    scale_fill_manual(
      values = fill_values,
      labels = fill_labels,
      name   = "Panel key",
      drop   = FALSE
    ) +
    theme(
      legend.position  = "bottom",
      legend.title     = element_text(size = rel(1), face = "bold"),
      legend.text      = element_text(size = rel(1))
    )

  if (has_loh) {
    p_main <- p_main +
      labs(caption = loh_caption) +
      theme(plot.caption = element_text(size = rel(1), colour = "grey40"))

    p_main <- add_event_symbols(p_main, event_tbl, band_ymin = 0, band_ymax = loh_band_h)
  }

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

  # snp_n (from the earlier "detecting chimeric reads" pass, chimera_functions.R
  # transition_pos/pos_count) is the count of already-known-chimeric reads whose
  # own allele-run boundary sits exactly at this SNP -- real read support that
  # exists independently of whether the run-pattern logic below can classify a
  # clean switch type. Used as the n_support fallback whenever this function
  # can't resolve a label (e.g. a nearby second island widens the window into
  # an unrecognised >3-run pattern), so real evidence isn't discarded just
  # because the *type* of switch is ambiguous.
  snp_n_fallback <- suppressWarnings(as.integer(pk$snp_n %||% NA_integer_))

  # All positions available from reads touching this peak on this chromosome
  avail_pos <- sort(unique(
    rt_df[read_id %in% touching_ids & as.character(chrom) == chr_name, pos]
  ))

  if (length(avail_pos) == 0)
    return(list(label = "undefined", seg_data = NULL,
                win_start = pk_start, win_end = pk_end, expanded = FALSE,
                n_support = snp_n_fallback))

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
                win_start = win_start, win_end = win_end, expanded = expanded,
                n_support = snp_n_fallback))
  }

  # Aggregate IS_REF by position within the window to get SNP_call
  read_win_df <- rt_df[
    read_id %in% touching_ids &
      as.character(chrom) == chr_name &
      pos >= win_start & pos <= win_end
  ]
  if (nrow(read_win_df) == 0)
    return(list(label = "undefined", seg_data = NULL,
                win_start = win_start, win_end = win_end, expanded = expanded,
                n_support = snp_n_fallback))

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

  # ── Classify by run pattern, collapsing single-SNP noise runs ───────────
  # A run of exactly one SNP position is likely a spurious near-fixed site
  # rather than a real recombination signal. Remove such runs and re-derive
  # the run sequence before pattern matching; seg_data is returned unchanged
  # so the plot still shows the full observed pattern.
  noise_runs <- peak_summary[, .N, by = run][N == 1L, run]
  if (length(noise_runs) > 0L) {
    clean_ps <- peak_summary[!run %in% noise_runs]
    if (nrow(clean_ps) > 0L) {
      clean_ps[, run2 := data.table::rleid(SNP_call)]
      runs <- as.character(clean_ps[, .(SNP_call = SNP_call[1L]), by = run2][order(run2), SNP_call])
    } else {
      runs <- as.character(seg_data$SNP_call)
    }
  } else {
    runs <- as.character(seg_data$SNP_call)
  }
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

  # ── Count reads individually showing the switch pattern ─────────────────
  # The run pattern above pools all touching reads by position and derives
  # one consensus sequence — it never checks whether any single read itself
  # crosses the junction. For "binary" (one real boundary) and
  # "internal_crossover" (a narrow consensus run inside an otherwise
  # heterozygous, unphased flank) labels, count touching reads whose own
  # zone_L (win_start..snp_p) and zone_R (snp_p..win_end) allele calls
  # differ — i.e. reads that are themselves chimeric across this peak,
  # rather than just population-level noise.
  #
  # "gene_conversion" (REF-ALT-REF / ALT-REF-ALT, three runs) can't use the
  # same differ-test: its flanks share the same allele by construction, so
  # state_L == state_R for a genuine GC read and the switch count would
  # always be zero. Support for a GC call instead comes from how many reads
  # actually traverse both flanking zones — i.e. fully span the LOH region
  # the peak sits in — which is the read-level evidence that this is a real
  # chimeric junction rather than population-level noise.
  # "undefined" (n_runs not in 1..3, or a run pattern that doesn't match any
  # recognised switch type -- e.g. a second nearby island widening the window
  # into a compound multi-run signal) keeps the snp_n fallback set above; it's
  # real per-position read evidence even when the switch type can't be typed.
  n_support <- snp_n_fallback
  if (label %in% c("binary", "internal_crossover", "gene_conversion")) {
    read_states <- read_win_df[, {
      sL <- classify_zone_state(pos, IS_REF, win_start, snp_p, zone_min_snps)
      sR <- classify_zone_state(pos, IS_REF, snp_p, win_end, zone_min_snps)
      .(state_L = sL, state_R = sR)
    }, by = read_id]
    n_support <- if (label == "gene_conversion") {
      sum(!is.na(read_states$state_L) & !is.na(read_states$state_R))
    } else {
      sum(!is.na(read_states$state_L) & !is.na(read_states$state_R) &
          read_states$state_L != read_states$state_R)
    }
  }

  list(
    label     = label,
    seg_data  = seg_data,
    win_start = win_start,
    win_end   = win_end,
    expanded  = expanded,
    n_support = n_support
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

  noise_runs <- peak_summary[, .N, by = run][N == 1L, run]
  if (length(noise_runs) > 0L) {
    clean_ps <- peak_summary[!run %in% noise_runs]
    if (nrow(clean_ps) > 0L) {
      clean_ps[, run2 := data.table::rleid(SNP_call)]
      runs <- as.character(clean_ps[, .(SNP_call = SNP_call[1L]), by = run2][order(run2), SNP_call])
    } else {
      runs <- as.character(seg_data$SNP_call)
    }
  } else {
    runs <- as.character(seg_data$SNP_call)
  }
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
# peaks_bridge_independent_tracts()
#   TRUE when pk_a already marks the EXIT boundary of an existing fixed LOH
#   segment lying to its own left (away from pk_b), AND pk_b already marks
#   the ENTRY boundary of an existing fixed LOH segment lying to its own
#   right (away from pk_a) -- i.e. both peaks already belong to their own,
#   separately-bounded tracts before this pair is even considered.
#
#   This matters because classify_edge_type()'s 3-zone classifier can't
#   tell "one continuous excursion away and back" (a single GC tract) apart
#   from "two independent same-direction excursions sitting close together"
#   (two separate GC tracts with a short recovery in between): both produce
#   the same state_L == state_R != state_M pattern (e.g. ALT-REF-ALT) and
#   so both get edge_type = "gene_conversion". When that happens between
#   two peaks that each already terminate their own complete tract, fusing
#   them via jaccard/edge_type alone (decide_fusion_mode()'s "automatic"
#   path) would transitively merge two independent LOH events into one
#   fusion group through igraph::components() -- direct LOH evidence
#   bridging the gap itself (loh_in_gap) is required instead.
#
#   win_start_a/win_end_a: peak A's read window (defaults to pos_a for both,
#     giving an exact-match check). Using the full window avoids false
#     negatives when the LOH segment boundary sits a few hundred bp inside
#     the peak window rather than exactly at the SNP anchor.
# ---------------------------------------------------------------------------
peaks_bridge_independent_tracts <- function(loh_segments, chr_name,
                                            pos_a, pos_b,
                                            win_start_a = pos_a, win_end_a = pos_a,
                                            win_start_b = pos_b, win_end_b = pos_b) {
  if (is.null(loh_segments) || nrow(loh_segments) == 0) return(FALSE)
  if (is.na(pos_a) || is.na(pos_b)) return(FALSE)
  segs <- loh_segments[
    as.character(chrom) == chr_name & loh_state %in% c("REF_fixed", "ALT_fixed")
  ]
  if (nrow(segs) == 0) return(FALSE)
  a_exits  <- any(segs$end   >= win_start_a & segs$end   <= win_end_a & segs$end   < pos_b)
  b_enters <- any(segs$start >= win_start_b & segs$start <= win_end_b & segs$start > pos_a)
  a_exits && b_enters
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
                               loh_in_gap = FALSE,
                               bridges_independent_tracts = FALSE) {

  if (edge_type %in% FUSION_HEURISTICS$unfusable_edge_types)
    return("none")

  # Both peaks already terminate their OWN separate fixed tract, and no LOH
  # bridges this particular gap: jaccard/edge_type alone are not enough
  # evidence to fuse them into one group (see peaks_bridge_independent_tracts()).
  # Falls to "supervised" rather than "none" so a human can still review and
  # approve it — it just can't fuse silently/automatically.
  if (bridges_independent_tracts && !loh_in_gap &&
      edge_type %in% FUSION_HEURISTICS$auto_edge_types)
    return(if (jaccard > 0) "supervised" else "none")

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
# label_snp_peaks_haplotypes()
#   Classifies each peak's haplotype run-pattern (binary / gene_conversion /
#   internal_crossover / undefined) via classify_peak_haplotype() -- the
#   per-peak read evidence alone, before any fusion or event calling. Peaks
#   that already carry a non-NA haplotype_label are left untouched, so this
#   can be called again downstream (e.g. by compute_peak_pairs()) at zero
#   extra cost once the caller has already labeled peaks.
#
#   Returns snp_peaks with haplotype_label and n_read_support columns added
#   (or overwritten for previously-unlabeled rows); all other rows/columns
#   are preserved as-is.
# ---------------------------------------------------------------------------
label_snp_peaks_haplotypes <- function(snp_peaks, rt_df, transition_pos,
                                       zone_min_snps = 2L) {
  out <- copy(snp_peaks)
  if (!"haplotype_label" %in% names(out)) out[, haplotype_label := NA_character_]
  if (!"n_read_support"  %in% names(out)) out[, n_read_support  := NA_integer_]
  if (nrow(out) == 0) return(out)

  out[, .row_idx := .I]
  for (.chr in unique(as.character(out$chrom))) {
    .chr_p <- out[as.character(chrom) == .chr & !is.na(snp_pos)][order(snp_pos)]
    for (.pi in seq_len(nrow(.chr_p))) {
      .pk <- .chr_p[.pi]
      if (!is.na(.pk$haplotype_label)) next
      .touching <- transition_pos[
        as.character(chrom) == .chr &
          pos >= .pk$peak_start & pos <= .pk$peak_end,
        unique(read_id)
      ]
      if (length(.touching) == 0) {
        # No reads touch this peak's window at all -- same "no evidence"
        # bucket classify_peak_haplotype() itself returns when it can't
        # resolve a label, so callers filtering on "undefined" catch both.
        out[.row_idx == .pk$.row_idx, haplotype_label := "undefined"]
        next
      }
      .hap <- classify_peak_haplotype(.pk, .chr, rt_df, .touching, zone_min_snps)
      out[.row_idx == .pk$.row_idx, `:=`(
        haplotype_label = .hap$label,
        n_read_support  = .hap$n_support
      )]
    }
  }
  out[, .row_idx := NULL]
  out[]
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
    return(list(peak_pairs = NULL, fused_peaks = NULL, snp_peaks = snp_peaks))

  peaks_dt <- copy(snp_peaks)
  peaks_dt <- peaks_dt[!is.na(snp_pos)]
  if (nrow(peaks_dt) == 0)
    return(list(peak_pairs = NULL, fused_peaks = NULL, snp_peaks = snp_peaks))
  peaks_dt[, peak_id := .I]
  peaks_dt[, chrom   := as.character(chrom)]

  # ── 0. Label any unlabeled peaks ─────────────────────────────────────────
  # label_snp_peaks_haplotypes() skips peaks that already carry a
  # haplotype_label (e.g. from the app's plot-building loop or the CLI's
  # prior step), so callers that do their own labeling first don't pay twice.
  peaks_dt <- label_snp_peaks_haplotypes(peaks_dt, rt_df, transition_pos, zone_min_snps)

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

      label_a <- pk_a$haplotype_label
      label_b <- pk_b$haplotype_label

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

      # Do both peaks already terminate their OWN separate fixed tract (one
      # exiting tract A, the other entering tract B), with no LOH bridging
      # THIS gap? See peaks_bridge_independent_tracts() for why that
      # combination blocks automatic fusion further down.
      bridges_independent_tracts <- if (!is.null(loh_segments)) {
        peaks_bridge_independent_tracts(loh_segments, chr_name,
                                        pk_a$snp_pos,    pk_b$snp_pos,
                                        win_start_a = pk_a$peak_start,
                                        win_end_a   = pk_a$peak_end,
                                        win_start_b = pk_b$peak_start,
                                        win_end_b   = pk_b$peak_end)
      } else {
        FALSE
      }

      if (!loh_in_gap && (is.na(gap_bp) || gap_bp > median_read_len)) next

      shared_reads <- intersect(reads_a, reads_b)

      # Pre-compute Jaccard so it is available for both the loh_crossover_mode
      # condition below and the standard-path recording.
      n_union_early <- length(union(reads_a, reads_b))
      pre_jaccard   <- if (n_union_early > 0) length(shared_reads) / n_union_early else 0

      # Activates whenever the LOH bridges the gap AND read overlap between the
      # two peaks is below the Jaccard threshold. A read that crosses the LOH
      # boundary only once (the expected crossover pattern) can never land in
      # both reads_a and reads_b, so the ideal case has zero shared reads and
      # low Jaccard. The original == 0 condition failed for small LOH regions
      # (gap < median read length) where a handful of long reads happen to span
      # both boundaries and enter shared_reads — with just 1–2 spanning reads
      # the 3-zone standard path cannot reliably call CO vs. NCO and returns
      # "independent_events", poisoning best_edge_type on both peaks. Switching
      # to jaccard < threshold extends LOH crossover mode to cover these small
      # LOH cases without disturbing pairs where the two peaks genuinely share a
      # large fraction of reads (high Jaccard = not an LOH crossing pattern).
      is_loh_crossover_mode <- loh_in_gap && pre_jaccard < jaccard_threshold

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
        jaccard  <- pre_jaccard

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

      # Determine fusion_mode for this pair.
      # is_loh_crossover_mode fires when LOH fills the gap AND no read is
      # chimeric at both peaks (expected: each read crosses the boundary once).
      # That description fits two very different situations:
      #   (a) !bridges_independent_tracts: the two peaks are the entry/exit
      #       boundaries of the SAME LOH tract — fuse them automatically.
      #   (b)  bridges_independent_tracts: each peak already exits its own
      #       independent LOH tract into this gap — these are the two flanks of
      #       a crossover point and must not be merged.
      fusion_mode <- if (is_loh_crossover_mode && !bridges_independent_tracts) {
        "automatic"
      } else if (is_loh_crossover_mode) {
        "none"
      } else {
        decide_fusion_mode(
          edge_type                  = edge_type,
          jaccard                    = jaccard,
          jaccard_threshold          = jaccard_threshold,
          loh_in_gap                 = loh_in_gap,
          bridges_independent_tracts = bridges_independent_tracts
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
        bridges_independent_tracts = bridges_independent_tracts,
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
      bridges_independent_tracts = logical(),
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

  # n_read_support is intentionally left off this aggregate: it's a per-peak
  # value already present on peaks_dt (set above, or carried in from the
  # caller), and merging it back in from here would just collide with that
  # column. Self-classifying peak classes (gene_conversion / internal_crossover)
  # are excluded from automatic fusion entirely
  # (FUSION_HEURISTICS$excluded_peak_classes), so their fusion_group is always
  # a singleton -- the per-peak value already equals what a group aggregate
  # would produce.
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
    lbl       <- peaks_dt[peak_id == pid, haplotype_label]
    own_label <- if (length(lbl) && !is.na(lbl)) lbl else NA_character_

    if (nrow(rel_pairs) == 0) {
      # Singleton — no pairs evaluated.  Use the peak's own haplotype_label so
      # callers can show something informative (e.g. "internal_crossover") instead
      # of NA.
      return(data.table(peak_id            = pid,
                        best_edge_type     = own_label,
                        best_fusion_mode   = NA_character_,
                        adjacent_pair_keys = NA_character_))
    }
    # Prefer auto-fusion pairs over supervised/none pairs.  Without this, a
    # supervised gene_conversion pair (edge_priority=1) can outrank an automatic
    # crossover pair (edge_priority=2), making the peak appear supervised even
    # though it is already committed to an auto-fusion group.
    auto_rel <- rel_pairs[fusion_mode == "automatic"]
    pool     <- if (nrow(auto_rel) > 0) auto_rel else rel_pairs
    pool[, priority := FUSION_HEURISTICS$edge_priority[edge_type]]
    best <- pool[which.min(priority)]

    # A pair's edge_type only stands in for this peak's OWN self-classification
    # when the pair actually fused it (fusion_mode == "automatic"). An
    # evaluated-but-unconfirmed ("supervised"/"none") pair describes this
    # peak's relationship to its *neighbor*, not this peak's own identity —
    # promoting it here would let the neighbor's verdict masquerade as a
    # self-classifying label (e.g. a genuinely "binary" peak reading as
    # "gene_conversion"), which lets single-peak rules fire on borrowed
    # evidence instead of leaving the peak available to rules that combine
    # both flanking peaks.
    best_edge <- if (nrow(auto_rel) > 0) best$edge_type else own_label

    data.table(
      peak_id            = pid,
      best_edge_type     = best_edge,
      best_fusion_mode   = best$fusion_mode,
      adjacent_pair_keys = paste(rel_pairs$pair_key, collapse = ";")
    )
  }

  edge_info <- rbindlist(lapply(peaks_dt$peak_id, get_peak_edge_info))
  peaks_dt  <- merge(peaks_dt, edge_info, by = "peak_id", all.x = TRUE)

  # Write haplotype labels back to a clean copy of snp_peaks so callers don't
  # need to re-run classify_peak_haplotype and the CLI can drop its standalone
  # labeling step.
  snp_peaks_out <- copy(snp_peaks)
  if (!"haplotype_label" %in% names(snp_peaks_out))
    snp_peaks_out[, haplotype_label := NA_character_]
  if (!"n_read_support" %in% names(snp_peaks_out))
    snp_peaks_out[, n_read_support := NA_integer_]
  label_map <- peaks_dt[, .(snp_pos, chrom, haplotype_label, n_read_support)]
  snp_peaks_out[label_map, on = .(snp_pos, chrom),
                `:=`(haplotype_label = i.haplotype_label,
                     n_read_support  = i.n_read_support)]

  list(
    peak_pairs  = pairs_dt,
    fused_peaks = peaks_dt,
    snp_peaks   = snp_peaks_out
  )
}
