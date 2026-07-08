# Changelog

All notable changes to ChimeraMapR are recorded here. Version numbers follow
`APP_VERSION` in `chimera_functions.R`, which is the single source of truth
read by `app.R` and `chimera_cli.R`.

## [Unreleased]

- Added an interstitial hemizygous-deletion call (`DELETION`, rendered as a
  `Δ` on the overview map). A HET-bounded fixed-allele LOH tract whose
  SNP-site read depth drops below `depth_drop` (default 0.60) of its higher
  HET flank is labeled a deletion of the missing homolog rather than a
  copy-neutral LOH tract. Uses the existing SNP-position coverage map only
  (no BAM/CNV/breakpoint modelling), so it is reported at `review`
  confidence. The higher of the two flanks is used as the diploid reference
  so a short depth-depressed (e.g. sub-telomeric) flank does not mask a real
  drop. Complements the existing terminal-deletion rule (R01).

## [0.7.0] - 2026-07-07

- Ported the beta-binomial EM + Viterbi LOH HMM to Rcpp/C++ for speed.
- Added aneuploidy detection: per-chromosome ploidy is estimated and shown
  as background shading on the overview plot.
- Extensive chaining-rule work on `loh_chain_analysis.R`: new
  `rule_one_sided_binary` (R11b, `GC_ONE_SIDED`) and `GC_UNRESOLVED` labels,
  gap-aware TEL-LOH chaining, exemptions so self-classifying peaks
  (gene_conversion / crossover / internal_crossover) aren't dropped or
  gated out by coverage/out-of-span checks, collapsing of interrupted or
  single-SNP-flight LOH regions, and disabling terminal-no-peak TCO calls.
- Added a `compound_binary` peak class (REF-ALT-REF-ALT / ALT-REF-ALT-REF
  patterns) that behaves like a binary peak for fusion/chaining but is
  reported separately in `--peak-list` output.
- Peak-fusion fixes: same-tract LOH pairs no longer incorrectly blocked
  from fusion, `peaks_bridge_independent_tracts()` fixed to use the peak
  window instead of the pad, and haplotype labeling moved into
  `compute_peak_pairs()` (simplifying the CLI).
- Added `del_rate_cutoff` as a SNP filtering parameter (default 10%,
  configurable under Advanced Settings).
- Added spanning-read counts (n-reads) to event/peak tables for support.
- Recombination Events tab now lists unresolved LOH regions/peaks as
  "Other Events" instead of silently dropping them; fixed three bugs
  affecting claimed/unclaimed event bookkeeping.
- Interface cleanup: help text on sidebar file inputs, dropped an unused
  legend on the Recombination Events tab, peaks below the qualifying
  height no longer shown in the displayed list, download-button layout
  fix.

## [0.6.0] - 2026-06-09

- Reworked `loh_chain_analysis.R` to handle long LOH tracts with crossovers,
  including CO-TERM / internal-CO classification at the single-peak level.
- Replaced Mclust GMM + RLE with a beta-binomial EM + Viterbi HMM for LOH
  calling; fixed the resulting namespace issue.
- Added LOH table display/output and LOH plots to each overview panel.
- Peak-fusion refactor (binary-only fusions, NA-peak filtering).
- Plot layout fixes.
- Added back the embedded gene-conversion rule and ASC bug fixes.
- Added depth calculation to help identify terminal deletions; peak-fusion
  fix.
- Fixed two event types over-extending their footprint on the chromosome.

## [0.5.0] - 2026-05-21

- First use of Whittaker-smoothed peak analysis.
- RDS download for the main graph and for peak/region plots (replacing
  oversized ggplot object downloads).
- Broke analysis functions out of `app.R` into `chimera_functions.R` so they
  can be sourced without Shiny (enables the CLI).
- Various interface and README updates.

Earlier versions (0.1.0–0.4.9) predate this changelog; see `git log` for
that history.
