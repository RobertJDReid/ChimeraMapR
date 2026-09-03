# Changelog

All notable changes to ChimeraMapR are recorded here. Version numbers follow
`APP_VERSION` in `chimera_functions.R`, which is the single source of truth
read by `app.R` and `chimera_cli.R`.

## [0.8.15] - 2026-09-03

### R12 corroborates its crossover claim at the tract

- `rule_loh_crossover` (R12) asserted `CO_GC` purely from a peak pair's
  `edge_type = "crossover"`, with no check that any read actually carried a
  homolog across the tract. That verdict comes from
  `classify_loh_crossover_edge()` over the two peaks' outer zones, and those
  zones are bounded by the *neighbouring peaks* — not by heterozygosity. When a
  zone lands inside another fixed tract every read reads the same allele there,
  so the L-R "crossing" it scores is the LOH map read back to itself with no
  phase content at all.
- On RAD5_15 `S288C_chrI` that produced a high-confidence `CO_GC` over the
  18.2 kb REF tract at 75,326–93,531. The pair backing it (74,978 ↔ 93,531)
  rests on **one** read, 73,740–94,351, whose left zone (71,797–74,978) is nine
  SNPs of the adjacent 2.7 kb ALT tract plus the single HET SNP at 74,978. The
  read is ALT across the ALT tract because the ALT tract is ALT, REF at 74,978,
  and REF onward — an uninformative read, scored as `ALT-REF` and therefore as a
  crossover.
- R12 now applies the same tract corroboration R10's `gene_conversion` branch
  has carried since 0.8.14, from `annotate_tract_read_support()` — whose
  flanking zones come from the adjacent chain tokens, so an uncallable flank
  yields no spanning reads instead of a fabricated state:
  `n_tract_switch > 0` → `CO_GC`, else `n_tract_return > 0` → `NCO_GC` (the
  reads say the tract closed, whatever the pair said), else `GC_UNRESOLVED`.
  The tract is real either way, so an unobserved outcome downgrades the call
  rather than discarding the event. Absent counts (no `full_read_loh`) leave the
  peak-based call untouched.
- The chrI call becomes `GC_UNRESOLVED` at "review": its left flank is the
  single SNP at 74,978, below `ZONE_CALL_HEURISTICS$min_evidence_snps`, so the
  tract has 0 junction-spanning reads and no CO/NCO outcome is observable. No
  other RAD5_15 event changed (17 events, 14 high / 3 review).

## [0.8.14] - 2026-09-02

### Gene-conversion outcomes are decided by reads crossing the tract

- R10 (`rule_peak_direct`) treated a `gene_conversion` peak as self-classifying
  and promoted whatever F token it was attached to straight to `NCO_GC`. But
  that label is read off a REF-ALT-REF (or ALT-REF-ALT) return pattern inside
  the *peak's* window, and `classify_peak_haplotype()` scores its `n_support`
  as the reads traversing that window — neither is tied to the tract the peak
  ends up attached to. Since `.attach_peaks()` binds any peak whose `snp_pos`
  falls within `peak_pad_bp` (200 bp) of a token boundary, a peak describing a
  short island could be bound to a tract many kb wide whose far junction no
  read ever reached.
- New `count_tract_junction_reads()` (`chimera_functions.R`) measures what the
  outcome call actually turns on. For a tract and its two flanking zones it
  classifies every read crossing **both** junctions: *return* (same homolog
  both flanks, and it is the one the tract is not — NCO), *switch* (a different
  homolog each side — CO), or *uninformative* (flanks already match the tract
  state: the unconverted homolog running straight through, which inside a fixed
  tract is most of the reads present).
- `annotate_tract_read_support()` attaches `n_tract_return` / `n_tract_switch`
  / `n_tract_spanning` to every F token, after `canonicalise()` since merging
  changes the spans the counts refer to. `full_read_loh` now flows into
  `run_chain_analysis()` and is passed by both `chimera_cli.R` and `app.R`; the
  counts are serialised into the chain step CSVs. When it is not supplied the
  fields are absent and every consumer falls back to the previous peak-based
  behaviour.
- R10 gates on the measurement: `n_tract_return > 0` → `NCO_GC`, else
  `n_tract_switch > 0` → `CO_GC`, else `GC_UNRESOLVED` (𝝤°) — the tract is real
  either way, so an unobserved outcome downgrades the call rather than
  discarding the event. `crossover` and `internal_crossover` peaks are exempt:
  a crossover signature is a single-junction observation by construction, with
  no return that has to span anything.
- The flanking zones come from the adjacent chain tokens rather than a fixed bp
  window either side: on RAD5_09 `S288C_chrII` the left flank of the
  91,510–98,709 tract *is* a 439 bp restored-HET island with another fixed
  tract immediately beyond it, and a fixed window would reach past the island
  into that neighbour and read the wrong state.
- A read counts as a *return* only if it positively carries the tract's allele
  at the in-tract SNPs it covers — its own calls, checked against the tract
  state from the LOH map, at `homog_frac` agreement. Both halves of that matter.
  Using the read's calls rather than a `classify_zone_state()` majority vote
  keeps short tracts observable: a tract can hold fewer SNPs than the zone
  caller needs, which would otherwise force a genuine 4 bp conversion on
  RAD5_07 `S288C_chrII` to "outcome unknown". Requiring the match to be
  *positive* is what keeps the count honest: treating an uncallable tract zone
  as permissive scores every opposite-homolog read that merely overlaps the
  tract as a return — roughly half the spanning reads by construction. On EV_10
  that alone manufactured eight `NCO_GC` calls, four of them "high confidence",
  on tracts with no chimeric peak anywhere near them; every supporting read was
  the *unconverted* homolog reading straight through a coverage dip. And the
  threshold is a fraction, not unanimity, because over a multi-kb tract one
  miscalled base is expected — the RAD5_09 `S288C_chrII` 91,510–98,709 call
  rests on a single read that is ALT at 35 of its 36 in-tract SNPs.

### New rule R11c — a tract's own reads can call it

- `rule_tract_read_evidence` fires on an F token whose reads settle the
  outcome, classifying from `n_tract_return` / `n_tract_switch` with
  `homog_frac` requiring the informative reads to agree; a mixed tract reports
  `AMBIGUOUS(mixed_tract_reads)` with the split in its notes.
- This is the fallback for a tract no peak-based rule can reach. R10 needs a
  self-classifying peak bound to the token, R11 needs a *binary* peak at each
  junction, R11b needs one binary peak and genuinely no peak opposite. A tract
  whose junction peak happens to be typed `gene_conversion` — because it
  describes a short restored-het island between two conversion patches —
  satisfies none of them, and fell through to `reconcile()` as an uncalled LOH
  even when reads plainly crossed it.
- It claims the LOH token but **no peak**. Peak exclusivity ("peaks are shared
  observations that back a single event") is right for calls whose evidence is
  the peak; here the evidence is the reads, and the junction peak may
  legitimately be the shared boundary between this tract and its neighbour,
  already consumed by the neighbour's call. Claiming the token still keeps one
  call per tract, while leaving the peak unclaimed lets a shared junction back
  both adjacent events. Placed last in `MOTIF_RULES` so the peak-based rules
  keep priority through the resolver's score.

### Confidence on read-supported calls

- Events carrying `support_kind = "tract_reads"` are gated on the new
  `tract_read_high_min` (4): one to three returning/switching reads → "review",
  more than three → "high". Tract-spanning reads are a far scarcer observation
  than peak-window traversals — crossing both junctions of a tract of any size
  outruns most reads — so a handful establishes the event without settling it.
- Applied ahead of the per-class rule but only for classes that could otherwise
  reach "high", so an `AMBIGUOUS` call is never promoted by having plenty of
  reads behind its ambiguity.

### Two ways one event was being reported twice

- `reconcile()` promoted unclaimed peaks by iterating `fused_peaks`, which
  holds one row **per peak** while `fused_pos_bp` is the fusion *group*
  aggregate (the mean of the constituent peaks' `snp_pos`). A conversion tract
  bracketed by two binary junction peaks is one fused pair sharing one
  position, so it emitted two identical `*_subres` events differing only in
  `n_support` - the two junctions' separate read counts. RAD5_09
  `S288C_chrVII` 435,397 (n=41) and 436,829 (n=37) both reported at 436,113,
  their midpoint; likewise `S288C_chrVIII` 418,921 (n=10) and 419,699 (n=11) at
  419,310. Four such pairs across the 43-sample set (also RAD5_07
  `S288C_chrIX` 202,877 and RAD5_13 `S288C_chrVII` 109,694).
- `compute_peak_pairs()` deliberately leaves `n_read_support` off the group
  aggregate, on the reasoning that self-classifying peak classes are excluded
  from fusion and so always form singletons. That holds for a peak whose own
  `haplotype_label` is `gene_conversion`; it does not hold here, where the
  *pair's* `edge_type` is `gene_conversion` while both *peaks* are `binary`,
  which is fusion-eligible. `reconcile()` now collapses the peak source to one
  row per position, combining support as the **minimum** across the group: an
  NCO is the claim that the tract closed, witnessed only by a read crossing
  both junctions, so the thinner junction bounds the achievable support.
- Separately, `.fired_claim_keys()` derived a candidate's claim keys from the
  tokens it consumed, but `.make_event()` can re-anchor an event off its
  matched token onto a peak's phase island. Two rules claiming different tokens
  could therefore describe the same interval with no collision visible to the
  resolver. On R187E_13 `S288C_chrI`, R10 matched the F token 28,705-122,023
  through the `internal_crossover` peak at 129,863, re-anchored onto that peak's
  island, and landed on 130,264-130,305 - the exact tract R11c had claimed on
  its own reads. Both committed, giving one crossover twice (n=39 and n=42).
  Claim keys now include each event's post-anchor footprint, so the two contend
  and `MOTIF_RULES` priority settles it in favour of the peak-backed call.

### Effect on the 43-sample RAD5-OE set

- 449 -> 451 events. Six new calls from R11c (three `NCO_GC`, two `CO_GC`, one
  `AMBIGUOUS(mixed_tract_reads)`), spread over five samples, and four surplus
  rows removed by the fused-pair collapse. No event was lost. Every new call
  sits on a tract the LOH map calls unambiguously fixed (`balance_mean`
  0.006-0.016 or 0.992-1.000, with a single 3-SNP tract at 0.113).
- One class change: RAD5_09 `S288C_chrII` 69,647-90,558 (20.9 kb), previously
  `NCO_GC` at high confidence with `n_support = 20`, now `GC_UNRESOLVED`. The
  `gene_conversion` peak at 90,610 describes the 439 bp island - its window is
  89,087-93,687 - but sat 52 bp inside the pad at that tract's right edge. No
  read crosses both of its junctions: the one read spanning the interval is the
  unconverted homolog, pure ALT throughout. The adjacent 91,510-98,709 tract,
  previously uncalled, is now `NCO_GC` on the single read that does return
  across it.
- Eleven `NCO_GC` calls move from "high" to "review", each resting on one to
  three returning reads. Several are large (18.5 kb, 16.7 kb, 20.8 kb), where
  thin junction-spanning support is expected rather than suspicious.
- Final split: 348 high / 103 review. No sample now reports the same interval
  twice.

### Known limitation exposed by this work

- R11c calls a tract from its reads alone, with no requirement that a chimeric
  peak sit anywhere near it. That is deliberate - it is the whole point of the
  rule - but it means the call inherits whatever the LOH map decided, and the
  LOH map is not equally trustworthy in every sample. On EV_10 all fifteen
  "fixed" segments have `balance_mean` between 0.72 and 0.87; none approaches 0
  or 1, and the sample carries no chimeric peaks at all. A sample with no
  genuine LOH gives the beta-binomial EM no fixed-state population to fit, and
  its fixed component appears to settle on the most allele-skewed positions
  available - typically local coverage dips. Compare RAD5_09, cleanly bimodal
  at 0.000-0.021 and 0.998-0.999.
- Nothing upstream guards against this: the local-deletion-rate SNP filter
  targets positions where reads register deletions, and Rd
  (`rule_interstitial_deletion`) needs `flank_depth_ratio < 0.60` where these
  dips sit at 0.62-0.85. The positive in-tract match now keeps such tracts from
  producing events (EV_10 goes from eight calls to none), so the problem is
  latent rather than active, but the LOH map itself is unchanged and no rule
  currently checks a tract's `balance_mean` before trusting it.

## [0.8.13] - 2026-08-05

- **Event labels moved below the LOH band**, on both the genome-wide overview
  map and the app's per-chromosome coverage tab. Symbols were drawn centred
  inside the coloured band, so on a short LOH region the glyph covered the
  region it was annotating and on a long one it read as text printed over a
  solid colour. They now hang just beneath the band, clear of the rectangles.
- To keep the vertical footprint close to what it was, the band itself is 20%
  shorter on both plots: 0.15 → 0.12 of the y ceiling on the overview, 0.10 →
  0.08 on the per-chromosome tab.
- The label lane below the band is reserved explicitly by lowering the y-axis
  floor. `geom_text` extents are physical rather than data units, so they
  neither grow the scale nor shrink with it; without a reserved lane the
  symbols are simply clipped at the panel edge. The lane is sized from the
  glyph height against the estimated panel height — on the overview that comes
  from the chromosome count (matching how the app and CLI size the export),
  and the symbol size drops from 5 to 4 above two chromosomes; the
  per-chromosome tab is a single panel and keeps size 5.
- `add_event_symbols()` gained `place = c("center", "below")` and `gap`. The
  default remains `"center"`, so any other caller is unaffected.
- **New plot palette**, collected into a single `CHIMERA_COLOURS` constant in
  `chimera_functions.R` so the overview map, the per-chromosome coverage tab
  and the read-level haplotype views can no longer drift apart. REF is
  `#3AA2FC` (was `dodgerblue`, `#1E90FF`) and ALT is `#A51301` (was
  `firebrick`, `#B22222`), applied to the LOH bands, the haplotype segment
  strips and the per-SNP read colours alike. SNP peak highlights are now
  `#FFC300` (they previously shared the REF blue, which made a peak marker
  hard to tell from an LOH call at a glance) and the uniform-coverage fit line
  is `#6D28D9` (it previously shared the ALT red). The two no longer share a
  colour: a peak marker and the fit it sits on are different claims. Amber and
  mid-green were both tried and dropped for the fit line — amber is ~9 ΔE from
  the ALT red even for normal colour vision and mid-green is within ~14 ΔE of
  the HET grey, while violet adds no confusable pair. Neither hue carries an
  alpha channel any more; the geoms set their own alpha, and an `#RRGGBB3F`
  line drawn at `alpha = 0.7` was compounding to ~17% opacity and disappearing
  into the pale 1N panel background. The overview fit line also goes from
  `linewidth` 0.6 to 0.8, matching the per-chromosome tab. Remaining colours —
  the ploidy panel backgrounds, HET grey, gridlines, boundary and sub-peak
  markers — are unchanged.
- Read-level SNP colours no longer come from `scale_color_viridis_d(option =
  "turbo", …)`; `IS_REF` is mapped through an explicit `scale_color_manual()`
  instead, with `NA` calls still drawn in `grey50` as before.

## [0.8.12] - 2026-08-05

- Interstitial hemizygous deletions are detected again. Two independent gates
  were suppressing them, and both had to be lifted before a known 5 kb deletion
  on `test_data/SYNv1_delta_test.csv.gz` (chrXI 500.2–505.3 kb) would call.
- **The SNP-level deletion-rate QC was erasing its own evidence.** A hemizygous
  deletion leaves one homolog absent, so at every SNP it spans ~50% of
  MAPQ-passing reads register a deletion instead of a base call — 36–46% across
  the SYNv1 block. That is the same per-SNP signature `del_rate_cutoff` (0.10)
  exists to catch in repeats/homopolymers, so all 38 informative SNPs were
  dropped, the region collapsed into a single 128 kb HET segment, and the
  chain caller saw an unremarkable `G` gap instead of a fixed tract. Nothing
  was drawn and nothing was called. The two cases differ in *distribution*, not
  in `del_frac`: an artifact hits a random subset of reads at isolated
  positions, whereas a deletion hits the same reads at every position across a
  contiguous block. New `find_hemizygous_del_blocks()` applies that test —
  a run of `del_block_min_snps` (5) consecutive over-cutoff SNPs spanning
  `del_block_min_bp` (500) with `del_block_read_coherence` (0.50) of its
  deletion calls coming from reads deleted across the block — and exempts
  qualifying blocks from exclusion. The SYNv1 block scores 0.99 coherence
  (21 reads deleted at all 38 positions; 7 stray single-position calls).
- **The deletion rule's depth test was biased upward by long flanks.**
  `.interstitial_flank_depth_ratio()` averaged real depth over the *entire*
  flanking token. Flank tokens run 50–100 kb and coverage drifts across them,
  so the 68 kb left flank averaged 51.2 against ~59 immediately beside the
  tract, putting the ratio at 0.62 and missing `depth_drop` (0.60) even once
  the SNPs were retained. It now measures each flank over the coverage segment
  *abutting* the tract — `compute_coverage_map()` has already located the
  changepoints, so that segment is the local diploid baseline by construction
  and no window width has to be guessed. SYNv1 chrXI: 31.7 against the
  abutting 491.0–499.8 kb segment (mean 54.5) → 0.58, a clear call. The `min()`
  of the two flanks is unchanged, so the sub-telomeric false-positive guard
  from 0.8.7 (RAD5_3 chrIII) still holds.
- Effect on results: `SYNv1_delta_test` now reports `DELETION S288C_chrXI
  500217–505314` (5.097 kb, `flank_depth_ratio=0.58`, review confidence)
  alongside its four pre-existing events, and the Δ symbol appears on the
  overview. Regression-checked on `test_data` chrI, chrII, chrIII, chrXI,
  chrXIII and chrXV: events tables byte-identical to 0.8.11, and the
  high-`del_frac` SNPs on those chromosomes (2–29 per chromosome, all isolated
  artifacts) are still excluded — no block qualified for the exemption on any
  real dataset.
- New parameters, exposed on both interfaces: `--del-block-min-snps`,
  `--del-block-min-bp`, `--del-block-read-coherence` in `chimera_cli.R`, and
  Min SNPs / Min Read Coherence inputs in the app sidebar. Setting
  `del_block_min_snps` arbitrarily high restores 0.8.11 behaviour exactly.
  `run_chimera_analysis()` also returns `$del_blocks` listing any exempted
  blocks.

## [0.8.11] - 2026-08-04

- Event symbols no longer depend on the LOH band. `build_overview_plot()` and
  the app's per-chromosome coverage plot both called `add_event_symbols()` from
  *inside* their LOH block, so a sample with no `REF_fixed`/`ALT_fixed` segment
  had every symbol silently dropped — the events were still called, still in
  the events table, and still visible on the individual peak plots, but the
  overview showed a bare coverage trace. `build_overview_plot()` compounded it
  with an early `return()` when the fill legend was empty, which on a 2N
  chromosome with no LOH bailed out before the symbol layer entirely. The
  symbol layer is now unconditional in both places and the fill-legend block is
  conditional rather than terminal.
- The classes most affected are exactly the ones that cannot produce an LOH
  tract by construction — `CO_GC_subres`, `NCO_GC_subres`,
  `CROSSOVER_NO_TRACT`, `GC_ONE_SIDED` — so the plot was hiding events
  precisely where the tract was too short for the LOH map to resolve.
- Sub-resolution calls are no longer automatically "review". `CO_GC_subres` and
  `NCO_GC_subres` are promoted from a self-classifying peak
  (`gene_conversion` / `internal_crossover`) whose per-read haplotype pattern
  is the authoritative evidence — which is why `classify_tract()` deliberately
  exempts these edge types from its `min_span` gate. Scoring them "review"
  regardless of depth meant a crossover with 60 supporting reads was flagged
  the same as one with 2. `build_event_table()` now applies the `min_span`
  floor itself: at or above it the call is "high", below it "review". The
  `_subres` suffix is unchanged and still records that no fixed tract was
  resolvable. Every other class keeps its previous confidence.
- Effect on results: on the synthetic crossover sets (`test_data/co_1bp.csv.gz`
  and `co_100bp.csv.gz`, 20 known crossovers each) all 20 events now plot and
  score high, where `co_1bp` previously plotted none and scored 0 high / 20
  review. `co_100bp` had been labelling correctly only by accident — a single
  junction resolved a 2 bp `REF_fixed` island, which flipped the LOH gate on
  and let all 20 symbols through. Regression-checked on `test_data` chrIII: 4
  events, LOH band, legend, caption and symbols unchanged, with
  `TERMINAL_DELETION` correctly still "review".

## [0.8.10] - 2026-08-03

- Per-read zone calls now require read evidence. `classify_zone_state()`
  resolved a zone by majority vote with no floor on the evidence behind it, and
  because the test was `majority >= 0.5` a dead tie (a read contributing 1 REF
  and 1 ALT SNP) always resolved to REF — a directional artefact that
  manufactures spurious REF calls and never spurious ALT ones, so it does not
  cancel out with depth. A zone now returns `NA` unless the read contributes at
  least `min_evidence_snps` SNPs to it *and* the vote clears `min_margin`,
  letting the existing `classifiable` filters drop the read as they already do
  when it has no SNPs in the zone at all. The new thresholds live in
  `ZONE_CALL_HEURISTICS` (3 SNPs, 0.15 margin — which keeps 3-SNP 2:1 splits at
  0.167 while rejecting ties).
- Edge calls are now weighted by read count. `classify_edge_type()` tested its
  motifs against `unique(patterns)`, the *set* of distinct L-M-R patterns, so
  one stray read vetoed an otherwise unanimous verdict. Patterns are now
  filtered to those held by at least `(1 - homog_frac)` of classifiable reads,
  and the survivors must account for `homog_frac` of the population before any
  positive call is made; otherwise the edge stays "ambiguous".
- `homog_frac` (0.80) moves into `FUSION_HEURISTICS` as the single source of
  truth and is threaded through the `compute_peak_pairs()` call sites. It was
  previously replicated as a literal function-argument default in three places.
  A `NULL` or out-of-range value now errors instead of silently collapsing
  every comparison to zero length and calling everything "ambiguous".
- Effect on results: clean gene conversions previously reported as
  `AMBIGUOUS(pair_edge:ambiguous)` now classify as `NCO_GC`. On the simulated
  swap sets every such AMBIGUOUS call resolves, with no events gained or lost
  and no unmatched calls; results are identical across `--min-run 1` and `2`,
  removing a sensitivity that previously changed the verdict on three tracts.
  Regression-checked on `test_data` chrXI and chrI (byte-identical) and RAD5_01
  (48 events, unchanged, its one ambiguous pair edge resolving to `NCO_GC`);
  crossover classes (`CO_GC`, `CO_TERM`, `TCO_CAPTURED_TCO`) are preserved.
  Remaining misses are tracts with ≤ 1 SNP, a resolution floor rather than a
  classification failure.

## [0.8.9] - 2026-08-03

- `chimera_cli.R` output flags are now combinable instead of mutually
  exclusive. `--peak-list`, `--overview-rds`, `--events-table` and the new
  `--peak-plot` each request one file, and any combination may be given in a
  single run (e.g. `--peak-list --peak-plot`). Previously choosing more than
  one was an error, so producing several outputs meant re-running the whole
  analysis per file.
- The default (no output flag) is now the three-file set `--peak-plot
  --peak-list --events-table`, where it was the overview PNG alone. Because the
  events table comes from the chain-based caller, a default run now also runs
  chain analysis and the PNG is the annotated overview (LOH band + event
  symbols). Pass `--peak-plot` alone for the previous peaks-only PNG behaviour.
- `--events-table` and `--chain-all` are no longer mutually exclusive. Given
  both, the events table is written to the `--output` path *and* as the
  `_step4_final_events.csv` step file, so each flag still delivers its
  documented output.
- `-o/--output` handling follows the multi-file outputs: a directory (or no
  `--output`) auto-names each file with its own suffix, a file path is used
  verbatim when exactly one output is selected, and with several outputs it is
  treated as a name stem.

## [0.8.8] - 2026-07-31

- Peak-to-SNP mapping is now deterministic. When several SNPs inside a peak
  interval tied on raw chimeric-read count *and* on distance to the smoothed
  peak position, `run_chimera_analysis()` broke the tie with an unseeded
  `sample()`, so repeat runs on identical input could assign different SNPs to
  the same peak — changing `snp_pos`/`snp_n` in the peak table and, through
  peak association, the downstream chain analysis. The remaining tie is now
  resolved by taking the leftmost candidate (`which.min(pos)`), chosen over
  first-row order so the result does not depend on the row order of the
  upstream coverage table. Candidates reaching this branch are equal in height
  and equidistant from the peak, so the choice among them is arbitrary; only
  its stability across runs changes.

## [0.8.7] - 2026-07-28

- `chimera_cli.R --peak-list` now writes only peaks that mapped to a qualifying
  SNP (raw chimeric-read count at or above `--min-peak-height`), matching the
  rows shown in the app's "Peak Summary" table. Previously the CLI exported
  every detected peak, including those whose interval held no SNP above the
  cutoff (`snp_pos = NA`, shown as "None above cutoff" in the app). The run
  summary gains a `Peaks above min height` line so the reported count and the
  CSV row count agree. Chain/fusion analysis is unaffected — it continues to
  receive the unfiltered peak table and applies its own filtering internally.

## [0.8.6] - 2026-07-13

- Single-chromosome coverage views now offer image and data downloads. Each
  per-chromosome tab in the "Chromosome plots" panel gains a **Download Plot
  Image (.png)** link and a **Download Plot Data (.rds)** link. The PNG renders
  the exact on-screen figure (shared `build_chr_cov_plot()` helper), and the
  `.rds` bundles the SNP coverage, fitted curve, peaks, and SNP peaks for the
  chromosome plus its LOH strip data (fixed-haplotype `REF_fixed`/`ALT_fixed`
  segments) and event labels when available, with `strain_ref`/`strain_alt`
  labels and metadata.
- Selected-region read plots now include their LOH strip data in the downloaded
  `.rds`. The LOH band that was drawn on the plot is persisted (`loh_data`)
  alongside `strain_ref`/`strain_alt` so the download reflects everything shown.
- Provisional terminal-crossover calling for peakless, gap-masked terminal LOH
  (`CO_TERM_PROBABLE`, rendered `TCO*`, `review` confidence). New rule
  `rule_terminal_loh_gapped_nopeak` (R02d) fires on `H {G:wide} [F]{tel}`: a
  terminal LOH tract reaching the telomere at diploid (2N) depth with no
  chimeric peak on its proximal boundary, where a wide SNP-desert gap (default
  `>= provisional_tco_min_gap_bp = 10000` bp) separates the bounding HET from
  the LOH. A recombination junction inside a gap that wide leaves no read able
  to span it, so no binary peak can form even though the tract is a genuine
  terminal crossover. The wide gap is the positive evidence for the missing
  peak; without it (HET abutting F) the tract stays `UNCATEGORIZED_LOH`, as the
  disabled `rule_terminal_no_peak` (R02b) intended.
- Fused gapped terminal crossovers are now surfaced and enumerated. When R02g
  (`rule_terminal_loh_gapped`) walks a binary-peaked proximal LOH across a wide
  gap to a telomere-reaching tract and the haplotype *switches* across the gap
  (e.g. `REF_fixed -> ALT_fixed`), the event is reported as the new class
  `CO_TERM_GAPPED` (`high` confidence) instead of a single `CO_TERM`. The number
  of fused terminal crossovers (proximal peak-confirmed junction + one per
  internal gap-switch) is carried in a new `n_fused` events-table column and
  rendered on the overview map as `N X TCO*`; the notes enumerate each
  transition. A wide gap that merely splits one *same-state* tract (no switch)
  stays a plain `CO_TERM`.
- Overview-map symbols: removed all `?` marks. `CO_TERM_PROBABLE` is now `TCO*`
  and `AMBIGUOUS(low_coverage)` is now `*`; the asterisk consistently marks a
  gap-inferred junction across the terminal-crossover family (`TCO*`,
  `N X TCO*`). To keep the two gene-conversion review classes distinct without
  a `?`, `GC_ONE_SIDED` is now a bold omicron with a superscript one and
  `GC_UNRESOLVED` a bold omicron with a degree sign.
- Phase-based rescue of `undefined` peaks. When a peak's consensus run pattern
  is a HET-bounded fixed-allele island (e.g. the `REF-HET-ALT-HET` signature a
  broad smoothed peak produces when it straddles an LOH island plus a
  neighbouring tract), `classify_peak_haplotype()` now phases every read
  spanning the island by its left- and right-flank majority allele instead of
  relying on population consensus — which reads HET on both flanks because the
  two reciprocal crossover orientations average out. A switch fraction
  `>= 0.80` is relabeled `internal_crossover` (crossover), `<= 0.20`
  `gene_conversion` (NCO); otherwise the peak stays `undefined`. This uses all
  reads (`full_read_loh`), so non-switching reads stay in the denominator —
  the ratio that separates a crossover from a fixed conversion patch on a het
  background. The island is located in the full-read consensus (not the
  chimeric subset), so only a genuine population-level fixed LOH tract — the
  same signal the LOH HMM sees — triggers a rescue; a fixed tract merely
  abutting a HET region is not misread as a het-bounded island. A rescued
  event's reported footprint is recalculated to that excised island tract
  rather than the broad smoothed peak window or a neighbouring LOH tract the
  motif scanner happened to anchor to. New `phase_call` / `phase_switch_frac`
  columns are surfaced in the peak table and the switch fraction is carried
  into the Recombination Events table. `full_read_loh` is threaded through
  `label_snp_peaks_haplotypes()` / `compute_peak_pairs()`; when absent the
  prior `undefined` behaviour is unchanged.
- Added an interstitial hemizygous-deletion call (`DELETION`, rendered as a
  `Δ` on the overview map). A HET-bounded fixed-allele LOH tract whose
  SNP-site read depth drops below `depth_drop` (default 0.60) of its HET
  flanks is labeled a deletion of the missing homolog rather than a
  copy-neutral LOH tract. Uses the existing SNP-position coverage map only
  (no BAM/CNV/breakpoint modelling), so it is reported at `review`
  confidence. The *lower* of the two flanks is used as the diploid reference,
  so the ratio only falls below threshold when the tract is depth-depressed
  against BOTH flanks — a single anomalously HIGH flank (common for short
  sub-telomeric HET islands, whose repetitive content inflates mapped depth)
  no longer manufactures a false deletion out of a full-depth tract. A
  resolved haplotype-switch junction peak (binary / crossover / gene
  conversion, with real spanning support) also vetoes the call: such a peak
  requires both homologs across the junction, which is incompatible with a
  hemizygous deletion, so the tract is deferred to the recombination rules
  (crossover / gene-conversion). Complements the existing terminal-deletion
  rule (R01).
- `rule_tel_adjacent_het_loh` (R02c, `CO_TERM_PROBABLE`) now recognises a
  subtelomeric MISALIGNMENT pileup interposed between an LOH tract and the
  telomere. Previously the interposed HET had to be within `tel_tol_bp` (5 kb)
  of the chromosome end for the tract to count as terminal; a wider island
  deferred the call to `rule_one_sided_binary` (`GC_ONE_SIDED`) — or, if the
  island's inflated depth dragged the tract's flank ratio below threshold, to a
  false `DELETION`. Yeast subtelomeres are repetitive, so reads from paralogous
  ends pile up there and manufacture a spurious HET island with anomalously
  HIGH depth. Such an island (depth ratio `>= subtel_misalign_depth`, default
  1.4, AND narrower than `subtel_misalign_max_bp`, default 25 kb — the width
  gate keeps a genuine wide diploid arm with globally elevated depth from being
  swallowed) is now treated as artifactual: the LOH tract behind it is called a
  probable terminal crossover (`CO_TERM_PROBABLE`, review confidence) rather
  than a one-sided gene conversion or a deletion. Thresholds calibrated on
  TEL-adjacent HET depth ratios across the test panel (real arms 0.77–1.26;
  thin subtelomeric pileups 1.4–1.9).

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
