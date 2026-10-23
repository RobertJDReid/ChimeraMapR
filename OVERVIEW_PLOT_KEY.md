# Overview plot key

Symbols and markings drawn on the genome-wide overview plot (`build_overview_plot()` in `chimera_functions.R`).

### Recombination event symbols

Event symbols are drawn in bold just below the LOH band, centred on the bp midpoint of each event's span. They come from `EVENT_SYMBOL_MAP` in `chimera_functions.R`; event classes not in the map are not drawn.

| Symbol | Event class(es) | Meaning |
|---|---|---|
| ✖ | `CO_GC` | Crossover with an associated gene conversion tract |
| ✖ | `CO_GC_subres` | Crossover with a conversion tract too short for the LOH map (some flanking SNPs show a fixed call) |
| ✖ | `CROSSOVER_NO_TRACT` | Crossover with no detectable tract (all flanking SNPs HET); called from peak evidence alone |
| Ο | `NCO_GC` | Non-crossover gene conversion |
| Ο | `NCO_GC_subres` | Non-crossover gene conversion with no flanking LOH tract (below LOH resolution) |
| Ο | `NCO_GC_in_terminal` | Non-crossover gene conversion sitting inside a terminal LOH tract |
| Ο | `NCO_GC_in_tract` | Small opposite-state fragment inside an interstitial LOH tract (patchy conversion, or a small LOH overlaid by a larger one) |
| Ο° | `GC_UNRESOLVED` | Confirmed gene-conversion tract (LOH flanked by two binary peaks), but NCO vs CO can't be determined |
| Ο¹ | `GC_ONE_SIDED` | One confirmed binary junction peak; the other boundary has no peak at all |
| ΟΟ | `DOUBLE_GC` | Two NCO gene conversions resolved within one token |
| TCO | `CO_TERM` | Terminal crossover (LOH running to the telomere) |
| TCO* | `CO_TERM_PROBABLE` | Provisional terminal crossover; its only junction is gap-masked (no peak) |
| N X TCO* | `CO_TERM_GAPPED` | N terminal crossovers fused across a SNP gap; the * marks the junction(s) inferred from the gap (falls back to plain "TCO*" if no count is available) |
| 2 X TCO | `TCO_CAPTURED_TCO` | A terminal crossover laid over an earlier terminal crossover |
| Δ | `TERMINAL_DELETION` | Terminal deletion |
| Δ | `DELETION` | Interstitial hemizygous deletion (HET-bounded, about half depth) |
| * | `AMBIGUOUS(low_coverage)` | Too few spanning reads to make a call (fewer than `min_span` = 3); shown as a marker for manual review |

### Other markings

| Mark | Meaning |
|---|---|
| Small black open circles | Read coverage at each SNP (y = number of reads) |
| Large pale-green filled circles with a black ring | SNPs inside a detected chimeric-read peak |
| Grey line | Uniform-coverage fit for each chromosome |
| Blue band under the x-axis | REF-fixed LOH (labelled with the REF strain name) |
| Red band under the x-axis | ALT-fixed LOH (labelled with the ALT strain name) |
| Pale yellow panel background, black border | Chromosome estimated at less than 2N ("<2N" in the legend) |
| Pale lilac panel background, black border | Chromosome estimated at more than 2N (">2N" in the legend) |
| White panel background | Diploid (2N) chromosome |
| Thin red vertical gridlines | Major x breaks |
| Thin black vertical gridlines | Minor x breaks, every 100 kb |
