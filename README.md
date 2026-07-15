![ChimeraMapR](images/logo.png)

_version 0.8.6_

An interactive **R Shiny app** for identifying **haplotype switches** in long-read DNA sequencing data, by tracking allele changes across known SNP positions, summarizing where chimeric reads cluster along each chromosome using Whittaker smoothing and peak detection, and classifying the underlying recombination events (crossovers, terminal crossovers, gene conversions) from the resulting LOH structure.

Chimeric reads contain sequence derived from more than one parental chromosome — for example, reads spanning a recombination breakpoint in a hybrid or cross. The app classifies each base call at known biallelic SNP positions as REF or ALT, then uses run-length encoding to detect reads where the allele identity switches and is sustained across multiple consecutive SNPs, marking a transition between parental haplotypes within a single read. A separate chain-analysis stage calls loss-of-heterozygosity (LOH) tracts along each chromosome and combines them with the detected peaks to label each event as a crossover, terminal crossover, gene conversion, or an ambiguous/review case.

The core input is a CSV of per-read base calls at biallelic SNP positions extracted from a BAM alignment file (read ID, base call, mapping quality, strand, and related fields — see [Features](#features) for the full set of required inputs). These bam file features are extracted using the [bamCol](https://github.com/RobertJDReid/bamCol) python app. _ChimeraMapR_ groups these calls by read, flags all reads containing haplotype switches, and plots per-read base calls at known SNP sites. A batch-mode [command-line interface](#command-line-interface) (`chimera_cli.R`) is also available for running the same analysis outside the Shiny app.

---

## Features

- Upload three input files: **read-level SNP calls (CSV)**, **SNP definition file (CSV or VCF)**, and **chromosome sizes (FAI)**
- Accepts SNP definitions as a simple CSV, a plain VCF, or a gzipped VCF (`.vcf.gz`)
- Detects chimeric reads where allele identity switches and is sustained across consecutive SNPs
- Creates an **overview per-chromosome plot** (counts + curve fit + peak markers), with aneuploid chromosomes shaded by estimated ploidy
- Reports a **peak summary table** (positions, boundaries, heights, and haplotype classification)
- Classifies each peak's underlying read pattern (`binary`, `compound_binary`, `gene_conversion`, `internal_crossover`, `undefined`) from the run of REF/ALT calls spanning it
- **Peak fusion**: evaluates adjacent peak pairs by shared-read Jaccard similarity and haplotype signature, auto-fusing, flagging for supervised review, or keeping independent as appropriate
- **LOH region calling**: a beta-binomial EM + Viterbi HMM (compiled C++ via Rcpp) segments each chromosome into REF-fixed / ALT-fixed / heterozygous tracts
- **Recombination event classification**: a rule-based chain analysis combines LOH tracts and peak calls to label each event (crossover, terminal crossover, gene conversion, one-sided/ambiguous, etc.) with a confidence level, surfaced in the **Recombination Events** tab
- Generates **per-peak read-level plots** showing all reads intersecting the peak SNP, plus fused-peak plots for post-fusion groups
- Chromosome order follows the input FAI file, preserving the original FASTA sequence order
- __Exports:__
  - Overview plot (PNG/RDS)
  - Peaks table, fused-peak table, peak-pair analysis, and post-fusion summary (CSV)
  - LOH regions table (CSV)
  - Recombination event table (CSV) and chain objects (RDS)
  - Whittaker curve fits (CSV)
  - Chimeric read IDs (TXT)

Exported read IDs can be used to filter the original alignment file, producing a BAM file of just the chimeric reads for use as a genome browser track.

---

## Installation

### Requirements
- R (recommended: ≥ 4.1)
- A C++ compiler toolchain (Xcode Command Line Tools on macOS, Rtools on Windows, or `build-essential` on Linux) — required to compile the LOH-calling HMM (`src/loh_hmm.cpp`) via Rcpp on first load
- R packages:
  - shiny ≥ 1.12.0
  - data.table ≥ 1.15.0
  - ggplot2 ≥ 3.4.0
  - pracma ≥ 2.4.4
  - R.utils _for some installs_
  - patchwork
  - Matrix
  - Rcpp
  - igraph
  - optparse (only needed to run `chimera_cli.R`)

_note: shiny must be ≥ 1.12.0 to avoid a region selection bug for individual chromosome plots_

Install packages in R:

```r
install.packages(c("shiny", "data.table", "ggplot2", "pracma", "Matrix", "Rcpp", "igraph", "optparse"))
```

---

## Running the app locally (recommended)

```bash
R -e "shiny::runApp('path/to/app.R', launch.browser = TRUE)"
```

or from within R:

```r
shiny::runApp("app.R")
```

On macOS/Linux, `chimera.sh` can be symlinked into a directory on your `PATH` (e.g. `/usr/local/bin/chimera`) to launch the app with a bare `chimera` command from any directory.

---

## Command-line interface

For batch processing without the Shiny UI, `chimera_cli.R` runs the same analysis pipeline from the terminal:

```bash
Rscript chimera_cli.R [options] <read_data> <snp_data> <chr_size.fai>
```

Three positional arguments are required, in order: the **read data** file (`.csv`/`.csv.gz`), the **SNP data** file (`.csv`/`.vcf`/`.vcf.gz`), and the **chromosome size** file (`.fai`).

#### Output-mode flags

Pick at most one; the default (no flag) is the PNG overview plot. These are mutually exclusive.

| Flag | Default | Output |
|---|:---:|---|
| _(none)_ | ✓ | Genome-wide overview plot (PNG) |
| `--peak-list` | | Detected peaks, including haplotype classification, as CSV |
| `--overview-rds` | | Overview plot object as RDS (re-plot later with `readRDS()` + `print()`) |
| `--events-table` | | Runs the chain-based caller and writes **only** the final events table as CSV (no intermediate step CSVs, no plot). Mutually exclusive with `--chain-all` |

#### Additional output flags

Layer extra outputs on top of the selected output mode.

| Flag | Description |
|---|---|
| `--coverage-map` | Also write the per-position coverage table and collapsed coverage segments as CSVs |
| `--chain-all` | Run the chain-based LOH recombination-event caller and write one CSV per pass (steps 0–4); in the default PNG mode this also annotates the overview plot with the LOH band and event symbols. Requires `loh_chain_analysis.R` alongside `chimera_functions.R` |

#### Analysis parameters

| Flag | Default | Description |
|---|:---:|---|
| `-n, --sample-name NAME` | `Sample_01` | Sample name used in output filenames and plot titles |
| `--mapq-cutoff INT` | `20` | Minimum MAPQ value; reads below this are excluded |
| `--baseq-cutoff INT` | `10` | Minimum base quality at the SNP position; calls below this are excluded |
| `--del-rate-cutoff FLOAT` | `0.10` | SNPs where more than this fraction of confidently-mapped reads register a deletion are excluded |
| `--min-run INT` | `2` | Minimum consecutive same-allele calls to count as a run |
| `--min-peak-height INT` | `10` | Minimum transition count at a SNP to qualify as a peak |
| `-l, --lambda FLOAT` | `1` | Whittaker smoothing penalty λ (lower = tighter fit) |

#### Chain-analysis tuning parameters

Apply only when `--chain-all` or `--events-table` is used.

| Flag | Default | Description |
|---|:---:|---|
| `--tel-tol KB` | `5` | Telomere tolerance in kb |
| `--merge-gap KB` | `5` | Same-state NA-gap merge threshold in kb |
| `--min-span INT` | `3` | Minimum spanning reads required for a read-based chain call |
| `--peak-pad BP` | `200` | Peak association padding in bp; a peak is attached to a token junction if its SNP position falls within this distance of the token boundary |

#### Output path

| Flag | Default | Description |
|---|:---:|---|
| `-o, --output PATH` | _cwd_ | Output file path or directory. If a directory (or omitted), a dated filename is auto-generated; the extension (`.png`/`.csv`/`.rds`) follows the output mode |

Run `Rscript chimera_cli.R --help` for the same list from the terminal.

---

## Input files

### 1) Read Data File (CSV / GZ)

Output from `bamCol.py` is used for read data at SNP positions. Gzipped files can be read without unzipping.

Required columns:

| Column | Description |
|---|---|
| `chrom` | Chromosome name |
| `pos` | Position on chromosome |
| `read_id` | Unique read identifier |
| `call` | Base call at this position |
| `mapq` | Read mapping quality score |
| `base_qual` | Base quality score at this position |
| `is_del` | Deletion flag (rows with `is_del == 1` are filtered out) |

Typical origin: per-read SNP pileups generated from minimap2 + bcftools or pysam.

---

### 2) SNP Definition File (CSV or VCF)

Defines which base is REF and which is ALT at each SNP position. Two formats are accepted:

#### Option A: CSV

A plain CSV with the following required columns:

| Column | Description |
|---|---|
| `CHROM` | Chromosome name |
| `POS` | Position on chromosome |
| `REF` | Reference allele |
| `ALT` | Alternate allele |

#### Option B: VCF (plain or gzipped)

Standard VCF format (v4.x) is accepted, including gzipped `.vcf.gz` files. The app automatically detects the format from the file extension. Only biallelic SNP sites are used — indels, multi-allelic sites, and any site where REF or ALT is not a single nucleotide are silently skipped.

VCF files can be generated directly from bcftools:

```bash
bcftools mpileup -f reference.fasta sample.bam | bcftools call -mv -Oz -o variants.vcf.gz
```

---

### 3) Chromosome Size File (FAI)

Standard FASTA index format generated by `samtools faidx`. Only the first two columns are used:

| Column | Description |
|---|---|
| column 1 | Chromosome name |
| column 2 | Chromosome length (bp) |

Chromosomes are displayed in the order they appear in this file, which matches the order of sequences in the original FASTA.

Generate with:

```bash
samtools faidx reference.fasta
```

---

## Analysis parameters

All parameters are set in the sidebar panel before running the analysis.

| Parameter | Default | Description |
|---|---|---|
| **Sample Name** | _(blank)_ | REF strain label used in exported file names and LOH band legends |
| **ALT Strain Name** | _(blank)_ | Optional ALT strain label for LOH band legends; generic REF/ALT labels are used if blank |
| **Minimum Run Length** | `2` | Minimum number of consecutive same-allele calls required to count as a sustained run. Increase for noisier or lower-coverage data |
| **Minimum Peak Height** | `10` | Minimum chimeric read count for a peak to be reported. A value of approximately half the median read depth is recommended |

### Advanced settings

Collapsed by default in the sidebar; also exposed as `chimera_cli.R` flags.

| Parameter | Default | Description |
|---|---|---|
| **Minimum MAPQ Value** | `20` | Minimum mapping quality score; reads below this threshold are excluded |
| **Base Quality Minimum** | `10` | Minimum base quality score at the SNP position; calls below this are excluded |
| **Max Local Deletion Rate** | `0.10` | SNPs where more than this fraction of confidently-mapped reads register a deletion (rather than a base call) are excluded — typically a sign of a nearby repeat/homopolymer destabilizing alignment |
| **Whittaker Lambda (λ)** | `1` | Smoothness penalty for the Whittaker smoother. Lower = tighter fit (preserves sharp peaks); higher = smoother curve |
| **Jaccard Threshold for Auto-Fusion** | `0.20` | Minimum Jaccard index of peak-associated reads to trigger automatic peak fusion |
| **Telomere Tolerance (Kb)** | `5` | LOH regions within this distance of a chromosome end are treated as terminal |
| **Same-State Merge Gap (Kb)** | `5` | NA gaps smaller than this between two same-state LOH runs are merged |
| **Min Spanning Reads for Read Call** | `3` | Minimum spanning reads required to make a read-based event call; below this, the result is Ambiguous |
| **Peak Association Padding (bp)** | `200` | A peak is attached to a token junction if its SNP position falls within this distance of the token boundary |
| **NCO-GC Homogeneity Fraction** | `0.80` | Fraction of spanning reads that must share the return pattern (e.g. REF-ALT-REF) to call a non-crossover gene conversion (NCO-GC) |

---

## Example input files

The following minimal examples illustrate the expected format of each input file. These are intentionally small and suitable for testing that the app runs end-to-end.

### Example 1: Read-level SNP calls (`example_reads.csv`)

```csv
chrom,pos,read_id,call,mapq,base_qual,is_del
S288C_chrI,1001,readA,A,60,35,0
S288C_chrI,1005,readA,G,60,34,0
S288C_chrI,1001,readB,A,55,30,0
S288C_chrI,1005,readB,A,55,31,0
S288C_chrI,1010,readB,G,55,33,0
S288C_chrI,1001,readC,G,48,29,0
S288C_chrI,1005,readC,G,48,32,0
S288C_chrI,1010,readC,G,48,30,0
```

In this example:
- `readA` switches from REF→ALT across SNPs (chimeric)
- `readB` switches from REF→ALT (chimeric)
- `readC` is consistently ALT (not chimeric)

### Example 2a: SNP definition table as CSV (`example_snps.csv`)

```csv
CHROM,POS,REF,ALT
S288C_chrI,1001,A,G
S288C_chrI,1005,A,G
S288C_chrI,1010,A,G
```

### Example 2b: SNP definition table as VCF (`example_snps.vcf`)

```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
S288C_chrI	1001	.	A	G	.	.	.
S288C_chrI	1005	.	A	G	.	.	.
S288C_chrI	1010	.	A	G	.	.	.
```

Both formats define the same REF/ALT alleles at each SNP position. The app detects the format automatically from the file extension.

### Example 3: Chromosome sizes (`example_genome.fai`)

```tsv
S288C_chrI	230218	6	60	61
```

Only the first two columns are used; additional `.fai` columns are ignored.

---

## Analysis overview

1. **Allele classification** — each base call is classified as REF, ALT, or OTHER by joining to the SNP definition file. Calls are first filtered by MAPQ, base quality, and local deletion-rate cutoffs
2. **Chimeric read detection** — run-length encoding identifies reads with sustained allele switches across consecutive SNP positions
3. **SNP-wise coverage** — chimeric reads are counted at each SNP position across all chromosomes
4. **Whittaker smoothing** — a smoothed curve is fit per chromosome
5. **Peak detection** — `pracma::findpeaks` identifies peaks in the smoothed signal above the user-defined minimum height, each classified by its underlying REF/ALT run pattern (`binary`, `compound_binary`, `gene_conversion`, `internal_crossover`)
6. **Peak fusion** — adjacent peaks are evaluated by shared-read Jaccard similarity and haplotype signature, then auto-fused, flagged for supervised review, or left independent
7. **LOH calling** — a beta-binomial EM + Viterbi HMM segments each chromosome into REF-fixed, ALT-fixed, and heterozygous tracts; chromosome-wide ploidy is estimated to flag aneuploid chromosomes
8. **Recombination event classification** — a rule-based chain analysis combines LOH tracts with (fused) peak calls to label each event as a crossover, terminal crossover, gene conversion, or an ambiguous/review case, with a confidence level
9. **Visualization** — per-peak and fused-peak plots show all chimeric reads spanning each event, coloured by REF (blue) vs ALT (red)

---

## Outputs

| File | Format | Description |
|---|---|---|
| Overview plot | PNG / RDS | Per-chromosome read counts with fitted curve, peak markers, and aneuploidy shading |
| Peak table / fused-peak table | CSV | Peak positions, boundaries, heights, and haplotype classification for all chromosomes |
| Peak-pair analysis / post-fusion summary | CSV | Jaccard scores and fusion decisions; blended pre/post-fusion peak classification |
| LOH regions | CSV | REF-fixed / ALT-fixed / heterozygous tract boundaries per chromosome |
| Recombination event table | CSV | Called events (crossover, terminal crossover, gene conversion, ambiguous, etc.) with confidence level |
| Chain objects | RDS | Intermediate chain-analysis objects for the Recombination Events tab |
| Curve fits | CSV | Fitted Whittaker curves and run parameters per chromosome |
| Chimeric read IDs | TXT | One read ID per line, for use with `samtools view -N` |

---

## License

MIT License © 2026 Robert J. D. Reid

Contributions welcome.
