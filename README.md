![ChimeraMapR](images/logo.png)

_version 0.1.0_

An interactive __R Shiny app__ for identifying **chimeric long reads** by tracking allele changes across SNP positions, then summarizing where chimeric reads cluster along each chromosome using LOESS smoothing and peak detection.

The app is designed for plotting per-read base calls at known SNP sites and is the
plotting function for the [bamCol](https://github.com/RobertJDReid/bamCol) python app.

---

## Features

- Upload three input files: **read-level SNP calls**, **SNP definition table**, and **chromosome sizes (FAI)**
- Detects chimeric reads where allele identify switches
- Creates an **overview per-chromosome plot** (counts + LOESS fit + peak markers)
- Reports a **peak summary table** (positions, boundaries, heights)
- Generates **per-peak read-level plots** showing all reads intersecting the peak SNP
- Exports:
  - overview plot (PNG)
  - peaks table (CSV)
  - chimeric read IDs (TXT)

---

Exported read IDs can be used to filter the original alignment file giving a bam file
of just the chimeric reads for use as a genome browser track.

## Installation

### Requirements
- R (recommended: ≥ 4.1)
- R packages:
  - shiny
  - tidyverse
  - pracma

Install packages in R:

```r
install.packages(c("shiny", "tidyverse", "pracma"))
```

---

## Running the app locally (reccomended)

```bash
R -e "shiny::runApp('path/to/app.R',launch.browser = TRUE)"
```

or from within R:

```r
shiny::runApp("app.R")
```

---

## Input files

### 1) Read Data File (CSV / GZ)

Output from `bamCol.py` is used for read data at SNP positions.
gzipped files can be read without unzipping.

Required columns:
- `chrom`
- `pos`
- `read_id`
- `call`
- `is_del`

Rows with `is_del == 1` are filtered out.

Typical origin: per-read SNP pileups generated from minimap2 + bcftools or pysam.

### 2) SNP Data File (CSV)

Required columns:
- `CHROM`
- `POS`
- `REF`
- `ALT`

### 3) Chromosome Size File (FAI)

Uses the first two columns:
- chromosome name
- chromosome length

---

## Example input files

The following minimal examples illustrate the expected format of each input file.
These are intentionally small and suitable for testing that the app runs end-to-end.

### Example 1: Read-level SNP calls (`example_reads.csv`)

```csv
chrom,pos,read_id,call,is_del
S288C_chrI,1001,readA,A,0
S288C_chrI,1005,readA,G,0
S288C_chrI,1001,readB,A,0
S288C_chrI,1005,readB,A,0
S288C_chrI,1010,readB,G,0
S288C_chrI,1001,readC,G,0
S288C_chrI,1005,readC,G,0
S288C_chrI,1010,readC,G,0
```

In this example:
- `readA` switches from REF→ALT across SNPs (chimeric)
- `readB` switches from REF→ALT (chimeric)
- `readC` is consistently ALT (not chimeric)

### Example 2: SNP definition table (`example_snps.csv`)

```csv
CHROM,POS,REF,ALT
S288C_chrI,1001,A,G
S288C_chrI,1005,A,G
S288C_chrI,1010,A,G
```

This table defines which bases are considered REF vs ALT at each SNP position.

### Example 3: Chromosome sizes (`example_genome.fai`)

```tsv
S288C_chrI	230218
```

Only the first two columns are used by the app; additional `.fai` columns are ignored.

---

## Analysis overview

- Allele classification per SNP (REF / ALT / OTHER)
- Chimeric read detection using allele switches across SNP positions
- SNP-wise chimeric read counts
- LOESS smoothing per chromosome
- Peak detection using `pracma::findpeaks`
- Per-peak read visualization showing all reads spanning the peak SNP

---

## Outputs

- Overview plot (PNG)
- Peak table (CSV)
- Chimeric read IDs (TXT)

---

## License

MIT License © 2026 Robert J. D. Reid

Contributions welcome.