# Mitochondrial DNA Analysis Pipeline

A reproducible end-to-end pipeline for mitochondrial DNA variant calling, heteroplasmy calculation, haplogroup classification, and functional annotation from whole-genome short read sequencing data (tested on illumina data).
This pipeline takes BAM files as input, extract mitochondrial reads and realign them again into the reference genome of interest (the revised cambridge reference sequence is used here).

---

## Overview

```
 BAM files
     │
     ▼
[1] Extract chrM BAM                (samtools)
     │
     ▼
[2] BAM → FASTQ                     (samtools)
     │
     ▼
[3] Align to rCRS + Pileup          (bwa mem + samtools mpileup + pileup_analysis.py)
     │                               ↳ heteroplasmy calculated here directly
     ▼
[4] Quality Filtering                (pileup_filter.R)
     │
     ▼
[5] Haplogroup Classification       (csv_to_vcf.py + Haplogrep3)
     │
     ▼
[6] Remove Haplogroup-Defining      (filter_haplo_variants.py)
    Variants
     │
     ▼
[7] Functional Annotation           (apply_annotation.py)
     │
     ▼
merged_annotated.csv
```

---

## Repository Structure

```
mtDNA_pipeline/
├── run_pipeline.sh                  # Main pipeline script (7 steps)
├── README.md                        # This file
├── environment.yml                  # Conda environment (Python + R)
├── scripts/
│   ├── pileup_analysis.py           # Parse samtools mpileup + calculate heteroplasmy → CSV output
│   ├── pileup_filter.R              # Quality filtering (you can change the parameters: ((Heteroplasmy > 0.05 & AF > 0.05) | AF > 0.08) &  (unmapped / Depth < 0.7)) + plots
│   ├── csv_to_vcf.py                # Convert filtered CSV → VCF for Haplogrep
│   ├── filter_haplo_variants.py     # Remove haplogroup-defining variants 
│   └── apply_annotation.py         # Merge annotation databases
├── reference/
│   └── rCRS.fasta                   # Place reference genome here (see below)
└── annotation_databases/
├── curated/ # examples of .csv annotation files 
│ ├── gnomAD_filtred.csv
│ ├── HelixMTdb_modified.csv
│ ├── MitImpact-db-3.1.3_modified.csv
│ ├── MitoMAP_RNA_tRNA.csv
│ ├── t_APOGEE_data_modified.csv
│ └── Filtre_data_bases.R
├── raw/ # if you want to use other databases, place them here
└── README.md # Description of annotation sources
```

---

## Dependencies

### System tools

| Tool | Version tested | Purpose |
|------|----------------|---------|
| [samtools](http://www.htslib.org/) | ≥ 1.17 | BAM manipulation, pileup |
| [BWA](https://github.com/lh3/bwa) | ≥ 0.7.17 | Read alignment |
| [Haplogrep3](https://github.com/genepi/haplogrep3) | ≥ 3.2 | Haplogroup classification |
| Python | ≥ 3.9 | Pileup parsing, annotation |
| R | ≥ 4.2 | Variant filtering, visualisation |

### Installation

#### Option A — Conda (recommended)

```bash
# Clone the repository
git clone https://github.com/WafaGhoul196/short_read_mtDNA_anlysis.git
cd short_read_mtDNA_anlysis

# Create and activate the environment
bash setup.sh
conda activate mtdna_pipeline
```

#### Option B — Manually
if problems with downloading haplogrep do it manually
**Haplogrep3:**
```bash
# Download the latest release
wget https://github.com/genepi/haplogrep3/releases/latest/download/haplogrep3.zip
unzip haplogrep3.zip -d haplogrep3
sudo mv haplogrep3/haplogrep3 /usr/local/bin/
haplogrep3 --version   # verify
```

**Python dependencies:**
```bash
pip install pandas
```

**R dependencies:**
```r
install.packages(c("dplyr", "tidyr", "ggplot2"))
```

---

## Reference Genome

The pipeline aligns against the **revised Cambridge Reference Sequence (rCRS)** of the human mitochondrial genome (NC_012920.1).

```bash
# Download from:  https://www.mitomap.org/MITOMAP/HumanMitoSeq
```

---

## Sample CSV Format

This file contains the sample IDs
The input CSV must have **Sequencing_number** in **column 6** (1-based), with a header row.

```
ID,Name,Date,Project,Status,Sequencing_number,...

```

---

## Annotation Databases

Here are the databases used for annotation:

| File | Source |
|------|--------|
| `gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv` | [gnomAD v3.1 mtDNA](https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna) |
| `HelixMTdb_20200327.tsv` | [HelixMTdb](https://helix.com/pages/mitochondrial-variant-database) |
| `MitImpact-db-3.1.3.txt` | [MitImpact](https://mitimpact.mcb2lab.org/) |
| `MutationsRNA MITOMAP Foswiki.csv` | [MITOMAP](https://www.mitomap.org/MITOMAP) |
| `t_APOGEE_2024.0.1.txt` | [t-APOGEE](https://mitimpact.mcb2lab.org/) |

These databases have been homogenized using this code: Filter_data_bases.R (this code is found in the `raw_databases/` folder)

You can annotate with other databases of interest. Just update the annotation code.
Each database file should contain at least a **Pos** column, and where relevant **Ref** and **Alt** columns for allele-specific annotation.

---

## Usage

### Full pipeline (all steps)

```bash
bash run_pipeline.sh \
  -c data_filtered.csv \
  -b /path/to/DRAGEN_results \
  -o /path/to/output \
  -r reference/rCRS.fasta \
  -t 16
```

### Run specific steps only

Steps are: `extract`, `fastq`, `pileup`, `filter`, `haplogrep`, `haplo_filter`, `annotate`

```bash
# Re-run only haplogroup filtering and annotation
bash run_pipeline.sh -c samples.csv -o results/ -s haplo_filter,annotate

# Re-run from haplogrep onward
bash run_pipeline.sh -c samples.csv -o results/ -s haplogrep,haplo_filter,annotate
```

### All options

```
-c  Path to sample CSV file        (required)
-b  DRAGEN results base directory  [default: /mnt/DRAGEN_pipeline_results/wgs]
-o  Output directory               [default: ./mtdna_results]
-r  Reference genome FASTA (rCRS)  [default: ./reference/rCRS.fasta]
-H  Haplogrep3 binary path         [default: haplogrep3]
-t  CPU threads                    [default: 8]
-s  Steps to run (comma-separated or 'all')
-h  Show help
```

---

## Output Files

```
<output_dir>/
├── pipeline.log                          # Full run log
├── <SAMPLE>_chrM.bam                     # Extracted mitochondrial BAM
├── <SAMPLE>_chrM.bam.bai
├── fastq/
│   ├── <SAMPLE>_chrM_R1.fq.gz
│   └── <SAMPLE>_chrM_R2.fq.gz
├── results/
│   ├── mitochondrial_<SAMPLE>.bam        # Re-aligned to rCRS
│   ├── mitochondrial_pileup_<SAMPLE>.txt # Raw pileup
│   ├── output_pileup_analysis_<SAMPLE>.csv               # Parsed pileup
│   ├── output_pileup_analysis_<SAMPLE>_heteroplasmy.csv  # Corrected heteroplasmy
│   ├── <SAMPLE>_data_filtred_heteroplasmy.csv            # Filtered variants
│   ├── all_patients_filtered_data.csv                    # All samples merged
│   └── plots/
│       ├── filter_stats_barplot.png
│       └── heteroplasmy_heatmap.png
├── haplogrep/
│   ├── <SAMPLE>_formatted.vcf
│   └── merged_haplogrep_results.txt      # Final haplogroup table
├── haplogroup_filtered/
│   ├── <SAMPLE>_no_haplo_variants.csv    # Variants without haplogroup-defining ones
│   └── <SAMPLE>_haplo_variants.csv       # Haplogroup-defining variants (removed)
└── annotation/
    └── merged_annotated.csv              # Final annotated variant table
```

---

## Heteroplasmy Formula

The formula is implimanted in `scripts/recalculate_heteroplasmy.py`:

```
Heteroplasmy = (A + T + G + C + ins + del − RefBase_count)
               ─────────────────────────────────────────────
                      A + T + G + C + ins + del
```

Note: only the first three major insertions and deletions are taken in account.

Thank you, Wafa :) 
