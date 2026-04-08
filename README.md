# VGP Genome Assembly Pipeline
### High-Quality Chromosome-Level Genome Assembly using PacBio HiFi, Bionano Optical Maps, and Illumina Hi-C Data

---

## Table of Contents
1. [Project Overview](#1-project-overview)
2. [Background](#2-background)
3. [Pipeline Overview](#3-pipeline-overview)
4. [Input Data](#4-input-data)
5. [Tools and Software](#5-tools-and-software)
6. [Pipeline Stages](#6-pipeline-stages)
7. [Results Summary](#7-results-summary)
8. [Repository Structure](#8-repository-structure)
9. [How to Reproduce](#9-how-to-reproduce)
10. [Key Concepts](#10-key-concepts)
11. [References](#11-references)

---

## 1. Project Overview

This project implements the **Vertebrate Genome Project (VGP)** assembly pipeline to generate a high-quality, near-error-free, gap-free, chromosome-level, haplotype-phased genome assembly. The pipeline was executed using **Galaxy**, a web-based bioinformatics platform, and uses a combination of three complementary sequencing technologies to produce the highest possible assembly quality.

**Organism:** *Saccharomyces cerevisiae* S288C (Baker's Yeast)  
**Reason for selection:** One of the most intensively studied eukaryotic model organisms, allowing precise evaluation of assembly quality against a well-known reference genome.  
**Data type:** Synthetic HiFi reads generated from the S288C genome, simulating a diploid organism.  
**Platform:** Galaxy (usegalaxy.org)  
**Assembly trajectory:** HiFi + Hi-C + Bionano (Trajectory D)

---

## 2. Background

### Why Genome Assembly?
Deciphering the complete genetic blueprint of an organism requires reconstructing its genome from millions of short sequencing reads. This process, known as **de novo genome assembly**, is one of the most computationally challenging problems in modern bioinformatics.

### Challenges in Genome Assembly
Two main factors determine the difficulty of genome assembly:

**1. Repetitive Elements**
- **Interspersed repeats:** Transposable elements (TEs) that occur at multiple locations throughout the genome
- **Tandem repeats (TRs):** Repeated sequences occurring at a single locus
- Repetitive elements constitute over one-third of mammalian genomes
- They are the main cause of assembly discontinuity and gaps

**2. Heterozygosity**
- Diploid organisms have two copies of each chromosome (one from each parent)
- Regions where the two copies differ are called **heterozygous regions**
- Assemblers must correctly separate these into two distinct haplotypes
- Failure leads to **haplotype switching errors** and false duplications

### The Vertebrate Genome Project (VGP)
The VGP is an international consortium (Genome 10K) with the goal of generating high-quality reference genome assemblies for every vertebrate species on Earth. The VGP pipeline combines:
- **PacBio HiFi reads** for accurate long-read assembly
- **Bionano optical maps** for long-range scaffolding
- **Hi-C chromatin data** for chromosome-level scaffolding and haplotype phasing

### Why PacBio HiFi Reads?
HiFi (High Fidelity) sequencing, introduced by PacBio in 2020, produces reads:
- **10–25 kbp** in length (much longer than Illumina's ~150 bp)
- **>99% accuracy** (Q20 or higher)
- Capable of spanning repetitive regions that short reads cannot resolve
- No amplification bias (uses native DNA)

---

## 3. Pipeline Overview

The VGP pipeline follows a modular structure with multiple quality control checkpoints:

```
┌─────────────────────────────────────────────────────────────────┐
│                    VGP ASSEMBLY PIPELINE                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  INPUT DATA                                                     │
│  ├── PacBio HiFi reads (3 FASTA files, 50x coverage)           │
│  ├── Illumina Hi-C reads (Forward + Reverse)                   │
│  └── Bionano optical maps (CMAP format)                        │
│                                                                 │
│  STAGE 1: HiFi READ PREPROCESSING                              │
│  └── Cutadapt → Remove adapter-containing reads                │
│                                                                 │
│  STAGE 2: GENOME PROFILING                                      │
│  ├── Meryl → K-mer counting (k=31)                             │
│  └── GenomeScope2 → Estimate genome size, heterozygosity       │
│                                                                 │
│  STAGE 3: CONTIG ASSEMBLY                                       │
│  └── Hifiasm (Hi-C mode) → Hap1 + Hap2 phased contigs         │
│                                                                 │
│  STAGE 4: ASSEMBLY QUALITY CONTROL                              │
│  ├── gfastats → Assembly statistics (N50, contig count, etc.)  │
│  ├── BUSCO → Gene completeness assessment                      │
│  └── Merqury → K-mer based quality value (QV)                 │
│                                                                 │
│  STAGE 5: BIONANO SCAFFOLDING                                   │
│  └── Bionano Hybrid Scaffold → Optical map integration         │
│                                                                 │
│  STAGE 6: HI-C SCAFFOLDING                                      │
│  ├── BWA-MEM2 → Map Hi-C reads to assembly                     │
│  ├── Filter and merge → Process chimeric reads                 │
│  ├── PretextMap → Generate Hi-C contact map (pre-YaHS)        │
│  ├── YaHS → Chromosome-level scaffolding                       │
│  ├── BWA-MEM2 → Remap Hi-C reads to final scaffolds           │
│  └── PretextMap + Pretext Snapshot → Final contact map         │
│                                                                 │
│  OUTPUT                                                         │
│  └── Chromosome-level haplotype-phased genome assembly         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 4. Input Data

| Data Type | Format | Files | Coverage | Source |
|---|---|---|---|---|
| PacBio HiFi reads | FASTA | 3 files | 50x | Zenodo (synthetic) |
| Illumina Hi-C Forward | FASTQ.GZ | 1 file | ~60x | Zenodo (SRR7126301_1) |
| Illumina Hi-C Reverse | FASTQ.GZ | 1 file | ~60x | Zenodo (SRR7126301_2) |
| Bionano optical maps | CMAP | 1 file | - | Zenodo |

**Data Source URLs:**
- HiFi reads: `https://zenodo.org/record/6098306/`
- Hi-C reads: `https://zenodo.org/record/5550653/`
- Bionano maps: `https://zenodo.org/records/5887339/`

---

## 5. Tools and Software

| Tool | Version | Purpose | Reference |
|---|---|---|---|
| **Cutadapt** | 4.4 | Adapter detection and read removal | Martin 2011 |
| **Meryl** | 1.3 | K-mer counting and database generation | Rhie et al. 2020 |
| **GenomeScope2** | 2.0 | Genome size and heterozygosity estimation | Ranallo-Benavidez et al. 2020 |
| **Hifiasm** | 0.19.8 | De novo assembly with Hi-C phasing | Cheng et al. 2021 |
| **gfastats** | 1.3.6 | Assembly graph manipulation and statistics | Formenti et al. 2022 |
| **BUSCO** | 5.5.0 | Gene completeness assessment | Simão et al. 2015 |
| **Merqury** | 1.3 | Reference-free quality value assessment | Rhie et al. 2020 |
| **Bionano Hybrid Scaffold** | 3.7.0 | Optical map-based scaffolding | Bionano Genomics |
| **BWA-MEM2** | 2.2.1 | Hi-C read alignment | Vasimuddin et al. 2019 |
| **Filter and merge** | 1.0 | Chimeric Hi-C read processing | Arima Genomics |
| **YaHS** | 1.2a.2 | Hi-C based chromosome scaffolding | Zhou et al. 2022 |
| **PretextMap** | 0.1.9 | Hi-C contact map generation | Harry 2022 |
| **Pretext Snapshot** | 0.0.3 | Contact map visualization | Harry 2022 |

---

## 6. Pipeline Stages

### Stage 1: HiFi Read Preprocessing (Cutadapt)
PacBio HiFi reads can contain adapter sequences that must be removed before assembly. Unlike NGS reads where adapters appear at read ends, HiFi adapters can appear anywhere within a read due to the nature of SMRT sequencing. Reads containing adapters are entirely discarded rather than trimmed, as their presence indicates lower quality data.

**Key parameters:**
- Maximum error rate: 0.1
- Minimum overlap: 35 bp
- Both adapter orientations checked
- Reads with adapters: discarded entirely

### Stage 2: Genome Profiling (Meryl + GenomeScope2)
Before assembly, we characterize the genome by analyzing k-mer frequencies in the raw reads. This provides estimates of genome size, heterozygosity, and repeat content, which guides downstream analysis.

**K-mer analysis workflow:**
1. Meryl counts all 31-mers in each HiFi file separately (parallelized)
2. Individual k-mer databases are merged using union-sum operation
3. A k-mer frequency histogram is generated
4. GenomeScope2 fits a statistical model to the histogram

**Key results from GenomeScope2:**
- Estimated haploid genome size: ~11.7 Mb
- Heterozygosity rate: 0.576%
- Diploid coverage: ~50x
- Haploid coverage peak: ~25x

### Stage 3: Genome Assembly (Hifiasm - Hi-C mode)
Hifiasm assembles the HiFi reads into contigs while using Hi-C data to phase the two haplotypes. The output is two separate assemblies: **Hap1** and **Hap2**, representing the two sets of chromosomes.

**How Hifiasm works:**
1. Builds an assembly graph from HiFi read overlaps
2. Identifies heterozygous bubbles in the graph
3. Uses Hi-C contact information to assign bubbles to haplotypes
4. Outputs two phased contig sets in GFA format

**GFA to FASTA conversion:**
The GFA (Graphical Fragment Assembly) format preserves graph information. gfastats converts GFA to FASTA for downstream tools.

### Stage 4: Assembly Quality Control

**gfastats - Assembly Statistics:**
Reports key metrics including:
- Number of contigs
- Total assembly length
- Largest contig size
- N50 (length at which 50% of assembly is in contigs of this size or larger)
- NG50 (N50 relative to estimated genome size)
- GC content

**BUSCO - Gene Completeness:**
Searches for universal single-copy orthologs expected to be present in the genome. Uses the Saccharomycetes lineage database.
- Complete single-copy: Good assembly
- Duplicated: Possible false duplications
- Fragmented/Missing: Possible assembly gaps

**Merqury - K-mer Quality Value:**
Compares k-mers in the reads to k-mers in the assembly to estimate:
- Quality Value (QV): Error rate of the assembly
- Completeness: Fraction of read k-mers present in assembly
- Phasing: How well k-mers are distributed between hap1 and hap2

### Stage 5: Bionano Optical Map Scaffolding
Bionano optical maps provide long-range physical information (up to megabases) that helps order and orient contigs into larger scaffolds. The technology labels specific sequence motifs along DNA molecules and images them to create a physical map.

**Hybrid scaffolding steps:**
1. Generate in silico maps from the contig sequences
2. Align in silico maps to Bionano optical maps
3. Identify and resolve conflicts between the two maps
4. Merge non-conflicting maps into hybrid scaffolds
5. Generate final AGP and FASTA scaffold files

### Stage 6: Hi-C Scaffolding (YaHS)
Hi-C data captures the 3D organization of chromatin in the nucleus. Regions that are physically close in 3D space (and therefore on the same chromosome) have higher interaction frequencies. This information is used to order and orient scaffolds into chromosome-level assemblies.

**Hi-C processing steps:**
1. Map forward and reverse Hi-C reads separately to the Bionano assembly using BWA-MEM2
2. Merge and filter chimeric reads using Arima protocol
3. Generate initial contact map with PretextMap (before YaHS)
4. Run YaHS to scaffold contigs into chromosomes
5. Remap Hi-C reads to YaHS scaffolds
6. Generate final contact map with PretextMap and Pretext Snapshot

**Why map Hi-C reads separately?**
Standard paired-end aligners assume a known insert size distribution. Hi-C ligation products can range from 1 bp to hundreds of megabases, violating this assumption. Therefore, each read is mapped independently.

---

## 7. Results Summary

### Assembly Statistics

| Metric | Hap1 | Hap2 |
|---|---|---|
| Number of contigs | 16 | 17 |
| Total length | ~11.3 Mb | ~12.2 Mb |
| Estimated genome size | 11.7 Mb | 11.7 Mb |

### Genome Profile (GenomeScope2)
| Parameter | Value |
|---|---|
| Estimated haploid size | ~11.7 Mb |
| Heterozygosity | 0.576% |
| Diploid coverage | ~50x |
| Model fit | >93% |

### Final Assembly Quality
| Metric | Result |
|---|---|
| Chromosome-level scaffolds | 16 (matching reference) |
| Contact map quality | Matches reference genome |
| Assembly completeness | High (BUSCO) |

### Key Result Images
- `results/genomescope/` - K-mer spectrum and genome model plots
- `results/busco/` - BUSCO completeness bar charts for hap1 and hap2
- `results/merqury/` - Spectra-CN and spectra-ASM k-mer plots
- `results/gfastats/` - Detailed assembly statistics table
- `results/pretext/` - Hi-C contact maps before and after YaHS scaffolding

---

## 8. Repository Structure

```
VGP-Genome-Assembly-Pipeline/
│
├── README.md                          ← This file
│
├── workflow/
│   └── VGP_assembly.ga                ← Galaxy workflow (importable)
│
├── results/
│   ├── genomescope/
│   │   ├── HiFi_synthetic_50x_01.png  ← Linear k-mer plot (sample 1)
│   │   ├── HiFi_synthetic_50x_02.png  ← Linear k-mer plot (sample 2)
│   │   ├── HiFi_synthetic_50x_03.png  ← Linear k-mer plot (sample 3)
│   │   └── *.txt                      ← GenomeScope summary files
│   │
│   ├── busco/
│   │   ├── busco_hap1_summary.png     ← BUSCO bar chart for hap1
│   │   ├── busco_hap1_summary.txt     ← BUSCO text summary for hap1
│   │   └── busco_hap2_summary.png     ← BUSCO bar chart for hap2
│   │
│   ├── merqury/
│   │   ├── spectra_cn_plot.png        ← Copy number spectrum
│   │   └── spectra_asm_plot.png       ← Assembly spectrum
│   │
│   ├── gfastats/
│   │   └── hap1_hap2_stats.txt        ← Assembly statistics table
│   │
│   └── pretext/
│       ├── contact_map_before_yahs.png ← Hi-C map before scaffolding
│       └── contact_map_after_yahs.png  ← Hi-C map after scaffolding
│
└── docs/
    └── methods.md                     ← Detailed step-by-step methods
```

---

## 9. How to Reproduce

### Requirements
- A Galaxy account at [usegalaxy.org](https://usegalaxy.org) or [usegalaxy.eu](https://usegalaxy.eu)
- Sufficient storage quota (the datasets are large)

### Steps

**1. Import the Galaxy Workflow**
- Log into Galaxy
- Go to **Workflow** menu
- Click **Import** 
- Upload the file: `workflow/VGP_assembly.ga`

**2. Upload Input Data**
Upload the following datasets to your Galaxy history:
```
HiFi reads (FASTA format):
https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_01.fasta
https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_02.fasta
https://zenodo.org/record/6098306/files/HiFi_synthetic_50x_03.fasta

Hi-C reads (FASTQSANGER.GZ format):
https://zenodo.org/record/5550653/files/SRR7126301_1.fastq.gz
https://zenodo.org/record/5550653/files/SRR7126301_2.fastq.gz

Bionano optical maps (CMAP format):
https://zenodo.org/records/5887339/files/bionano.cmap
```

**3. Organize Data**
- Create a collection from the 3 HiFi FASTA files named `HiFi data`
- Rename Hi-C files as `Hi-C_dataset_F` and `Hi-C_dataset_R`

**4. Run the Workflow**
- Go to **Workflow** menu
- Click **Run** on `VGP_assembly.ga`
- Select your input datasets
- Click **Run Workflow**

**5. Follow the Tutorial**
For detailed step-by-step guidance, refer to the official GTN tutorial:
> https://gxy.io/GTN:T00039

---

## 10. Key Concepts

| Term | Definition |
|---|---|
| **Contig** | A contiguous gapless sequence assembled from reads |
| **Scaffold** | One or more contigs joined by gap sequences (Ns) |
| **N50** | Contig length at which 50% of the assembly is covered |
| **Haplotype** | One copy of a chromosome pair |
| **Phasing** | Assigning contigs to their correct haplotype of origin |
| **K-mer** | A substring of length k from a sequence |
| **BUSCO** | Benchmarking Universal Single-Copy Orthologs |
| **Hi-C** | Chromatin conformation capture sequencing |
| **Optical map** | Physical map of DNA using labeled restriction sites |
| **False duplication** | Assembly error where one region appears twice |
| **Purging** | Removing falsely duplicated sequences from assembly |

---

## 11. References

1. Rhie, A., et al. (2021). Towards complete and error-free genome assemblies of all vertebrate species. *Nature*, 592, 737–746.
2. Cheng, H., et al. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods*, 18, 170–175.
3. Zhou, C., et al. (2022). YaHS: yet another Hi-C scaffolding tool. *Bioinformatics*, 39(1).
4. Ranallo-Benavidez, T. R., et al. (2020). GenomeScope 2.0 and Smudgeplots for reference-free profiling of polyploid genomes. *Nature Communications*, 11, 1432.
5. Rhie, A., et al. (2020). Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. *Genome Biology*, 21, 245.
6. Simão, F. A., et al. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. *Bioinformatics*, 31(19), 3210–3212.
7. Formenti, G., et al. (2022). Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs. *Bioinformatics*, 38(17).
8. Wenger, A. M., et al. (2019). Accurate circular consensus long-read sequencing improves variant detection and assembly of a human genome. *Nature Biotechnology*, 37, 1155–1162.

---

## Author
**Abid Hussain**  
Bioinformatics Pipeline Implementation  
Vertebrate Genome Project Tutorial - Galaxy Training Network  
April 2026
EOF
