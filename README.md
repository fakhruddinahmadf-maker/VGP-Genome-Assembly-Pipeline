# Chromosome-Level Genome Assembly using the VGP Pipeline

A Galaxy-based bioinformatics project implementing the Vertebrate Genome Project (VGP) workflow to produce a haplotype-phased, chromosome-scale genome assembly from multi-platform sequencing data.

---

## About This Project

**Organism:** *Saccharomyces cerevisiae* S288C (Baker's Yeast)  
**Platform:** Galaxy (usegalaxy.org)  
**Assembly strategy:** HiFi + Hi-C + Bionano optical maps (VGP Trajectory D)  
**Data:** Synthetic PacBio HiFi reads simulating a diploid yeast genome

The goal was to reconstruct both haplotypes of a diploid genome from scratch, validate quality at multiple checkpoints, and produce chromosome-level scaffolds — without relying on any reference genome during assembly.

---

## Why These Three Technologies?

| Technology | What It Contributes |
|---|---|
| PacBio HiFi | Long, accurate reads (10–25 kb, >99% accuracy) for contig assembly |
| Bionano optical maps | Mega-base scale physical information for ordering contigs |
| Illumina Hi-C | 3D chromatin contact data for chromosome-level scaffolding |

Each technology addresses a different limitation. HiFi reads span repetitive regions that short reads cannot resolve. Optical maps bridge gaps between distant contigs. Hi-C data anchors scaffolds to chromosome-scale structures.

---

## Pipeline at a Glance

```
HiFi reads (3 files)
    │
    ▼
Cutadapt ──── Remove adapter-containing reads
    │
    ▼
Meryl + GenomeScope2 ──── K-mer profiling → genome size & heterozygosity
    │
    ▼
Hifiasm (Hi-C mode) ──── Phased contig assembly → Hap1 + Hap2
    │
    ▼
QC checkpoint: gfastats / BUSCO / Merqury
    │
    ▼
Bionano Hybrid Scaffold ──── Optical map integration
    │
    ▼
BWA-MEM2 → YaHS ──── Hi-C scaffolding to chromosome level
    │
    ▼
PretextMap / Pretext Snapshot ──── Final contact map visualisation
```

---

## Input Data

All datasets are publicly available on Zenodo.

| File | Format | Description |
|---|---|---|
| HiFi_synthetic_50x_01–03.fasta | FASTA | PacBio HiFi reads, 50x coverage |
| SRR7126301_1.fastq.gz | FASTQ.GZ | Hi-C forward reads |
| SRR7126301_2.fastq.gz | FASTQ.GZ | Hi-C reverse reads |
| bionano.cmap | CMAP | Bionano optical map |

**Download links:**
- HiFi: https://zenodo.org/record/6098306/
- Hi-C: https://zenodo.org/record/5550653/
- Bionano: https://zenodo.org/records/5887339/

---

## Tools Used

| Tool | Version | Role |
|---|---|---|
| Cutadapt | 4.4 | Adapter removal |
| Meryl | 1.3 | K-mer counting |
| GenomeScope2 | 2.0 | Genome size estimation |
| Hifiasm | 0.19.8 | De novo phased assembly |
| gfastats | 1.3.6 | Assembly statistics |
| BUSCO | 5.5.0 | Gene completeness |
| Merqury | 1.3 | Quality value assessment |
| Bionano Hybrid Scaffold | 3.7.0 | Optical map scaffolding |
| BWA-MEM2 | 2.2.1 | Hi-C alignment |
| YaHS | 1.2a.2 | Chromosome scaffolding |
| PretextMap + Snapshot | 0.1.9 / 0.0.3 | Contact map generation |

---

## Key Results

**Genome profile (GenomeScope2)**
- Estimated haploid size: ~11.7 Mb
- Heterozygosity: 0.576%
- Coverage: ~50x diploid (~25x per haplotype)

**Contig assembly (Hifiasm)**
- Hap1: 16 contigs, ~11.3 Mb total
- Hap2: 17 contigs, ~12.2 Mb total

**Final assembly**
- 16 chromosome-level scaffolds (matches the expected yeast karyotype)
- High BUSCO completeness scores for both haplotypes
- Hi-C contact maps show clean diagonal banding consistent with correct chromosome assignment

---

## Repository Structure

```
.
├── README.md
├── workflow/
│   └── VGP_assembly.ga          ← Importable Galaxy workflow
├── results/
│   ├── genomescope/             ← K-mer plots and summary files
│   ├── busco/                   ← BUSCO bar charts (hap1 and hap2)
│   ├── merqury/                 ← Spectra-CN and spectra-ASM plots
│   ├── gfastats/                ← Assembly statistics table
│   └── pretext/                 ← Hi-C contact maps (before/after YaHS)
└── docs/
    └── methods.md               ← Detailed step-by-step methods
```

---

## How to Run

**Requirements:** A Galaxy account at usegalaxy.org or usegalaxy.eu with sufficient storage quota.

1. **Import the workflow** — Go to *Workflow → Import* and upload `workflow/VGP_assembly.ga`

2. **Upload input data** — Add all six input files to your Galaxy history. Create a collection from the three HiFi FASTA files named `HiFi data`. Rename the Hi-C files to `Hi-C_dataset_F` and `Hi-C_dataset_R`.

3. **Run** — Go to *Workflow → Run*, select your inputs, and click *Run Workflow*.

For a full step-by-step walkthrough, follow the official Galaxy Training Network tutorial: https://gxy.io/GTN:T00039

---

## Glossary

**Contig** — A gapless sequence assembled directly from reads.  
**Scaffold** — Contigs joined together with gap sequences (Ns).  
**N50** — The contig length at which 50% of the total assembly is covered.  
**Phasing** — Separating reads/contigs into their correct haplotype of origin.  
**BUSCO** — A benchmark measuring how many expected single-copy genes are present and complete.  
**Hi-C** — A sequencing method that captures physical contacts between genomic regions in the nucleus.

---

## References

- Rhie et al. (2021). *Nature*, 592, 737–746.
- Cheng et al. (2021). *Nature Methods*, 18, 170–175.
- Zhou et al. (2022). *Bioinformatics*, 39(1).
- Ranallo-Benavidez et al. (2020). *Nature Communications*, 11, 1432.
- Rhie et al. (2020). *Genome Biology*, 21, 245.
- Simão et al. (2015). *Bioinformatics*, 31(19), 3210–3212.
