# Detailed Methods: VGP Genome Assembly Pipeline

## Table of Contents
1. [Platform and Environment](#1-platform-and-environment)
2. [Data Acquisition](#2-data-acquisition)
3. [Step 1: HiFi Read Preprocessing](#3-step-1-hifi-read-preprocessing)
4. [Step 2: K-mer Counting with Meryl](#4-step-2-k-mer-counting-with-meryl)
5. [Step 3: Genome Profiling with GenomeScope2](#5-step-3-genome-profiling-with-genomescope2)
6. [Step 4: Genome Assembly with Hifiasm](#6-step-4-genome-assembly-with-hifiasm)
7. [Step 5: GFA to FASTA Conversion](#7-step-5-gfa-to-fasta-conversion)
8. [Step 6: Assembly Statistics with gfastats](#8-step-6-assembly-statistics-with-gfastats)
9. [Step 7: Completeness Assessment with BUSCO](#9-step-7-completeness-assessment-with-busco)
10. [Step 8: K-mer Quality Assessment with Merqury](#10-step-8-k-mer-quality-assessment-with-merqury)
11. [Step 9: Bionano Optical Map Scaffolding](#11-step-9-bionano-optical-map-scaffolding)
12. [Step 10: Hi-C Read Mapping](#12-step-10-hi-c-read-mapping)
13. [Step 11: Initial Hi-C Contact Map](#13-step-11-initial-hi-c-contact-map)
14. [Step 12: YaHS Chromosome Scaffolding](#14-step-12-yahs-chromosome-scaffolding)
15. [Step 13: Final Hi-C Contact Map](#15-step-13-final-hi-c-contact-map)
16. [Quality Control Summary](#16-quality-control-summary)

---

## 1. Platform and Environment

All analyses were performed using the **Galaxy** bioinformatics platform (usegalaxy.org/usegalaxy.eu), a web-based open-source platform that allows accessible and reproducible bioinformatics analysis without requiring command-line expertise.

**Galaxy History Name:** VGP_Genome_Assembly  
**Analysis Date:** 2026  
**Organism:** *Saccharomyces cerevisiae* S288C (synthetic diploid reads)

---

## 2. Data Acquisition

### 2.1 PacBio HiFi Reads
Three synthetic HiFi read files were uploaded to Galaxy in FASTA format. These reads simulate 50x coverage of a diploid yeast genome. After upload, all three files were organized into a **Galaxy dataset collection** named `HiFi data` to allow parallel processing.

| File | Format | Dataset # |
|---|---|---|
| HiFi_synthetic_50x_01.fasta | FASTA | Collection item 1 |
| HiFi_synthetic_50x_02.fasta | FASTA | Collection item 2 |
| HiFi_synthetic_50x_03.fasta | FASTA | Collection item 3 |

**Source:** `https://zenodo.org/record/6098306/`

### 2.2 Illumina Hi-C Reads
Two paired-end Hi-C read files were uploaded in FASTQSANGER.GZ format and renamed for clarity.

| File | Renamed to | Dataset # |
|---|---|---|
| SRR7126301_1.fastq.gz | Hi-C_dataset_F | 99 |
| SRR7126301_2.fastq.gz | Hi-C_dataset_R | 100 |

**Source:** `https://zenodo.org/record/5550653/`

### 2.3 Bionano Optical Maps
One Bionano CMAP file was uploaded in CMAP format.

| File | Format | Dataset # |
|---|---|---|
| bionano.cmap | CMAP | 230 |

**Source:** `https://zenodo.org/records/5887339/`

---

## 3. Step 1: HiFi Read Preprocessing

### Tool: Cutadapt (Galaxy version 4.4+galaxy0)

### Why this step?
PacBio HiFi reads are generated using SMRT sequencing, where a DNA polymerase reads a circular template multiple times. This process can introduce adapter sequences at unpredictable positions within the reads (not just at the ends). Reads containing adapter sequences are of lower quality and must be removed before assembly.

Unlike standard adapter trimming (which cuts the adapter off), we **discard the entire read** if an adapter is found, because their presence signals the read is unreliable.

### Parameters Used

| Parameter | Value | Reason |
|---|---|---|
| Read type | Single-end | HiFi reads are treated as single-end |
| Input | HiFi_collection | All 3 HiFi FASTA files as a collection |
| Adapter 1 name | First adapter | PacBio SMRTbell adapter sequence |
| Adapter 1 sequence | ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT | Forward adapter |
| Adapter 2 name | Second adapter | PacBio SMRTbell adapter (reverse) |
| Adapter 2 sequence | ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT | Reverse adapter |
| Maximum error rate | 0.1 | Allow up to 10% mismatches when detecting adapters |
| Minimum overlap | 35 bp | Adapter must overlap read by at least 35 bp to be flagged |
| Check reverse complement | Yes | Check both orientations of adapter |
| Discard trimmed reads | Yes | Remove entire read if adapter found |

### Output
- **Dataset #105:** `HiFi_collection (trimmed)` — A collection of 3 adapter-free FASTA files

---

## 4. Step 2: K-mer Counting with Meryl

### Tool: Meryl (Galaxy version 1.3+galaxy6)

### Why this step?
Before assembling the genome, we analyze the frequency distribution of all k-mers (substrings of length k) in the raw reads. This k-mer spectrum gives us information about:
- Estimated genome size
- Level of heterozygosity
- Repeat content
- Sequencing error rate
- Sequencing depth/coverage

### Sub-step 2a: Count K-mers per File

Each HiFi file in the collection is processed separately (parallelized) to generate individual k-mer databases.

| Parameter | Value | Reason |
|---|---|---|
| Operation | Count: count occurrences of canonical k-mers | Count how many times each k-mer appears |
| Input | HiFi_collection (trimmed) | Use adapter-removed reads |
| K-mer size | 31 | Long enough to be unique, short enough to tolerate errors |

**Output:** `meryldb` — A collection of 3 k-mer count databases (Dataset #109)

**Why k=31?**
- 31-mers are long enough that most are unique in a ~12 Mb genome
- Short enough to be robust to sequencing errors
- Standard choice for genomes under 10 Gb

### Sub-step 2b: Merge K-mer Databases

The three individual databases are merged into one using union-sum operation, which combines counts across all files.

| Parameter | Value |
|---|---|
| Operation | Union-sum: return k-mers in any input, sum the counts |
| Input | Collection meryldb (all 3 databases) |

**Output:** `Merged meryldb` (Dataset #113)

### Sub-step 2c: Generate K-mer Histogram

A frequency histogram is generated from the merged database for use by GenomeScope2.

| Parameter | Value |
|---|---|
| Operation | Generate histogram dataset |
| Input | Merged meryldb |

**Output:** `meryldb histogram` (Dataset #117)

---

## 5. Step 3: Genome Profiling with GenomeScope2

### Tool: GenomeScope (Galaxy version 2.0+galaxy2)

### Why this step?
GenomeScope2 uses the k-mer frequency histogram to fit a statistical model and estimate key genome properties. It uses nonlinear least-squares optimization to fit a mixture of negative binomial distributions to the histogram.

### Parameters Used

| Parameter | Value | Reason |
|---|---|---|
| Input histogram | meryldb histogram | The k-mer frequency histogram from Meryl |
| Ploidy | 2 | Diploid organism (two copies of each chromosome) |
| K-mer length | 31 | Must match the k-mer size used in Meryl |
| Output summary | Yes | Get text summary of estimates |
| Create testing.tsv | Yes | Save model parameters |

### Outputs (Datasets #121–#126)
| Output | Description |
|---|---|
| Linear plot | K-mer frequency vs coverage |
| Log plot | Logarithmic version of linear plot |
| Transformed linear plot | Frequency × coverage vs coverage (enhances higher-order peaks) |
| Transformed log plot | Log version of transformed plot |
| Model | Detailed model fitting report |
| Summary | Estimated genome properties |

### Key Results Interpretation
The k-mer profile shows a **bimodal distribution**:
- **First peak (~25x coverage):** Heterozygous k-mers (present in only one haplotype)
- **Second peak (~50x coverage):** Homozygous k-mers (present in both haplotypes)

| Estimated Parameter | Value |
|---|---|
| Haploid genome size | ~11.7 Mb |
| Heterozygosity rate | 0.576% |
| Model fit | >93% |
| Diploid coverage | ~50x |

---

## 6. Step 4: Genome Assembly with Hifiasm

### Tool: Hifiasm (Galaxy version 0.19.8+galaxy0)

### Why this step?
Hifiasm is a state-of-the-art de novo assembler specifically designed for PacBio HiFi reads. In **Hi-C phased mode**, it uses Hi-C contact information to separate the two haplotypes during assembly, producing two fully phased assemblies (hap1 and hap2).

### How Hifiasm works:
1. Builds a string graph from all-vs-all HiFi read overlaps
2. Identifies heterozygous bubbles (regions where haplotypes differ)
3. Uses Hi-C read pair information to phase bubbles — reads that map to the same chromosome will have more Hi-C contacts
4. Outputs separate graphs for each haplotype

### Parameters Used

| Parameter | Value | Reason |
|---|---|---|
| Assembly mode | Standard | Standard diploid assembly |
| Input reads | HiFi_collection (trimmed) | Use adapter-removed HiFi reads |
| Hi-C partition | Specify | Enable Hi-C phasing mode |
| Hi-C R1 reads | Hi-C_dataset_F | Forward Hi-C reads |
| Hi-C R2 reads | Hi-C_dataset_R | Reverse Hi-C reads |

### Outputs (Datasets #148, #149)
| Output | Dataset # | Description |
|---|---|---|
| Hi-C hap1 balanced contig graph | 148 | GFA graph for haplotype 1 |
| Hi-C hap2 balanced contig graph | 149 | GFA graph for haplotype 2 |

Both outputs are renamed with tags:
- `Hap1 contigs graph` (#hap1 tag)
- `Hap2 contigs graph` (#hap2 tag)

---

## 7. Step 5: GFA to FASTA Conversion

### Tool: gfastats (Galaxy version 1.3.6+galaxy0)

### Why this step?
Hifiasm outputs assembly graphs in **GFA (Graphical Fragment Assembly)** format, which encodes both sequences and their relationships (overlaps). Most downstream QC tools require standard **FASTA** format. gfastats converts GFA to FASTA while preserving all sequence information.

### Parameters Used

| Parameter | Value |
|---|---|
| Input GFA files | Hap1 contigs graph AND Hap2 contigs graph |
| Tool mode | Genome assembly manipulation |
| Output format | FASTA |
| Generate initial set of paths | Yes (toggled on) |

### Outputs (Datasets #171, #172)
- `Hap1 contigs FASTA` (#hap1 tag)
- `Hap2 contigs FASTA` (#hap2 tag)

---

## 8. Step 6: Assembly Statistics with gfastats

### Tool: gfastats (Galaxy version 1.3.6+galaxy0)

### Why this step?
Before and after each assembly step, we collect standard assembly statistics to monitor quality and improvement. gfastats provides comprehensive metrics about the assembly.

### Sub-step 6a: Generate Statistics

| Parameter | Value |
|---|---|
| Input files | Hap1 contigs graph AND Hap2 contigs graph (GFA) |
| Tool mode | Summary statistics generation |
| Expected genome size | 11747160 (from GenomeScope2) |
| Thousands separator | No |

**Outputs:** `Hap1 stats` (Dataset #159), `Hap2 stats` (Dataset #160)

### Sub-step 6b: Join Hap1 and Hap2 Statistics

Tool: **Column join (Galaxy version 0.0.3)**
- Joins Hap1 stats and Hap2 stats side by side for easy comparison
- **Output:** `gfastats on hap1 and hap2 (full)` (Dataset #161)

### Sub-step 6c: Filter to Contigs Only

Tool: **Search in textfiles (Galaxy version 1.1.1)**
- Removes lines containing "scaffold" (not relevant at contig stage)
- Match type: Don't match / case insensitive / basic regex
- **Output:** `gfastats on hap1 and hap2 contigs` (Dataset #162)

### Key Statistics Explained

| Statistic | Definition |
|---|---|
| # contigs | Total number of contigs |
| Total length | Sum of all contig lengths |
| Largest contig | Length of the longest contig |
| N50 | Half the assembly is in contigs of this size or larger |
| NG50 | N50 calculated relative to estimated genome size |
| GC content | Percentage of G and C bases |

---

## 9. Step 7: Completeness Assessment with BUSCO

### Tool: Busco (Galaxy version 5.5.0+galaxy0)

### Why this step?
BUSCO (Benchmarking Universal Single-Copy Orthologs) assesses assembly completeness by searching for genes that are expected to be present exactly once in a complete genome of a given taxonomic group. The presence, absence, or duplication of these genes reveals assembly quality.

### Parameters Used

| Parameter | Value | Reason |
|---|---|---|
| Sequences to analyze | Hap1 contigs FASTA AND Hap2 contigs FASTA | Evaluate both haplotypes |
| Lineage data source | Use cached lineage data | Use pre-downloaded database |
| Mode | Genome assemblies (DNA) | We are assessing genome assembly |
| Gene predictor | Use Metaeuk | Modern gene predictor for eukaryotes |
| Lineage selection | Select lineage | Manually choose database |
| Lineage | Saccharomycetes | Yeast-specific BUSCO database |
| Outputs | Short summary text AND summary image | Both text and visual output |

### Outputs (Datasets #173–#175 for hap1, #176–#178 for hap2)
| Output | Description |
|---|---|
| Short summary text | Percentage breakdown of BUSCO gene categories |
| Summary image | Visual bar chart of BUSCO results |
| Full table | Detailed per-gene results |

### BUSCO Categories Explained
| Category | Meaning |
|---|---|
| Complete Single-copy (C,S) | Gene found once — ideal |
| Complete Duplicated (C,D) | Gene found multiple times — possible false duplication |
| Fragmented (F) | Gene partially found — possible assembly gap |
| Missing (M) | Gene not found — possible missing sequence |

---

## 10. Step 8: K-mer Quality Assessment with Merqury

### Tool: Merqury (Galaxy version 1.3+galaxy3)

### Why this step?
Merqury provides a reference-free way to assess assembly quality by comparing k-mers in the raw reads to k-mers in the assembly. It calculates:
- **QV (Quality Value):** Estimated base-level accuracy of the assembly
- **Completeness:** What fraction of read k-mers are present in the assembly
- **Spectra plots:** Visual representation of k-mer distribution across assemblies

### Parameters Used

| Parameter | Value |
|---|---|
| Evaluation mode | Default mode |
| K-mer counts database | Merged meryldb |
| Number of assemblies | Two assemblies |
| First genome assembly | Hap1 contigs FASTA |
| Second genome assembly | Hap2 contigs FASTA |

### Outputs (Datasets #179–#181)
| Output | Description |
|---|---|
| Stats collection | Completeness statistics per assembly |
| Plots collection | Spectra-CN and spectra-ASM plots |
| QV stats collection | Quality value statistics |

### Plot Interpretation

**Spectra-CN plot (spectra_cn_plot.png):**
- Grey region (left): K-mers only in reads — sequencing errors or missing assembly sequence
- Red area: K-mers present once in assembly — heterozygous regions
- Blue area: K-mers present twice in assembly — homozygous regions or duplications

**Spectra-ASM plot (spectra_asm_plot.png):**
- Shows which k-mers are in hap1 only, hap2 only, or shared
- Large shared peak at diploid coverage (50x): Good — homozygous regions correctly shared
- Haploid peaks split between hap1 and hap2: Good — heterozygous regions correctly phased

---

## 11. Step 9: Bionano Optical Map Scaffolding

### Tool: Bionano Hybrid Scaffold (Galaxy version 3.7.0+galaxy3)

### Why this step?
After contig assembly, we have many separate contigs that need to be ordered and oriented relative to each other. Bionano optical maps provide long-range physical information that can span gaps between contigs, allowing them to be joined into larger scaffolds.

### What are Bionano optical maps?
Bionano's Saphyr system labels specific sequence motifs along long native DNA molecules (hundreds of kilobases to megabases in length) using fluorescent labels. The pattern of labels creates a unique "barcode" that can be mapped to the genome. This physical map is independent of sequencing and can resolve structural variations and repetitive regions.

### Parameters Used

| Parameter | Value | Reason |
|---|---|---|
| NGS FASTA | Hap1 contigs FASTA | The assembly to scaffold |
| BioNano CMAP | bionano.cmap (Dataset #230) | The optical map |
| Configuration mode | VGP mode | Optimized settings for VGP pipeline |
| Genome maps conflict filter | Cut contig at conflict | Resolve conflicts by cutting |
| Sequences conflict filter | Cut contig at conflict | Resolve conflicts by cutting |

### Outputs (Datasets #231–#235)
| Output | Description |
|---|---|
| NGScontigs scaffold NCBI trimmed (#231) | Contigs successfully scaffolded with Bionano |
| NGScontigs not scaffolded trimmed (#232) | Contigs that could not be scaffolded |
| Hybrid scaffold report (#233) | Summary of scaffolding results |
| Conflicts (#234) | Detected conflicts between sequence and optical map |
| AGP file (#235) | Assembly description file |

### Concatenation Step
Tool: **Concatenate datasets**
- Combines scaffolded contigs (#231) and unscaffolded contigs (#232)
- **Output:** `Hap1 assembly bionano` (Dataset #236) — Complete assembly with Bionano scaffolding

### Bionano Statistics (gfastats)
After scaffolding, gfastats is run again on the Bionano assembly to compare statistics before and after scaffolding.
- **Output:** `Bionano stats` (Dataset #237)

---

## 12. Step 10: Hi-C Read Mapping

### Tool: BWA-MEM2 (Galaxy version 2.2.1+galaxy1)

### Why this step?
Before using Hi-C data for scaffolding, the Hi-C reads must be aligned to the current assembly. Hi-C reads are mapped **separately** (forward and reverse reads in independent jobs) because the insert size distribution of Hi-C ligation products cannot be predicted (it ranges from 1 bp to hundreds of megabases), unlike standard paired-end sequencing.

### Forward Read Mapping

| Parameter | Value |
|---|---|
| Reference genome | Hap1 assembly bionano (Dataset #236) |
| Build index | Yes (build from history) |
| Read type | Single-end |
| Input reads | Hi-C_dataset_F (forward reads) |
| Set read groups | Do not set |
| Analysis mode | Simple Illumina mode |
| BAM sorting | Sort by read names (QNAME) |

**Output:** `BAM forward` (Dataset #240)

### Reverse Read Mapping

| Parameter | Value |
|---|---|
| Reference genome | Hap1 assembly bionano (Dataset #236) |
| Build index | Yes (build from history) |
| Read type | Single-end |
| Input reads | Hi-C_dataset_R (reverse reads) |
| Set read groups | Do not set |
| Analysis mode | Simple Illumina mode |
| BAM sorting | Sort by read names (QNAME) |

**Output:** `BAM reverse` (Dataset #241)

### Why sort by read name?
Sorting by read name (QNAME) groups forward and reverse reads from the same Hi-C pair together, which is required for the next filtering step.

### Merging and Filtering

Tool: **Filter and merge chimeric reads from Arima Genomics (Galaxy version 1.0+galaxy1)**

This tool follows the Arima Genomics Hi-C data processing protocol:
- Merges forward and reverse BAM files
- Filters out chimeric reads (reads that span ligation junctions)
- Produces a clean BAM file ready for contact map generation

| Parameter | Value |
|---|---|
| First set of reads | BAM forward (Dataset #240) |
| Second set of reads | BAM reverse (Dataset #241) |

**Output:** `BAM Hi-C reads` (Dataset #242)

---

## 13. Step 11: Initial Hi-C Contact Map

### Tools: PretextMap + Pretext Snapshot

### Why this step?
Before running YaHS scaffolding, we generate a Hi-C contact map to visualize the current state of the assembly. This baseline map helps us evaluate the improvement after scaffolding.

### PretextMap (Galaxy version 0.1.9+galaxy0)

| Parameter | Value |
|---|---|
| Input BAM | BAM Hi-C reads (Dataset #242) |
| Sort by | Don't sort |

**Output:** `PretextMap output` (Dataset #243)

### Pretext Snapshot (Galaxy version 0.0.3+galaxy1)

| Parameter | Value |
|---|---|
| Input Pretext map | PretextMap output (Dataset #243) |
| Output format | PNG |
| Show grid | Yes |

**Output:** `Pretext Snapshot on dataset 243` (Dataset #244) — Collection of 18 PNG images

### Contact Map Interpretation
- **Blocks along the diagonal:** High contact frequency within chromosomes (expected)
- **Off-diagonal signals:** Contacts between different chromosomes or misassemblies
- **17 blocks visible:** Corresponding to 17 scaffolds in Bionano assembly
- **Signals near diagonal:** Indicate correctly ordered contigs

---

## 14. Step 12: YaHS Chromosome Scaffolding

### Tool: YaHS (Galaxy version 1.2a.2+galaxy1)

### Why this step?
YaHS (Yet Another Hi-C Scaffolding tool) uses the Hi-C contact information to order and orient scaffolds into chromosome-level assemblies. Unlike other Hi-C scaffolding tools, YaHS does not require the user to specify the number of chromosomes.

### How YaHS works:
1. Breaks contigs at positions with insufficient Hi-C coverage (possible assembly errors)
2. Creates a contact matrix by splitting contigs into chunks
3. Calculates joining scores based on normalized contact frequencies between contig pairs
4. Builds a scaffolding graph with contigs as nodes
5. Simplifies the graph (removes ambiguous edges, resolves bubbles)
6. Traverses the graph to assemble chromosomal scaffolds
7. Performs multiple rounds at decreasing resolution (increasing chunk sizes)

### Parameters Used

| Parameter | Value | Reason |
|---|---|---|
| Input contig sequences | Hap1 assembly bionano (Dataset #236) | The Bionano-scaffolded assembly |
| Alignment file | BAM Hi-C reads (Dataset #242) | Hi-C alignments to contigs |
| Restriction enzyme | Enter specific sequence | Use actual enzyme sequence |
| Restriction enzyme sequence | CTTAAG | HindIII recognition site (used in Hi-C library prep) |

### Outputs (Datasets #263–#267)
| Output | Description |
|---|---|
| AGP initial break files (#263) | Assembly description before scaffolding |
| AGP break files (#264) | Assembly description of breaks made |
| AGP scaffolding results (#265) | Assembly description of scaffolding |
| Final scaffolds AGP (#266) | Final assembly description |
| YaHS Scaffolds FASTA (#267) | **Final chromosome-level assembly in FASTA format** |

---

## 15. Step 13: Final Hi-C Contact Map

### Tools: BWA-MEM2 → Filter and merge → PretextMap → Pretext Snapshot

### Why this step?
After YaHS scaffolding, we remap the Hi-C reads to the final scaffolded assembly to generate a new contact map. Comparing this map to the pre-YaHS map allows us to evaluate the improvement in chromosome-level organization.

### Forward Read Remapping (BWA-MEM2)

| Parameter | Value |
|---|---|
| Reference | YaHS Scaffolds FASTA (Dataset #267) |
| Read type | Single-end |
| Input | Hi-C_dataset_F |
| BAM sorting | Sort by read names |

**Output:** `BAM forward YaHS` (Dataset #275)

### Reverse Read Remapping (BWA-MEM2)

| Parameter | Value |
|---|---|
| Reference | YaHS Scaffolds FASTA (Dataset #267) |
| Read type | Single-end |
| Input | Hi-C_dataset_R |
| BAM sorting | Sort by read names |

**Output:** `BAM reverse YaHS` (Dataset #276)

### Filter and Merge

| Parameter | Value |
|---|---|
| First set of reads | BAM forward YaHS (Dataset #275) |
| Second set of reads | BAM reverse YaHS (Dataset #276) |

**Output:** `BAM Hi-C reads YaHS` (Dataset #277)

### PretextMap

| Parameter | Value |
|---|---|
| Input BAM | BAM Hi-C reads YaHS (Dataset #277) |
| Sort by | Don't sort |

**Output:** `PretextMap output YaHS` (Dataset #278)

### Pretext Snapshot

| Parameter | Value |
|---|---|
| Input Pretext map | PretextMap output YaHS (Dataset #278) |
| Output format | PNG |
| Show grid | Yes |

**Output:** `Pretext Snapshot YaHS` (Dataset #279) — Final contact map collection

### Final Contact Map Interpretation
Comparing before (Dataset #244) and after (Dataset #279) YaHS:
- After YaHS: Cleaner diagonal blocks representing individual chromosomes
- Inversions resolved: Off-diagonal signals are corrected
- 16 chromosomal scaffolds: Matching the reference S. cerevisiae genome
- Contact pattern matches the reference genome contact map

---

## 16. Quality Control Summary

| Stage | Tool | Key Metrics |
|---|---|---|
| Raw reads | Cutadapt | Reads with adapters removed |
| Genome profile | GenomeScope2 | Size ~11.7 Mb, Het 0.576% |
| Contig assembly | gfastats | N50, contig count, total length |
| Gene completeness | BUSCO | % complete, duplicated, missing |
| K-mer quality | Merqury | QV score, completeness % |
| Bionano scaffolding | gfastats | Scaffold N50, scaffold count |
| Final assembly | Pretext Snapshot | Chromosome-level contact map |

### Final Assembly Comparison with Reference

| Metric | Our Assembly | S. cerevisiae Reference |
|---|---|---|
| Total length | ~11.7 Mb | ~12.1 Mb |
| Chromosomes | 16 | 16 + mitochondrial |
| Assembly quality | Chromosome-level | Chromosome-level |
| Contact map | Matches reference | Reference |

The final assembled genome shows near-identical structure to the S288C reference genome, demonstrating the effectiveness of the VGP pipeline in producing high-quality chromosome-level assemblies.
EOF
