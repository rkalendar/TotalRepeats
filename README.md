# TotalRepeats

**an integrated genome-wide platform for rapid de novo identification, classification, annotation, comparative analysis, and visualization of repetitive elements**

[![Java](https://img.shields.io/badge/Java-25+-orange.svg?logo=openjdk)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-blue.svg)]()
[![Online Tool](https://img.shields.io/badge/Try%20Online-TotalRepeats-green.svg)](https://primerdigital.com/tools/repeats.html)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)

---

## Table of Contents

- [Overview](#overview)
- [What's New](#whats-new)
- [Key Features](#key-features)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Basic Syntax](#basic-syntax)
  - [Memory Configuration](#memory-configuration)
  - [Common Options](#common-options)
  - [Advanced Options](#advanced-options)
- [Detailed Option Reference](#detailed-option-reference)
- [Input/Output Formats](#inputoutput-formats)
- [Example Workflows](#example-workflows)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Author & Contact](#author--contact)
- [Citation](#citation)
- [License](#license)
- [Related Tools](#related-tools)

---

## Overview

TotalRepeats is a universal, *de novo* tool for genome-wide identification, classification, annotation, visualization, and comparison of repetitive DNA elements. It efficiently detects a wide spectrum of repeats, including:

- **Mobile genetic elements** — transposons, retrotransposons (LINEs, SINEs, LTR elements, DNA transposons)
- **Tandem arrays** — microsatellites (SSRs), minisatellites, telomeric repeats, centromeric satellite DNA
- **Gene duplication** — duplicated genes, rearrangements, and copy-number variants
- **Large-scale structural variations** — segmental duplications, rearrangements
- **Unique sequence identification** — detection of sequences that occur only once across the entire genome and within individual chromosomes
 
The tool is particularly well-suited for comparative genomics, evolutionary biology, structural variation studies, and bioinformatics research. Masking, clustering, and annotation are fully **multithreaded**, scaling efficiently from personal laptops to HPC/supercomputing clusters.

---

## What's New

Latest additions to TotalRepeats:

- **Comparative analyses are no longer capped at ~2 GB.** The `-combine`, `-combine2`, and `-combinemask` modes previously failed with `OutOfMemoryError: Requested string length exceeds VM limit` once the *combined* length of all input sequences passed ~2.1 GB — Java's hard limit on the size of a single array/`String`. They now concatenate the inputs **virtually** and address them with **64-bit (long) coordinates**, so a joint analysis scales to pangenome-sized data (tens of Gb) without ever materializing a single oversized string. Each individual sequence/chromosome is still limited to 2 GB; the *total* across files is now effectively bounded only by available RAM.
- **Per-file reports in every comparative mode.** `-combine`, `-combine2`, and `-combinemask` now each emit an individual GFF3 annotation, PNG, and SVG per input file (named after that file), *in addition to* the combined report. Each genome can be inspected on its own while the joint clustering is preserved — the same family keeps the same `ClusterID`, colour, row, and reference label in both the per-file and combined views.
- **Pangenome analysis — core / accessory / unique — across all comparative modes.** Every combine run now writes a pangenome report classifying each repeat family by how many sequences it occurs in: **core** (present in all), **accessory** (present in some), and **unique** (present in one). It includes the family-frequency spectrum, per-sequence statistics, a pairwise shared-family matrix (Jaccard similarity), and a machine-readable presence/absence matrix. See [`-combine`](#-combine--comparative-analysis).
- **Combined image is now SVG.** The joint (cross-file) visualization is written as scalable SVG, which has no pixel-size ceiling and stays sharp at pangenome scale. Per-file images are still produced as both PNG and SVG.
- The parallel clustering is now **deterministic and reproducible** — repeated runs on the same input produce identical cluster assignments and `ClusterID`s, independent of how many worker threads are used.
- **Performance improvements.** Faster low-complexity/SSR masking and sequence clustering — k-mer encoding and reduced memory traffic in the hot loops — with identical results.

---

## Key Features

| Feature | Description |
|---|---|
| 🔍 **Rapid *de novo* detection** | K-mer-based identification and classification of all repeat classes — no prior repeat library required |
| 🧬 **Comparative analysis** | Analyse multiple genomes or assemblies simultaneously with synchronized repeat clustering |
| 🔀 **Polymorphism detection** | Identify interspecific repeat polymorphisms without whole-genome alignment |
| 📚 **External annotation** | Annotate repeats against Repbase, Dfam/FamDB, or custom species-specific libraries |
| 📊 **Publication-ready visualization** | Generate PNG and SVG images of repeat landscapes, annotated maps, and comparative views |
| ⚡ **Multithreaded performance** | Parallel processing across all cores for masking, clustering, and annotation; multithreaded clustering is **deterministic** — identical, reproducible results on every run |
| 🧩 **Pangenome analysis** | Core / accessory / unique repeat-family classification with a presence/absence matrix across genomes — produced by every comparative mode (`-combine`, `-combine2`, `-combinemask`) |
| 🗄️ **Pangenome-scale inputs** | Combined comparative analysis is not limited by the 2 GB single-sequence cap — inputs are concatenated virtually and addressed with 64-bit coordinates |

---

## How It Works

```
┌──────────────────────────────────────────────────────────────────────────┐
│                       TotalRepeats Pipeline                              │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  1. INPUT                    2. MASKING                                  │
│  ┌────────────────┐          ┌─────────────────────────────┐             │
│  │ FASTA files    │────────▶ │ K-mer decomposition        │             │
│  │ (single or dir)│          │ Identify repeated k-mers    │             │
│  └────────────────┘          │ Merge overlapping regions   │             │
│                              └─────────────┬───────────────┘             │
│                                            │                             │
│  3. CLUSTERING                             ▼                             │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │ • Group repeat copies into families by sequence identity │            │
│  │ • Classify: tandem (STR), unclassified (UCRP),           │            │
│  │   classified (CRP) via external library                  │            │
│  │ • Multithreaded k-mer vector comparison                  │            │
│  └─────────────────────────┬────────────────────────────────┘            │
│                            │                                             │
│  4. ANNOTATION             ▼               5. OUTPUT                     │
│  ┌──────────────────────────┐    ┌──────────────────────────┐            │
│  │ Match clusters against   │───▶│ GFF3 annotations        │            │
│  │ Repbase / Dfam / custom  │    │ Masked FASTA             │            │
│  │ library (optional)       │    │ PNG / SVG visualizations │            │
│  └──────────────────────────┘    │ Summary statistics       │            │
│                                  └──────────────────────────┘            │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

1. **Input** — Reads one or more FASTA files (or an entire directory of genomes).
2. **Masking** — Decomposes sequences into k-mers, identifies over-represented k-mers, and masks repeated regions.
3. **Clustering** — Groups repeat copies into families by sequence similarity; classifies each family as tandem repeat, unclassified interspersed, or annotated (if a library is provided).
4. **Annotation** *(optional)* — Matches repeat clusters against Repbase, Dfam/FamDB, or a custom database to assign superfamily and class labels.
5. **Output** — Produces GFF3 annotations, soft-masked FASTA, publication-ready images, and summary statistics.

---

## Requirements

| Requirement | Details |
|---|---|
| **Java** | Version **25** or higher ([download](https://www.oracle.com/java/technologies/downloads/)) |
| **OS** | Windows, Linux, or macOS — platform-independent |
| **RAM** | Depends on genome size (see [Memory Configuration](#memory-configuration)) |

> **Verify your Java version:**
> ```bash
> java -version
> ```
> If you see a version below 25, update Java or install via [Conda](#option-2-using-conda).

**Setting `JAVA_HOME` / `PATH`:** <https://www.java.com/en/download/help/path.html>

---

## Installation

### Option 1: Direct Download

```bash
# 1. Clone the repository (or download the JAR directly)
git clone https://github.com/rkalendar/TotalRepeats.git
cd TotalRepeats

# 2. Confirm Java 25+ is available
java -version

# 3. Run — no build step required
java -jar dist/TotalRepeats.jar genome.fasta -out=/results/output
```

### Option 2: Using Conda

If Java 25 is not installed system-wide:

```bash
# Create a dedicated environment
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n java25 openjdk=25
conda activate java25

# Verify
java -version
```

---

## Quick Start

```bash
# Show usage / help (also printed when no input file is given)
java -jar TotalRepeats.jar --help

# Analyse a single genome
java -jar TotalRepeats.jar genome.fasta

# Analyse all genomes in a directory
java -jar TotalRepeats.jar /path/to/genomes/

# Large genome with memory allocation
java -Xms32g -Xmx64g -jar TotalRepeats.jar large_genome.fasta

# Annotate repeats with an external library
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.ref

# Compare multiple genomes side by side
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -combine
```

---

## Usage

### Basic Syntax

```bash
java [JVM flags] -jar TotalRepeats.jar <input_file_or_directory> [options]
```

- **Input** can be a single FASTA file or a directory containing multiple genomes.
- **Options** are space-separated and follow the input path.
- No additional dependencies are required.
- **Need help?** Run with `-help`, `--help`, `-h`, `-?`, `/?`, or `/h` (or with no arguments at all) to print the full usage guide and exit. The usage guide is also shown if the supplied input file or folder does not exist.

### Memory Configuration

When working with large genomes, allocate additional heap memory using JVM flags:

| Chromosome Size | Recommended RAM | JVM Flags |
|---|---|---|
| < 100 MB | Default | None needed |
| 100–500 MB | 16 GB | `-Xms8g -Xmx16g` |
| 500 MB – 2 GB | 32 GB | `-Xms16g -Xmx32g` |
| > 2 GB | 64+ GB | `-Xms32g -Xmx64g` |

**What the flags mean:**

- **`-Xms`** (initial heap) — Memory pre-allocated at startup. Avoids costly resizing during execution.
- **`-Xmx`** (maximum heap) — Upper memory limit. Prevents `OutOfMemoryError` for large genomes.

> **Tip:** Set `-Xms` to roughly one-quarter of `-Xmx` for a good balance between startup speed and memory efficiency.

> **Comparative runs (`-combine`, `-combine2`, `-combinemask`):** size the heap to the *combined* length of all input files, not the largest single file. These modes now handle a combined size above 2 GB (they no longer build one giant string), so the only practical limit is the heap you provide — allocate accordingly for pangenome-scale inputs.

### Common Options

| Option | Description | Default |
|---|---|---|
| `kmer=N` | K-mer size for repeat detection (range: 9–21) | `19` |
| `sln=N` | Minimum repeat block length (bp) | `80` |
| `image=WxH` | Output image dimensions (pixels) | Auto |
| `imgx=N` | Image width scaling factor (1 = compressed, 20 = stretched) | `5` |
| `flanks=N` | Extend each repeat by N bases on both sides | `0` |
| `-seqshow` | Include repeat sequences in GFF3 output | Off |
| `-maskonly` | Generate only the masked FASTA (skip clustering/annotation) | Off |
| `-normal` | Preventing multithreaded clustering | Off |
| `-help` | Show the usage guide and exit (aliases: `--help`, `-h`, `-?`, `/?`, `/h`) | — |

### Advanced Options

| Option | Description |
|---|---|
| `-combine` | Comparative analysis — synchronized repeat clustering across multiple files; per-file reports + combined report + pangenome (see [Comparing the comparative modes](#comparing-the-comparative-modes)) |
| `-combine2` | Comparative analysis where all target sequences are treated as a single combined search space — repeat identification is performed over the whole set, then clustered jointly. Produces per-file reports and a pangenome report (see [Comparing the comparative modes](#comparing-the-comparative-modes)) |
| `-combinemask` | Comparative analysis with mask-based detection — synchronized clustering across files for masking-mode comparison and pangenome studies; per-file reports + combined report + pangenome (no per-file masked FASTA) |
| `-homology` | Mask all homologous regions to highlight unique sequences between genomes |
| `-amask` | Perform masking via pairwise sequence alignment (instead of k-mer-based) |
| `-readmask` | Import existing mask files for clustering, annotation, and visualization |
| `-readgff` | Import GFF annotation files for direct visualization |
| `-extract` | Split a multi-entry FASTA into individual files (one per sequence) |
| `-maskscomp` | Compare masked outputs from different tools or pipelines |
| `-out=<path>` | Path to output folder (default: current folder) |
| `-lib=<path>` | Annotate repeats with an external library (Repbase, Dfam, or custom FASTA) |

---

## Detailed Option Reference

### `kmer=` — K-mer Size

Sets the k-mer length for repeat masking (range: 9–21). Shorter k-mers detect shorter repeats but increase runtime and noise; longer k-mers are more specific but miss short elements.

| K-mer Range | Best For | Examples |
|---|---|---|
| 9–12 | Very short repeats | CRISPR spacers |
| 13–17 | Small chromosomes, organellar genomes | Mitochondria, plastids, bacteria |
| 18–21 | Eukaryotic chromosomes (recommended) | Human, plant, fungal genomes |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta kmer=21
```

### `sln=` — Minimum Repeat Block Length

Sets the minimum length (bp) for a repeat region to be included in the output. Shorter values capture small elements (e.g., Alu ~300 bp, tRNA-SINEs ~100 bp) at the cost of more fragmented results.

| Value | Effect |
|---|---|
| 20–60 | Captures very short repeats (CRISPR spacers, MITE elements) |
| 60–80 | Good balance for most analyses (default) |
| 100+ | Filters out short elements; focuses on longer repeats |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta sln=100
```

### Multithreaded Clustering

Enables multithreaded clustering. The pairwise k-mer-vector comparison that dominates the clustering step is spread across CPU cores, giving a large speed-up on genomes with many repeat families. 
The parallel clustering is **deterministic**: repeated runs on the same input produce identical cluster assignments and `ClusterID`s, and the result does not depend on the number of worker threads. Output is therefore fully reproducible — across repeated runs and across machines with different core counts — which matters for published analyses. Changing the thread count affects only runtime, not the result.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta 
```

**Controlling the number of threads.** Multithreaded clustering uses Java's common Fork/Join worker pool, which by default draws on all available CPU cores. To cap the worker-thread count — for example on a shared cluster node, or to leave cores free for other jobs — set the standard JVM property *before* `-jar`:

```bash
# Limit clustering to 8 worker threads
java -Djava.util.concurrent.ForkJoinPool.common.parallelism=8 -Xms16g -Xmx32g -jar \
    TotalRepeats.jar genome.fasta
```

Setting the value to 1 forces single-threaded execution or use the `-normal` flag to prevent multithreaded clustering.

```bash
# Limit clustering to 1 worker threads
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -normal
```

### `-lib=` — External Repeat Library

Annotate detected repeat families against a curated database. TotalRepeats does **not** ship with libraries — download them separately:

| Source | URL | Notes |
|---|---|---|
| **Dfam / FamDB** | <https://www.dfam.org/releases/current/families/FamDB> | Open access; HMM-based profiles |
| **Repbase** | <https://www.girinst.org/> | Requires registration; curated consensus sequences |
| **Custom** | — | Build your own FASTA library for a species of interest |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.ref
```

### Comparing the comparative modes

All three comparative modes cluster repeats *jointly* across the input files, and all three now produce per-file reports, a combined report, and a pangenome report. They differ mainly in how repeats are detected and in which per-file files are written:

| Mode | Repeat detection | Combined image | Per-file outputs | Pangenome | Typical use |
|---|---|---|---|---|---|
| `-combine` | k-mer masking **per sequence**, then joint clustering | SVG | GFF3 · PNG · SVG · masked FASTA | ✅ | Compare homologous chromosomes / assemblies |
| `-combine2` | sequences treated as **one** combined search space for detection | SVG | GFF3 · PNG · SVG · masked FASTA | ✅ | Treat a set of sequences as a single search space |
| `-combinemask` | mask-based detection per sequence, then joint clustering | SVG | GFF3 · PNG · SVG | ✅ | Masking-mode comparison / pangenome studies |

> All three handle a *combined* input larger than 2 GB — see [What's New](#whats-new). The combined (cross-file) image is SVG in every mode; per-file images are PNG **and** SVG.

### `-combine` — Comparative Analysis

Performs synchronized repeat clustering across multiple sequences. All files in the input directory are processed together, enabling direct cross-genome comparison of repeat families.

**Ideal for:**
- Comparing homologous chromosomes across species or strains
- Detecting repeat polymorphisms without whole-genome alignment
- Analysing long-read assemblies (e.g., ONT, PacBio HiFi)

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -combine
```

**Output files produced by a `-combine` run:**

| File | Description |
|---|---|
| `<file>.gff`, `<file>.png`, `<file>.svg` | Per-file annotation and images for each input sequence |
| `<file>.msk` | Per-file soft-masked FASTA |
| `<report>.gff`, `<report>.svg` | Combined annotation and image (SVG) across all sequences |
| `<report>_pangenome.txt` | Pangenome summary (see below) |
| `<report>_pangenome.tsv` | Repeat-family presence/absence matrix |

**Pangenome report.** Using the synchronized clustering, each repeat family (CRP cluster) is mapped back to the sequences it occurs in and classified as:

- **core** — present in *all* analyzed sequences (shared content);
- **accessory** (shell) — present in some but not all sequences;
- **unique** (cloud) — present in a single sequence (differing content).

`<report>_pangenome.txt` reports the core/accessory/unique breakdown, the family-frequency spectrum (number of families present in exactly *k* sequences), per-sequence statistics (CRP / UCRP / STR content), the genome-specific (unique) families, and a pairwise shared-family matrix with Jaccard similarity. `<report>_pangenome.tsv` is a machine-readable presence/absence matrix — one row per family, one column per sequence — suitable for downstream analysis and plotting. The `ClusterID` column matches the GFF3 `ClusterID`, so the reports cross-reference each other.

The pangenome report is generated automatically for any multi-sequence comparative run (`-combine`, `-combine2`, or `-combinemask`).

### `-combinemask` — Pangenome Analysis

Performs pangenome-scale comparative analysis with mask-based detection. Enables synchronized clustering and annotation across assemblies, and supports comparing different masking strategies or parameters.

Like `-combine`, it now also writes an **individual GFF3 + PNG + SVG per input file**, a **combined GFF + SVG**, and a **pangenome report** (`_pangenome.txt` / `_pangenome.tsv`). It does not write per-file masked FASTA. See [Comparing the comparative modes](#comparing-the-comparative-modes).

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -combinemask
```

### `-homology` — Homology Masking

Masks all homologous regions (both within and between input sequences) to reveal **unique, non-shared sequence**. This inverts the usual repeat-finding question: instead of "what's repeated?", it answers "what's unique?"

**Use cases:**
- Identifying lineage-specific sequence between closely related genomes
- Comparing human chromosomes X and Y
- Discovering strain-specific insertions or deletions

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -homology
```

### `imgx=` — Image Width Scaling

Controls the horizontal resolution of output images (PNG and SVG). Higher values produce wider, more detailed images suitable for zooming into specific regions.

| Value | Effect |
|---|---|
| `1` | Maximum compression — overview of entire genome |
| `10` | Default — good for most chromosome-level views |
| `20–30` | Stretched — detailed inspection of individual regions |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta imgx=20
```

### `-readmask` — Import Existing Masks

Imports pre-computed masking files (e.g., from RepeatMasker or a previous TotalRepeats run) for clustering, annotation, and visualization without re-running the masking step.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar masked_genome.fasta -readmask
```

### `-readgff` — Import GFF for Visualization

Imports an existing GFF annotation file and generates visualizations directly, bypassing repeat detection.

```bash
java -Xms16g -jar TotalRepeats.jar annotation.gff -readgff
```

### `-extract` — Split Multi-Entry FASTA

Splits a multi-entry FASTA file into separate files (one per sequence record). Useful for preparing input when you need per-chromosome processing.

```bash
java -jar TotalRepeats.jar multi_entry.fasta -extract
```

### `-maskscomp` — Compare Masking Tools

Compares masked outputs from different software tools (e.g., TotalRepeats vs. RepeatMasker). Requires a **tab-delimited** control file listing sequence/mask pairs:

```bash
java -jar TotalRepeats.jar runfiles.txt -maskscomp
```

**Example `runfiles.txt`:**

```
NC_017844.1.fasta	NC_017844.1.fasta.msk
NC_017849.1.fasta	NC_017849.1.fasta.msk
NC_017850.1.fasta	NC_017850.1.fasta.msk
```

---

## Input/Output Formats

### Input: FASTA

TotalRepeats accepts standard FASTA format or plain-text nucleotide sequences. Multi-entry FASTA files are supported. Each individual sequence (chromosome / record) may be up to 2 GB; in the comparative modes (`-combine`, `-combine2`, `-combinemask`) the *combined* size across all files is not limited by that cap.

```
>sequence_id Optional description
ATCGATCGATCGATCGATCGATCG...
```

### Output Overview

| Output File | Format | Description |
|---|---|---|
| **Repeat annotations** | GFF3 | Coordinates, classification, and cluster assignments for every repeat element |
| **Masked genome** | FASTA | Soft-masked sequence (repeats in lowercase) |
| **Repeat landscape** | PNG / SVG | Publication-ready visualization of repeat distribution |
| **Summary statistics** | Text | Genome-wide repeat content, class breakdown, and family counts |
| **Per-file reports** *(combine modes)* | GFF3 / PNG / SVG | Individual annotation and images for each input file, alongside the combined output |
| **Pangenome report** *(combine modes)* | Text + TSV | Core / accessory / unique repeat-family classification and presence/absence matrix |

### GFF3 Output Format

Results are saved in [GFF3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), a nine-column, tab-delimited standard:

| Column | Field | Description |
|---|---|---|
| 1 | **Seqid** | Source file or sequence identifier |
| 2 | **Type** | Repeat class: `STR` (tandem), `UCRP` (unclassified interspersed), or `CRP` (classified/annotated) |
| 3 | **ClusterID** | Repeat family cluster identifier |
| 4 | **Start** | Start position (1-based) |
| 5 | **End** | End position (inclusive) |
| 6 | **Length** | Repeat element length (bp) |
| 7 | **Strand** | `+` (forward) or `-` (complement) |
| 8 | **Phase** | Not used (`.`) |
| 9 | **Attributes** | Repeat sequence (if `-seqshow` is enabled) |

---

## Example Workflows

### Basic Repeat Analysis

<details>
<summary><strong>Windows</strong></summary>

```bash
# Single file
java -jar C:\TotalRepeats\TotalRepeats.jar C:\Genomes\NC_014637.fasta

# Custom k-mer and minimum length
java -jar C:\TotalRepeats\TotalRepeats.jar C:\Genomes\ kmer=16 sln=30

# Extract repeat sequences with 100 bp flanking regions
java -jar C:\TotalRepeats\TotalRepeats.jar C:\Genomes\NC_014637.fasta -seqshow flanks=100
```

</details>

<details>
<summary><strong>Linux / macOS</strong></summary>

```bash
# Single file
java -jar TotalRepeats.jar /data/genomes/NC_014637.fasta

# Custom k-mer and minimum length
java -jar TotalRepeats.jar /data/genomes/ kmer=16 sln=30

# Extract repeat sequences with 100 bp flanking regions
java -jar TotalRepeats.jar /data/genomes/NC_014637.fasta -seqshow flanks=100
```

</details>

### Large Genome Analysis

```bash
# Human T2T genome with Repbase annotation
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/T2T-CHM13v2.0/ \
    -lib=/data/libraries/humsub.ref

# Tasmanian devil genome
java -Xms16g -Xmx64g -jar TotalRepeats.jar /data/genomes/Sarcophilus_harrisii/ -out=/results/output

# Very large genome (Iberian ribbed newt, ~20 Gb) — stretched image output
java -Xms16g -Xmx64g -jar TotalRepeats.jar /data/genomes/Pleurodeles_waltl/ imgx=20
```

### Comparative Genomics

```bash
# Compare bacterial strains — synchronized repeat clustering
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -combine

# Compare bacterial strains — synchronized repeat clustering
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -combine2

# Compare bacterial strains — synchronized repeat clustering
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -combinemask

# Identify unique regions between related fungal genomes
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Pyricularia_oryzae/ -homology
```

### Benchmarking Against Other Tools

```bash
# Compare TotalRepeats masking vs. RepeatMasker output
java -jar TotalRepeats.jar comparison_list.txt -maskscomp
```

---

## Troubleshooting

### No input file (usage is printed instead of results)

If you launch TotalRepeats with no arguments, with an explicit help flag, or with a first argument that is not an existing file or folder, the program prints the usage guide and exits without analysing anything:

```bash
# Explicit help request — prints usage and exits
java -jar TotalRepeats.jar --help

# Missing input (options only) — prints usage
java -jar TotalRepeats.jar -kmer=18 -seqshow

# Wrong path — prints "input file or folder not found" + usage
java -jar TotalRepeats.jar /no/such/genome.fasta
```

Make sure the **input file or folder is the first argument**, before any options:

```bash
java -jar TotalRepeats.jar genome.fasta -kmer=18 -seqshow   # correct
```

The recognised help flags are `-help`, `--help`, `-h`, `-?`, `/?`, `/h`, and `/help` (case-insensitive).

### Java version errors

Ensure Java 25+ is installed and active:

```bash
java -version
```

If multiple Java versions are installed, set `JAVA_HOME` explicitly:

```bash
export JAVA_HOME=/path/to/jdk-25
export PATH="$JAVA_HOME/bin:$PATH"
```

### `OutOfMemoryError`

Increase heap allocation according to the [Memory Configuration](#memory-configuration) table:

```bash
java -Xms16g -Xmx64g -jar TotalRepeats.jar genome.fasta
```

> **Tip:** If your machine has limited RAM, try `-maskonly` to skip clustering — the masking step uses significantly less memory.

> **`Requested string length exceeds VM limit` in a comparative run:** this specific error came from the old `-combine`/`-combine2`/`-combinemask` code building a single string of the whole concatenation, which Java caps at ~2.1 GB. It is resolved — these modes now concatenate virtually and use 64-bit coordinates. If you still hit a plain `OutOfMemoryError` (heap), raise `-Xmx` per the table above.

### No repeats detected

- **Check k-mer size:** For small genomes or short repeats, reduce `kmer` (e.g., `kmer=12`).
- **Lower minimum length:** Reduce `sln` (e.g., `sln=30`) to capture shorter elements.
- **Verify input:** Ensure the FASTA file contains valid nucleotide sequence.

### Output images are too small or too compressed

Increase the image scaling factor:

```bash
java -jar TotalRepeats.jar genome.fasta imgx=30
```

Or specify exact dimensions:

```bash
java -jar TotalRepeats.jar genome.fasta image=4000x2000
```

### Slow clustering on large genomes

Enable multithreaded clustering:

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta
```

By default this uses all available CPU cores. To cap the worker-thread count (e.g. on a shared node), add `-Djava.util.concurrent.ForkJoinPool.common.parallelism=N` before `-jar`. The clustering result is deterministic, so the thread count changes only the runtime, not the output.

---

## FAQ

**Q: Does TotalRepeats require a pre-built repeat library?**
No. TotalRepeats performs *de novo* detection using k-mer analysis — no external library is needed. However, you can optionally provide a library (Repbase, Dfam, or custom) via `-lib=` to classify and annotate the detected repeat families.

**Q: What organisms does TotalRepeats support?**
Any organism with a DNA sequence. TotalRepeats is organism-agnostic and has been used on viral, bacterial, fungal, plant, and animal genomes.

**Q: How does TotalRepeats compare to RepeatMasker?**
RepeatMasker relies on a curated library (Repbase/Dfam) and performs homology-based detection — it can only find repeats represented in its database. TotalRepeats works *de novo*, detecting all repeated sequences regardless of prior annotation. The two tools are complementary: use TotalRepeats for comprehensive discovery and RepeatMasker for precise classification. The `-maskscomp` option lets you compare their outputs directly.

**Q: Can I process multiple genomes at once?**
Yes. Point TotalRepeats at a directory containing FASTA files, and it will process all of them. Add `-combine` (or `-combine2` / `-combinemask`) for synchronized cross-genome repeat clustering — each of these produces a per-file report and image for every input, a combined report and image, and a pangenome report (core / accessory / unique families plus a presence/absence matrix). See [Comparing the comparative modes](#comparing-the-comparative-modes).

**Q: What is the maximum genome size TotalRepeats can handle?**
TotalRepeats has been tested on genomes exceeding 20 Gb (e.g., *Pleurodeles waltl*). A single sequence (one FASTA record) is capped at 2 GB by Java's array/`String` limit, but in the comparative modes the *combined* size across files is not — those modes concatenate virtually and use 64-bit coordinates. The practical limiting factor is available RAM, so allocate sufficient heap via `-Xms` / `-Xmx`.

**Q: My `-combine` run crashed with `Requested string length exceeds VM limit` — is that fixed?**
Yes. That error came from concatenating all inputs into one string, which Java caps at ~2.1 GB. The comparative modes (`-combine`, `-combine2`, `-combinemask`) now build the concatenation virtually and address it with 64-bit coordinates, so the combined analysis runs at pangenome scale. Any remaining memory issue would be an ordinary heap `OutOfMemoryError` — raise `-Xmx`.

**Q: Can I import results from other tools?**
Yes. Use `-readmask` to import masked FASTA files from RepeatMasker or other tools, and `-readgff` to import GFF annotations for visualization.

**Q: What does the `-quick` flag do exactly?**
It enables multithreaded pairwise comparison during the clustering step, which significantly speeds up processing for large genomes with many repeat families. The parallel clustering is deterministic — repeated runs produce identical cluster assignments and `ClusterID`s regardless of the thread count, so results stay reproducible. The flag `-fast` is an alias for `-quick`. To cap the number of worker threads, add `-Djava.util.concurrent.ForkJoinPool.common.parallelism=N` before `-jar` (see the `-quick` / `-fast` entry under [Detailed Option Reference](#detailed-option-reference)).

---

## Author & Contact

**Ruslan Kalendar**
📧 ruslan.kalendar@helsinki.fi

🌐 **Online version:** <https://primerdigital.com/tools/repeats.html>
💻 **Source code:** <https://github.com/rkalendar/TotalRepeats>

---

## Citation

If you use TotalRepeats in your research, please cite the software and the online resource:

> Kalendar, R. *TotalRepeats: an integrated genome-wide platform for rapid de novo identification, classification, annotation, comparative analysis, and visualization of repetitive elements* [Computer software]. University of Helsinki. <https://github.com/rkalendar/TotalRepeats>

```bibtex
@software{Kalendar_TotalRepeats,
  author = {Kalendar, Ruslan},
  title  = {{TotalRepeats: an integrated genome-wide platform for de novo identification, classification, annotation, comparative analysis, and visualization of repetitive elements}},
  url    = {https://github.com/rkalendar/TotalRepeats},
  note   = {Online version: https://primerdigital.com/tools/repeats.html}
}
```

If a peer-reviewed article describing TotalRepeats is available, please cite it as well.

---

## License

This project is distributed under the terms of the [GNU General Public License v3.0](LICENSE.txt).

---

## Related Tools

| Tool | Description |
|---|---|
| **[Repeater](https://github.com/rkalendar/Repeater)** | Pairwise sequence alignment tool |
| **[NCBI RefSeq GenBank Downloader](https://github.com/rkalendar/genbanktools)** | Batch-download GenBank records by gene or accession |

---

<p align="center">
  <em>TotalRepeats — Genome-wide repeat discovery, from laptop to supercomputer.</em>
</p>
