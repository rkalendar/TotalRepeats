# TotalRepeats

**A universal genome-wide tool for rapid *de novo* identification, classification, annotation, comparative analysis, and visualization of repetitive elements**

[![Java](https://img.shields.io/badge/Java-25+-orange.svg?logo=openjdk)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-blue.svg)]()
[![Online Tool](https://img.shields.io/badge/Try%20Online-TotalRepeats-green.svg)](https://primerdigital.com/tools/repeats.html)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)

---

## Table of Contents

- [Overview](#overview)
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
- **Genes duplication** — genes duplication, rearrangement and copy-number variants
- **Large-scale structural variations** — segmental duplications, rearrangements
- **Unique sequence identification** — detection of sequences that occur only once across the entire genome and within individual chromosomes
 
The tool is particularly well-suited for comparative genomics, evolutionary biology, structural variation studies, and bioinformatics research. Masking, clustering, and annotation are fully **multithreaded**, scaling efficiently from personal laptops to HPC/supercomputing clusters.

---

## Key Features

| Feature | Description |
|---|---|
| 🔍 **Rapid *de novo* detection** | K-mer-based identification and classification of all repeat classes — no prior repeat library required |
| 🧬 **Comparative analysis** | Analyse multiple genomes or assemblies simultaneously with synchronized repeat clustering |
| 🔀 **Polymorphism detection** | Identify interspecific repeat polymorphisms without whole-genome alignment |
| 📚 **External annotation** | Annotate repeats against Repbase, Dfam/FamDB, or custom species-specific libraries |
| 📊 **Publication-ready visualization** | Generate PNG and SVG images of repeat landscapes, annotated maps, and comparative views |
| ⚡ **Multithreaded performance** | Parallel processing across all cores for masking, clustering, and annotation |
| 🧩 **Pangenome support** | Cross-tool benchmarking and pangenome-scale repeat comparisons via `-combinemask` |

---

## How It Works

```
┌──────────────────────────────────────────────────────────────────────────┐
│                       TotalRepeats Pipeline                             │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  1. INPUT                    2. MASKING                                  │
│  ┌────────────────┐          ┌─────────────────────────────┐             │
│  │ FASTA files    │────────▶│ K-mer decomposition         │             │
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
git clone https://github.com/<username>/TotalRepeats.git
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
| `-quick` | Enable multithreaded clustering for faster processing | Off |

### Advanced Options

| Option | Description |
|---|---|
| `-combine` | Pangenome comparative analysis — synchronized repeat clustering across multiple files |
| `-homology` | Mask all homologous regions to highlight unique sequences between genomes |
| `-combinemask` | Compare multiple masking files for cross-tool benchmarking or pangenome studies |
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

### `-combine` — Comparative Analysis

Performs synchronized repeat clustering across multiple sequences. All files in the input directory are processed together, enabling direct cross-genome comparison of repeat families.

**Ideal for:**
- Comparing homologous chromosomes across species or strains
- Detecting repeat polymorphisms without whole-genome alignment
- Analysing long-read assemblies (e.g., ONT, PacBio HiFi)

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -combine
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

### `-combinemask` — Pangenome Analysis

Performs pangenome-scale comparative analysis using multiple masking files. Enables synchronized clustering and annotation across assemblies, and supports benchmarking different masking tools or parameters.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -combinemask
```

### `imgx=` — Image Width Scaling

Controls the horizontal resolution of output images (PNG and SVG). Higher values produce wider, more detailed images suitable for zooming into specific regions.

| Value | Effect |
|---|---|
| `1` | Maximum compression — overview of entire genome |
| `5` | Default — good for most chromosome-level views |
| `10–20` | Stretched — detailed inspection of individual regions |

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

TotalRepeats accepts standard FASTA format or plain-text nucleotide sequences. Multi-entry FASTA files (up to 2 GB per file) are supported.

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
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -quick
```

---

## FAQ

**Q: Does TotalRepeats require a pre-built repeat library?**
No. TotalRepeats performs *de novo* detection using k-mer analysis — no external library is needed. However, you can optionally provide a library (Repbase, Dfam, or custom) via `-lib=` to classify and annotate the detected repeat families.

**Q: What organisms does TotalRepeats support?**
Any organism with a DNA sequence. TotalRepeats is organism-agnostic and has been used on viral, bacterial, fungal, plant, and animal genomes.

**Q: How does TotalRepeats compare to RepeatMasker?**
RepeatMasker relies on a curated library (Repbase/Dfam) and performs homology-based detection — it can only find repeats represented in its database. TotalRepeats works *de novo*, detecting all repeated sequences regardless of prior annotation. The two tools are complementary: use TotalRepeats for comprehensive discovery and RepeatMasker for precise classification. The `-maskscomp` option lets you compare their outputs directly.

**Q: Can I process multiple genomes at once?**
Yes. Point TotalRepeats at a directory containing FASTA files, and it will process all of them. Add `-combine` for synchronized cross-genome repeat clustering.

**Q: What is the maximum genome size TotalRepeats can handle?**
TotalRepeats has been tested on genomes exceeding 20 Gb (e.g., *Pleurodeles waltl*). The limiting factor is available RAM — allocate sufficient heap memory via `-Xms` and `-Xmx` flags.

**Q: Can I import results from other tools?**
Yes. Use `-readmask` to import masked FASTA files from RepeatMasker or other tools, and `-readgff` to import GFF annotations for visualization.

**Q: What does the `-quick` flag do exactly?**
It enables multithreaded pairwise comparison during the clustering step. This significantly speeds up processing for large genomes with many repeat families, at no cost to accuracy.

---

## Author & Contact

**Ruslan Kalendar**
📧 ruslan.kalendar@helsinki.fi

🌐 **Online version:** <https://primerdigital.com/tools/repeats.html>


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
