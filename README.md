# TotalRepeats

**A universal genome-wide tool for rapid *de novo* identification, classification, annotation, comparative analysis, and visualization of repetitive elements**

[![Java](https://img.shields.io/badge/Java-25+-orange.svg)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-blue.svg)]()
[![Online Tool](https://img.shields.io/badge/Try%20Online-TotalRepeats-green.svg)](https://primerdigital.com/tools/repeats.html)

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
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
- [Author & Contact](#author--contact)
- [Citation](#citation)

---

## Overview

TotalRepeats is a universal, *de novo* tool for genome-wide identification, classification, annotation, visualization, and comparison of repetitive DNA elements. It efficiently detects a wide spectrum of repeats, including:

- **Mobile genetic elements** — transposons, retrotransposons
- **Tandem arrays** — microsatellites, minisatellites, telomeres, satellite DNA
- **Large-scale structural variations** — duplications, rearrangements

The tool is particularly well-suited for comparative genomics, evolutionary biology, structural variation studies, and bioinformatics research. Masking, clustering, and annotation are fully multithreaded, scaling efficiently from personal computers to HPC/supercomputing clusters.

---

## Key Features

| Feature | Description |
|---------|-------------|
| 🔍 **Rapid Identification** | Fast *de novo* detection and classification of repetitive elements |
| 🧬 **Comparative Analysis** | Analyze multiple genomes or assemblies simultaneously |
| 🔀 **Polymorphism Detection** | Identify interspecific polymorphisms without whole-genome alignment |
| 📚 **External Annotation** | Annotate repeats using Repbase, Dfam/FamDB, or custom datasets |
| 📊 **Visualization** | Generate publication-ready images (PNG/SVG) |
| ⚡ **Multithreaded** | Efficient parallel processing for large-scale analyses |

---

## Requirements

| Requirement | Details |
|-------------|---------|
| **Operating System** | Windows, Linux, or macOS |
| **Java Version** | Java 25 or higher |
| **RAM** | Depends on genome size (see [Memory Configuration](#memory-configuration)) |

**Download Java:** https://www.oracle.com/java/technologies/downloads/

**Set Java Path:** https://www.java.com/en/download/help/path.html

---

## Installation

### Option 1: Direct Download

1. Download `TotalRepeats.jar` from the `dist` directory
2. Place it in your preferred location
3. Ensure Java 25+ is installed and in your PATH

### Option 2: Using Conda

```bash
# Add conda-forge channel and set priority
conda config --add channels conda-forge
conda config --set channel_priority strict

# Create environment with OpenJDK 25
conda create -n java25 openjdk=25

# Activate environment
conda activate java25

# Verify installation
java -version
```

---

## Quick Start

```bash
# Basic analysis of a single FASTA file
java -jar TotalRepeats.jar genome.fasta

# Analyze all genomes in a directory
java -jar TotalRepeats.jar /path/to/genomes/

# Large genome with memory allocation
java -Xms32g -Xmx64g -jar TotalRepeats.jar large_genome.fasta
```

---

## Usage

### Basic Syntax

```bash
java -jar TotalRepeats.jar <input_file_or_directory> [options]
```

Input can be a single sequence file (FASTA/plain text) or a directory containing multiple genomes. No additional dependencies are required.

### Memory Configuration

| Genome Size | Recommended RAM | JVM Flags |
|-------------|-----------------|-----------|
| < 100 MB | Default | None needed |
| 100–500 MB | 16–32 GB | `-Xms8g -Xmx16g` |
| 500 MB – 2 GB | 64 GB | `-Xms32g -Xmx64g` |
| > 2 GB | 128+ GB | `-Xms64g -Xmx128g` |

**JVM Memory Parameters:**

- **`-Xms`** (Initial Heap Size): Memory allocated at startup. Helps avoid allocation delays for large genomes.
- **`-Xmx`** (Maximum Heap Size): Upper limit of heap memory. Set this to prevent `OutOfMemoryError`.

**Example for large genomes:**
```bash
java -Xms32g -Xmx128g -jar TotalRepeats.jar /data/genomes/large_genome/
```

### Common Options

| Option | Description | Default |
|--------|-------------|---------|
| `kmer=N` | k-mer size for repeat detection (9–21) | 19 |
| `sln=N` | Minimum repeat block length | 80 |
| `image=WxH` | Output image dimensions | Auto |
| `imgx=N` | Image width scaling (1=compressed, 20=stretched) | 5 |
| `flanks=N` | Extend repeats by N bases on each side | 0 |
| `-seqshow` | Include repeat sequences in output | — |
| `-maskonly` | Generate only the masked output file | — |

### Advanced Options

| Option | Description |
|--------|-------------|
| `-combine` | Pangenome comparative analysis with synchronized repeat clustering across multiple files |
| `-homology` | Mask homologous regions to highlight unique sequences between genomes |
| `-combinemask` | Compare multiple masking files for cross-tool benchmarking |
| `-quick` | Enable multithreaded clustering for faster processing |
| `-amask` | Perform masking via pairwise sequence alignment |
| `-readmask` | Import existing mask files for clustering/annotation |
| `-readgff` | Import GFF files for visualization |
| `-extract` | Split multi-entry FASTA into separate files |
| `-maskscomp` | Compare masked files from different tools |
| `-lib=PATH` | Use external repeat library (Repbase, Dfam, or custom) |

---

## Detailed Option Reference

### `kmer=` — K-mer Size

Sets the minimum k-mer value for repeat sequence masking (range: 9–21).

- Use **9–12** for very short repeats (e.g., CRISPR arrays)
- Use **18** for short chromosomes
- Use **19–21** for eukaryotic chromosomes (recommended)

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta kmer=21
```

### `sln=` — Minimum Sequence Length

Sets the minimum length for sequences included in the analysis.

- Use **60–80** to detect short repeats like Alu elements (~100 nt)
- Use **>100** to filter out short repeats

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta sln=100
```

### `-lib=` — External Repeat Library

Annotate repeats using external databases. Libraries are not included with TotalRepeats but can be obtained from:

- **Dfam/FamDB:** https://www.dfam.org/releases/current/families/FamDB
- **Repbase:** https://www.girinst.org/

You can also create custom libraries for specific species.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.ref
```

### `-combine` — Comparative Analysis

Performs synchronized repeat clustering across multiple sequences for comparative genomics. Ideal for:

- Comparing homologous chromosomes from different species or strains
- Detecting repeat polymorphisms without full-genome alignment
- Analyzing long-read assemblies (e.g., ONT sequencing)

Place all target sequences in one directory:

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -combine
```

### `-homology` — Homology Masking

Masks all homologous regions (both within and between sequences) to highlight unique sequences. Useful for:

- Identifying unique regions between closely related genomes
- Comparing human chromosomes X and Y
- Discovering sequence heterogeneity and polymorphisms

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -homology
```

### `-combinemask` — Pangenome Analysis

Performs pangenome comparative analyses using multiple masking files. Supports:

- Synchronized repeat clustering and annotation
- Cross-tool benchmarking
- Comparing outputs from different assemblies or algorithms

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -combinemask
```

### `imgx=` — Image Scaling

Controls the width of output images (PNG and SVG formats).

- **1** = Maximum compression
- **5** = Default
- **20** = Maximum stretch for detailed analysis

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta imgx=20
```

### `-readmask` — Import Mask Files

Imports existing masking files for repeat clustering, annotation, and visualization.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar masked_genome.fasta -readmask
```

### `-readgff` — Import GFF Files

Imports GFF annotation files for direct visualization.

```bash
java -Xms16g -jar TotalRepeats.jar annotation.gff -readgff
```

### `-extract` — Split FASTA

Divides a multi-entry FASTA file into separate files (one per entry).

```bash
java -jar TotalRepeats.jar multi_entry.fasta -extract
```

### `-maskscomp` — Compare Mask Files

Compares masked files from different software tools. Requires a tab-delimited file listing paired inputs:

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

### Input: FASTA Format

TotalRepeats accepts sequences in FASTA format or plain text. Multi-entry FASTA files (up to 2 GB) are supported.

```
>sequence_id Optional description
ATCGATCGATCGATCG...
```

### Output: GFF3 Format

Results are saved in [GFF3 format](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), a nine-column, tab-delimited file:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | Seqid | Source file name |
| 2 | Type | `STR` (tandem repeat), `UCRP` (unclassified), or `CRP` (classified) |
| 3 | ClusterID | Cluster identifier |
| 4 | Start | Start position (1-based) |
| 5 | End | End position (inclusive) |
| 6 | Length | Repeat element length |
| 7 | Strand | `+` (forward) or `-` (complement) |
| 8 | Phase | Not used (`.`) |
| 9 | Attributes | Sequence (if `-seqshow` is enabled) |

---

## Example Workflows

### Basic Repeat Analysis (Windows)

```bash
# Single file analysis
java -jar C:\TotalRepeats\TotalRepeats.jar C:\Genomes\NC_014637.fasta

# Custom parameters
java -jar C:\TotalRepeats\TotalRepeats.jar C:\Genomes\ kmer=16 sln=30

# Extract repeat sequences with flanking regions
java -jar C:\TotalRepeats\TotalRepeats.jar C:\Genomes\NC_014637.fasta -seqshow flanks=100
```

### Large Genome Analysis (Linux)

```bash
# Human genome with annotation
java -Xms32g -Xmx128g -jar TotalRepeats.jar /data/genomes/T2T-CHM13v2.0/ \
    -lib=/data/libraries/humsub.ref

# Tasmanian devil genome
java -Xms32g -Xmx128g -jar TotalRepeats.jar /data/genomes/Sarcophilus_harrisii/

# Very large genome (Iberian ribbed newt)
java -Xms64g -Xmx256g -jar TotalRepeats.jar /data/genomes/Pleurodeles_waltl/ -imgx=20
```

### Comparative Genomics

```bash
# Compare bacterial strains
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -combine

# Identify unique regions between genomes
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Pyricularia_oryzae/ -homology
```

---

## Author & Contact

**Ruslan Kalendar**  
📧 ruslan.kalendar@helsinki.fi

**Online Version:** https://primerdigital.com/tools/repeats.html

---

## Citation

If you use TotalRepeats in your research, please cite:

> *[Citation information to be added]*

---

## License

*[License information to be added]*

---

## Related Tools

- [Repeater](https://github.com/rkalendar/Repeater) — Pairwise sequence alignment tool
