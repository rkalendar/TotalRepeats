# TotalRepeats

**an integrated genome-wide platform for rapid de novo identification, classification, annotation, comparative analysis, and visualization of repetitive elements**

[![Java](https://img.shields.io/badge/Java-26+-orange.svg?logo=openjdk)](https://www.oracle.com/java/technologies/downloads/)
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
  - [Clustering measures](#clustering-measures)
  - [Benchmark: sensitivity and speed on simulated genomes](#benchmark-sensitivity-and-speed-on-simulated-genomes)
  - [Choosing a measure](#choosing-a-measure)
  - [Comparing the comparative modes](#comparing-the-comparative-modes)
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

TotalRepeats is a *de novo* tool for genome-wide identification, classification, annotation, visualization, and comparison of repetitive DNA elements. It detects a broad spectrum of repeat types, including:

- **Mobile genetic elements** — retrotransposons (LINEs, SINEs, LTR elements) and DNA transposons
- **Tandem arrays** — microsatellites (SSRs), minisatellites, telomeric repeats, centromeric satellite DNA
- **Gene duplication** — duplicated genes and copy-number variants
- **Large-scale structural variations** — segmental duplications, rearrangements

The tool is intended for comparative genomics, evolutionary biology, structural-variation studies, and general bioinformatics analysis. Masking, clustering, and annotation are **multithreaded** and run across the available CPU cores on workstation and HPC hardware alike; measured runtimes are given under [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes).

In addition, the `-homology` mode inverts the search to report **unique (single-copy) sequence** — regions occurring only once across the whole genome or within an individual chromosome (see [`-homology`](#-homology--homology-masking)).

---

## What's New

Latest additions to TotalRepeats:

- **More robust command-line parsing.** Options are now parsed token by token and matched exactly, so an input path — or a `-lib=`/`-ref=`/`-out=` value — that happens to contain a mode name (e.g. a folder named `joint_genomes`) no longer accidentally switches on that mode. Flag and option names are case-insensitive and accept an optional leading dash (`joint`, `-joint`, and `--joint` are equivalent), and an unrecognised option now prints a warning instead of being silently ignored.
- **The 4-mer composition profile (`-vector`) is now the default clustering measure.** It needs almost no additional memory (no per-block k-mer index) and, unlike `-contain`, can group copies that have diverged past sharing any exact k-mer; it is slower than the previous implementation of the measure (see [Cost](#comparing-the-two-measures)). `-vector` is still accepted (it now selects the default); pass `-contain` to cluster by asymmetric k-mer containment instead, which places a short element with the longer repeat that contains it, at a higher memory cost. See [Clustering measures](#clustering-measures).

- **`-vector` rebuilt: ~5× more verified homology, with most spurious groupings eliminated.** The composition measure used one similarity threshold at every block length, even though the score chance alone produces varies enormously with length — and with the genome. That was wrong in both directions at once: checked against evidence independent of it (exact 25-mer sharing, adjudicated by Smith-Waterman against a shuffled null), it missed most real homology **and only 52% of what it grouped was real**. It now compares raw 4-mer counts by cosine against a chance floor calibrated per length scale from the run's own data, and matches a long seed window-by-window at the candidate's own scale. On `MF782455` (1.4 Mb; 765 non-SSR blocks, median block length 190 bp; default parameters) this raises clustered blocks from 205 to 436, and members sharing exact 25-mer evidence with their cluster seed from 67 to 354; in a random sample of 50 seed–member pairs adjudicated by Smith-Waterman, the fraction confirmed homologous rose from 52% (26/50) to 98% (49/50). Memory is unchanged. It can now also find copies too diverged to share any exact k-mer, which `-contain` cannot. See [Clustering measures](#clustering-measures).
- **Comparative analyses are no longer capped at ~2.1 GB.** The `-collate`, `-joint`, and `-combinemask` modes previously failed with `OutOfMemoryError: Requested string length exceeds VM limit` once the *combined* length of all input sequences passed Java's 2³¹−1 byte (≈2.1 GB) hard limit on a single array or `String`. They now concatenate the inputs **virtually** and address them with **64-bit (long) coordinates**, so a joint analysis scales to pangenome-sized data (tens of Gb) without ever materializing a single oversized string. Each individual sequence/chromosome remains capped at ~2.1 GB; the *total* across files is bounded only by available RAM.
- **Per-file reports in every comparative mode.** `-collate`, `-joint`, and `-combinemask` now each emit an individual annotation table, PNG, and SVG per input file (named after that file), *in addition to* the combined report. Each genome can be inspected on its own while the joint clustering is preserved — the same family keeps the same `ClusterID`, colour, row, and reference label in both the per-file and combined views.
- **Pangenome analysis — core / accessory / unique — across all comparative modes.** Every combine run now writes a pangenome report classifying each repeat family by how many sequences it occurs in: **core** (present in all), **accessory** (present in some), and **unique** (present in one). It includes the family-frequency spectrum, per-sequence statistics, a pairwise shared-family matrix (Jaccard similarity), and a machine-readable presence/absence matrix. See [`-collate`](#-collate--comparative-analysis).
- **Combined image is now SVG.** The joint (cross-file) visualization is written as scalable SVG, which has no pixel-size ceiling and stays sharp at pangenome scale. Per-file images are still produced as both PNG and SVG.
- **Multithreaded clustering is now the default.** Clustering runs across all CPU cores out of the box; pass `-normal` to force the single-threaded path instead. The multithreaded path is **deterministic and reproducible** — repeated runs on the same input produce identical cluster assignments and `ClusterID`s, independent of how many worker threads are used.
- **Performance improvements.** Low-complexity/SSR masking and sequence clustering were optimised (k-mer encoding and reduced memory traffic in the inner loops). Cluster assignments and `ClusterID`s are unchanged relative to the previous release on the test sets in `test/`; measured runtimes are reported under [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes).

---

## Key Features

| Feature | Description |
|---|---|
| 🔍 ***De novo* detection** | k-mer-based identification and classification of repeat classes across the spectrum listed in the [Overview](#overview) — no prior repeat library required |
| 🧬 **Comparative analysis** | Analyse multiple genomes or assemblies simultaneously with synchronized repeat clustering |
| 🔀 **Polymorphism detection** | Identify repeat presence/absence polymorphisms between genomes (intra- or interspecific) without whole-genome alignment |
| 📚 **External annotation** | Annotate repeats against Repbase, Dfam/FamDB, or custom species-specific libraries |
| 📊 **Vector and raster visualization** | PNG and SVG images of repeat landscapes, annotated maps, and comparative views, with adjustable dimensions (`image=`, `imgx=`) |
| ⚡ **Multithreaded performance** | Parallel processing across all cores for masking, clustering, and annotation (multithreaded by default; `-normal` forces single-threaded); multithreaded clustering is **deterministic** — identical, reproducible results on every run |
| 🧩 **Pangenome analysis** | Core / accessory / unique repeat-family classification with a presence/absence matrix across genomes — produced by every comparative mode (`-collate`, `-joint`, `-combinemask`) |
| 🗄️ **Pangenome-scale inputs** | Combined comparative analysis is not limited by the ~2.1 GB single-sequence cap — inputs are concatenated virtually and addressed with 64-bit coordinates |

---

## How It Works

```
┌──────────────────────────────────────────────────────────────────────────┐
│                       TotalRepeats Pipeline                              │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  1. INPUT                    2. MASKING                                  │
│  ┌────────────────┐          ┌─────────────────────────────┐             │
│  │ FASTA files    │────────▶ │ k-mer decomposition         │             │
│  │ (single or dir)│          │ Identify repeated k-mers    │             │
│  └────────────────┘          │ Merge overlapping regions   │             │
│                              └─────────────┬───────────────┘             │
│                                            │                             │
│  3. CLUSTERING                             ▼                             │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │ • Group repeat copies into families by similarity        │            │
│  │ • Classify: tandem (STR), unclassified (UCRP),           │            │
│  │   classified (CRP) via external library                  │            │
│  │ • Multithreaded comparison: 4-mer composition (default), │            │
│  │   or k-mer containment with -contain                     │            │
│  └─────────────────────────┬────────────────────────────────┘            │
│                            │                                             │
│  4. ANNOTATION             ▼               5. OUTPUT                     │
│  ┌──────────────────────────┐    ┌──────────────────────────┐            │
│  │ Match clusters against   │───▶│ Repeat annotation tables │            │
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
5. **Output** — Produces repeat annotation tables, soft-masked FASTA, PNG and SVG images, and summary statistics.

---

## Requirements

| Requirement | Details |
|---|---|
| **Java** | Version **26** or higher ([download](https://www.oracle.com/java/technologies/downloads/)) |
| **OS** | Windows, Linux, or macOS — platform-independent |
| **RAM** | Depends on genome size (see [Memory Configuration](#memory-configuration)) |

> **Verify your Java version:**
> ```bash
> java -version
> ```
> If you see a version below 26, update Java or install via [Conda](#option-2-using-conda).

**Setting `JAVA_HOME` / `PATH`:** <https://www.java.com/en/download/help/path.html>

---

## Installation

### Option 1: Direct Download

```bash
# 1. Clone the repository (or download the JAR directly)
git clone https://github.com/rkalendar/TotalRepeats.git
cd TotalRepeats

# 2. Confirm Java 26+ is available
java -version

# 3. Run — no build step required
java -jar dist/TotalRepeats.jar genome.fasta -out=/path/to/results/
```

### Option 2: Using Conda

If Java 26 is not installed system-wide:

```bash
# Create a dedicated environment
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n java26 openjdk=26
conda activate java26

# Verify
java -version
```

---

## Quick Start

```bash
# Show usage / help (also printed when no input file is given)
java -jar TotalRepeats.jar --help

# Save the usage guide to a file instead
java -jar TotalRepeats.jar --help > help.txt

# Analyse a single genome
java -jar TotalRepeats.jar genome.fasta

# Analyse all genomes in a directory
java -jar TotalRepeats.jar /path/to/genomes/ -out=/path/to/report/

# Large genome with memory allocation
java -Xms32g -Xmx64g -jar TotalRepeats.jar large_genome.fasta

# Annotate repeats with an external library
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.ref

# Compare multiple genomes side by side
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -collate
```

---

## Usage

### Basic Syntax

```bash
java [JVM flags] -jar TotalRepeats.jar <input_file_or_directory> [options]
```

- **Input** can be a single FASTA file or a directory containing multiple genomes. It must be the **first** argument.
- **Options** follow the input path and may be given in any order. Each argument is parsed as a single token and matched exactly: names are case-insensitive and a leading dash is optional (`joint`, `-joint`, and `--joint` are equivalent). An unrecognised option is reported on stderr and skipped.
- **Output folder** is set with `-out=<path>` (default: the input file's folder). There is no positional output argument.
- No additional dependencies are required.
- **Need help?** Run with `-help`, `--help`, `-h`, `--h`, `-?`, `--?`, `/?`, `/h`, `/help`, `help`, `?`, `-usage`, `--usage`, or `/usage` (case-insensitive), or with no arguments at all, to print the full usage guide and exit. The usage guide is also shown if the supplied input file or folder does not exist.

### Memory Configuration

When working with large genomes, allocate additional heap memory using JVM flags:

| Largest single sequence (bases) | Recommended RAM | JVM flags |
|---|---|---|
| < 100 Mb | JVM default heap | None |
| 100–500 Mb | 16 GB | `-Xms8g -Xmx16g` |
| 500 Mb – 2 Gb | 32 GB | `-Xms16g -Xmx32g` |
| > 2 Gb | 80 GB or more (heap plus OS and JVM overhead) | `-Xms32g -Xmx64g` |

**What the flags mean:**

- **`-Xms`** (initial heap) — Memory pre-allocated at startup. Avoids costly resizing during execution.
- **`-Xmx`** (maximum heap) — Upper memory limit. Prevents `OutOfMemoryError` for large genomes.

> **Tip:** Set `-Xms` to roughly half of `-Xmx`, as in the table above, so the heap is pre-sized close to the working set and costly resizing during the run is avoided.

> **Comparative runs (`-collate`, `-joint`, `-combinemask`):** size the heap to the *combined* length of all input files, not the largest single file. The combined length is not subject to the ~2.1 GB single-sequence cap, because the inputs are concatenated virtually and addressed with 64-bit coordinates; available heap is then the main constraint, while each individual sequence remains capped.

### Common Options

| Option | Description | Default |
|---|---|---|
| `kmer=N` | K-mer size for repeat detection (range: 9–21) | `19` |
| `sln=N` | Minimum repeat block length (bp) | `80` |
| `image=WxH` | Output image dimensions (pixels) | Auto |
| `imgx=N` | Image width scaling factor (range: 1–30; 1 = maximum compression, 20–30 = stretched) | `10` |
| `flanks=N` (alias `flangs=N`) | Extend each repeat block by N bp on both sides | `0` |
| `gap=N` | Maximum gap (bp) when merging adjacent repeated k-mers into one block | `2 × kmer` |
| `-seqshow` | Include repeat sequences in the annotation table | Off |
| `-nossr` | Disable simple-sequence-repeat (SSR) detection | On (detection enabled) |
| `-maskonly` | Generate only the masked FASTA (skip clustering/annotation) | Off |
| `-normal` | Use single-threaded clustering (multithreaded is the default) | Off |
| `-vector` | Cluster by 4-mer composition profile, compared by cosine against a chance floor calibrated per length scale from the run's own data; long seeds are matched window-by-window. Requires no per-block k-mer index, and on the sequences tested groups diverged copies that share no exact k-mer. See [Clustering measures](#clustering-measures) | Default measure |
| `-contain` | Cluster by asymmetric k-mer containment: groups size-disparate homologous blocks (a short element within a longer one) and reverse-complement copies. Builds a per-block k-mer index costing ~8 bytes per distinct k-mer per block, so size `-Xmx` accordingly. No effect in `-homology`. See [Clustering measures](#clustering-measures) | Not selected |
| `-help` | Show the usage guide and exit (aliases: `--help`, `-h`, `-?`, `/?`, `/h`, `/help`, and others — see [Troubleshooting](#troubleshooting)) | — |

`-vector` and `-contain` are mutually exclusive; if both are supplied, `-contain` takes precedence. The leading dash is optional and names are case-insensitive: `-kmer=19`, `kmer=19` and `--KMER=19` are equivalent.

### Advanced Options

| Option | Description |
|---|---|
| `-collate` | Comparative analysis — synchronized repeat clustering across multiple files. Produces per-file reports, a combined report, and a pangenome report (see [Comparing the comparative modes](#comparing-the-comparative-modes)) |
| `-joint` | Comparative analysis where all target sequences are treated as a single combined search space — repeat identification is performed over the whole set, then clustered jointly. Produces per-file reports and a pangenome report (see [Comparing the comparative modes](#comparing-the-comparative-modes)) |
| `-combinemask` | Comparative analysis with mask-based detection — synchronized clustering across files for masking-mode comparison and pangenome studies. Produces per-file reports, a combined report, and a pangenome report; no per-file masked FASTA |
| `-homology` | Mask all homologous regions to highlight unique sequences between genomes |
| `-amask` | Perform masking via pairwise sequence alignment (instead of k-mer-based) |
| `-readmask` | Import existing mask files for clustering, annotation, and visualization |
| `-readgff` | Import GFF annotation files for direct visualization |
| `-extract` | Split a multi-entry FASTA into individual files (one per sequence) |
| `-maskscomp` | Compare masked outputs from different tools or pipelines |
| `-out=<path>` | Path to output folder (default: the input file's folder) |
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
| 20–59 | Captures very short repeats (CRISPR spacers, miniature inverted-repeat transposable elements (MITEs)) |
| 60–80 | Good balance for most analyses (default) |
| 100+ | Filters out short elements; focuses on longer repeats |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta sln=100
```

### Multithreaded Clustering

Clustering is **multithreaded by default** — no flag is required. The pairwise block comparison that dominates the clustering step — the 4-mer cosine under the default measure, or k-mer containment under `-contain` — is distributed across the available CPU cores, shortening wall-clock time on genomes with many repeat families (measured timings: see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)). This applies to both [clustering measures](#clustering-measures).

The multithreaded clustering is **deterministic**: for a given input and parameter set, repeated runs produce identical cluster assignments and `ClusterID`s, and the result does not depend on the number of worker threads. Determinism follows from the design — the parallel stage writes only to private per-candidate slots and the assignments are applied in index order — and was verified by comparing outputs across runs with differing worker-thread counts. The thread count affects runtime only.

```bash
# Multithreaded by default — no flag needed
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta
```

**`-normal` — single-threaded clustering.** Pass `-normal` to run the clustering on a single thread instead of in parallel. This is slower on large genomes, but is useful for debugging or when you want the single-threaded code path.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -normal
```

**Capping threads without `-normal`.** To keep multithreading but limit how many worker threads it uses — for example on a shared cluster node, or to leave cores free for other jobs — set the standard JVM property *before* `-jar`. This bounds Java's common Fork/Join pool, which the parallel clustering uses:

```bash
# Keep multithreading, but limit it to 8 worker threads
java -Djava.util.concurrent.ForkJoinPool.common.parallelism=8 -Xms16g -Xmx32g -jar \
    TotalRepeats.jar genome.fasta
```

Setting that value to `1` runs the parallel code path on a single thread; `-normal` instead selects a separate single-threaded code path. Both end up using one thread, so either works when you need to disable parallelism.

### Clustering measures

Two different measures can group repeat blocks into families, and they answer different questions. **The 4-mer composition profile (`-vector`) is the default**; `-contain` selects asymmetric k-mer containment.

| Flag | Measure | Groups |
|---|---|---|
| *(none)* or `-vector` | 4-mer composition, compared by cosine against a per-run calibrated chance floor | **Diverged** copies that no longer share exact k-mers; blocks of differing length via scale-matched windows, though strongly size-disparate containment may be missed |
| `-contain` | Asymmetric containment on exact k-mers | Blocks of very different length, including a short element inside a long one |

#### `-contain` — asymmetric k-mer containment

A block is absorbed into a (longer) seed when a large share of its k-mers *occur inside* that seed, measured on **canonical** k-mers so that a copy inserted in the reverse-complement orientation is still found and labelled on the minus strand:

```
C(A, B) = |kmers(A) ∩ kmers(B)| / |kmers(A)|   — how much of A occurs in B
```

Ported from [GeneDistance](https://github.com/rkalendar/GeneDistance). Because it asks "how much of the short block is contained in the long one?" rather than "how similar are their overall compositions?", it places a short element with the longer repeat that contains it.

- The k-mer length is chosen automatically from the median block length (an odd value in the range 11–23), independent of the `kmer=` used for repeat *detection*.
- Costs roughly 8 bytes per distinct k-mer per block, so give the JVM enough heap (`-Xmx`) for large combined runs.

#### `-vector` — the 4-mer composition profile (the default)

Each block becomes a raw 4-mer count vector, and two blocks are compared by the **cosine** of those vectors. A candidate joins a seed when the cosine clears a **chance floor calibrated from the run's own data**, separately for each length scale.

Three design changes distinguish the current measure from the previous one, each addressing a specific failure mode:

- **A chance floor, per length scale.** The score that chance alone produces depends strongly on block length, and on the genome: two *unrelated* 250 bp windows of `MF782455` reach a cosine that a single flat threshold cannot tell from real homology, while at 2 kb chance never comes close. It is also not a constant of the algorithm — at 600 bp the chance level is far higher on a compositionally homogeneous *E. coli* replicon than on `MF782455`. So the floor is **measured per run**, from windows of the input sequence screened for homology, rather than hard-coded. This is the profile's counterpart to the containment measure's chance control, which it previously lacked entirely.
- **Scale-matched windows.** When a seed is much longer than a candidate, the candidate is compared against the best **window** of the seed at its own scale, instead of against the seed's whole composition. This is what lets it see a short element that is a copy of only *part* of a longer block — real repeat homology is usually partial, and averaging it over the whole block dilutes it away.
- **Cosine over raw counts.** The previous measure compared ratios only across 4-mers that two blocks *shared*, discarding which 4-mers were **absent** — the main signal for a short block — and it ignored any 4-mer seen 3 times or fewer, which left 40% of blocks on `MF782455` (defaults, `kmer=19`, `sln=80`) with too little to compare at all.

Its principal advantage remains memory: it requires no per-block k-mer index, so it is preferable where memory is the binding constraint — and it is the only measure that resolves homology too diverged to share exact k-mers. Each individual comparison is also *cheaper* than before (cosine is a single merge pass, whereas the previous measure was quadratic in the shared support). The remaining overhead is the per-run chance-floor calibration and the window-by-window scan of a long seed; since that scan was optimised to score candidates against dense window vectors, `-vector` is the faster of the two measures on large block sets (see [Benchmark: sensitivity and speed on simulated genomes](#benchmark-sensitivity-and-speed-on-simulated-genomes)).

#### Comparing the two measures

The two measures group genuinely different blocks. `-vector` (the default) is fast and needs almost no memory; `-contain` costs ~8 bytes per distinct k-mer per block but catches size-disparate containment (a short element inside a long one) that composition can miss.

```bash
# Default clustering (4-mer composition profile) across a pangenome
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -joint

# Group size-disparate homologous repeats by k-mer containment instead
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -joint -contain
```

Measured with default parameters (`kmer=19`, `sln=80`) on `MF782455` (a 1.4 Mb viral genome; 765 non-SSR blocks, median block length 190 bp), and on the 16 plastid genomes in `test/6` run with `-collate`:

| Mode | MF782455 | 16 genomes, `-collate` |
|---|---|---|
| *(default)* `-vector` | 436 blocks / 74 clusters | 176 / 41 |
| `-contain` | 465 / 92 | 179 / 31 |

**Blocks clustered is not by itself a measure of quality** — a loose threshold inflates it with noise — so each grouping was checked against evidence independent of both measures: whether a member shares exact 25-mers with its family's seed (the probability that a given 25-mer matches at a given position under a uniform base model is 4⁻²⁵), adjudicated by Smith-Waterman local alignment (match +2, mismatch −3, gap open/extend −5/−2, both orientations) against a null built by shuffling the seed, calling a pair homologous when it scores at least twice its own shuffled null. On `MF782455`:

| Mode | blocks clustered | members with 25-mer evidence | pairs confirmed by Smith-Waterman (random sample) |
|---|---|---|---|
| `-vector`, previous release | 205 | 67 | **52%** (26/50) |
| *(default)* `-vector`, now | 436 | **354** | **98%** (49/50) |
| `-contain` | 465 | 365 | 100% (50/50) |

The last two columns measure different things: the middle column counts members sharing exact 25-mer evidence with their seed, while the right-hand column is the Smith-Waterman verdict over a random sample of seed–member pairs. The percentages are over that sample, not over the block count.

The previous `-vector` was wrong in both directions at once: it missed most real homology *and* nearly half of what it did group showed no homology at all — those pairs aligned no better than against shuffled sequence. Both faults had the same root cause, a single similarity threshold applied at every block length even though the chance level varies enormously with length. On `MF782455` the rebuilt measure recovers about five times as many verified homologous members (67 → 354) at close to `-contain`'s precision; the size of the gain is genome-dependent (see the 7.5 Mb chromosome below).

The two measures remain complementary rather than nested. `-vector` is now the *more* sensitive of the two on some data: on `NC_014649` it groups 75 blocks to `-contain`'s 48, and 40 member–seed pairs drawn at random from those 75 were all confirmed homologous by alignment (40/40), because composition still recognises copies diverged past sharing any exact k-mer — only about half of the 75 share a 25-mer with their seed. The gain is not confined to small genomes: on a 7.5 Mb chromosome of *Pyricularia oryzae* (`ch3_MG8`), `-vector` goes from 384 blocks at 50% confirmed to **600 at 88% confirmed**.

**Cost.** Calibrating the floor and scanning windows is not free, but the scan has since been optimised: scoring a candidate against **dense** seed-window vectors by direct indexing, rather than by a sorted sparse merge, made the scan 4–6× faster at byte-identical output. Timings recorded before that optimisation (`MF782455` ~2 s to ~7 s, and ~7 s to ~17 s on the 7.5 Mb chromosome) therefore overstate the present cost and are superseded; current end-to-end timings for both measures are given under [Benchmark: sensitivity and speed on simulated genomes](#benchmark-sensitivity-and-speed-on-simulated-genomes). Memory is unchanged: the calibration builds k-mer sets for a few hundred sample windows only, never a per-block index, so `-vector` retains its small, near-constant memory footprint. `-contain` remains available where containment of size-disparate blocks is specifically required.

> **Note.** The previous release quoted figures for a 9 × *E. coli* set that is not included in this repository (`-vector` 2415 / 386, `-contain` 5835 / 733). Those `-vector` numbers predate the rebuild and no longer apply; `-contain` is unaffected. Re-run that set to refresh them.

#### Benchmark: sensitivity and speed on simulated genomes

The figures above come from real genomes, where the true repeat content is not known exactly. To measure sensitivity directly, the two measures were also run on simulated genomes carrying a planted ground truth: 20–30 families of 500–2000 bp elements (see the tables below), each inserted 50–60 times at a fixed per-base substitution rate and separated by unique random spacers of 300–400 bp drawn with equal base frequencies (≈50% GC). Assembled genome sizes were 1.2–2.4 Mb. Each condition used a fixed random seed, so every run is reproducible.

Two properties make the comparison between measures controlled. First, repeat **detection is identical** for both measures — only the clustering step differs — so both are scored on precisely the same set of blocks. Second, homology is scored **by position within the element**, corrected for strand: masking fragments an element into pieces, and two non-overlapping pieces of the same element are different sequences that should not be grouped. Scoring on family membership alone would penalise both measures for a decision that is in fact correct.

Scores are pairwise: *recall* is the fraction of truly homologous block pairs placed in one cluster, and *precision* the fraction of pairs placed in one cluster that are truly homologous.

**Runtime** (seconds, wall clock, one run per condition — values are not averaged; 4 CPU cores, JDK 26, default parameters, measured on the current build after the dense-window optimisation described above. `-maskonly` isolates the shared detection stage, so the clustering columns are derived as total − `-maskonly` and are estimates):

| Condition | Blocks | `-maskonly` | `-vector`, total | `-contain`, total | Clustering, `-vector` | Clustering, `-contain` |
|---|---|---|---|---|---|---|
| 2% divergence, 800 bp | 1556 | 1.0 | 2.3 | 1.8 | 1.3 | **0.9** |
| 8% divergence, 500 bp | 3762 | 1.3 | 4.9 | 7.5 | **3.6** | 6.2 |
| 10% divergence, 2000 bp | 7680 | 1.7 | 7.8 | 22.3 | **6.2** | 20.7 |
| 10% divergence, 50% reverse-complement | 3126 | 1.0 | 3.6 | 6.8 | **2.6** | 5.8 |
| 15% divergence, 800 bp | 1379 | 1.0 | 3.1 | 2.4 | 2.1 | **1.4** |

**Sensitivity** (position-aware pairwise scores over the same block set; "clusters" is the number of clusters resolved, against the planted family count):

| Condition | Families | `-vector`: clustered / precision / recall / clusters | `-contain`: clustered / precision / recall / clusters |
|---|---|---|---|
| 2% divergence, 800 bp | 20 | 99.4% / 0.81 / 0.95 / 33 | 99.3% / 0.81 / **0.99** / 20 |
| 8% divergence, 500 bp | 30 | 59.2% / 0.79 / **0.093** / 418 | 50.6% / 0.71 / 0.044 / 439 |
| 10% divergence, 2000 bp | 20 | 61.9% / 0.68 / **0.089** / 1062 | 34.3% / 0.82 / 0.022 / 896 |
| 10% divergence, 50% reverse-complement | 20 | 76.7% / 0.64 / **0.249** / 312 | 34.6% / 0.86 / 0.023 / 376 |
| 15% divergence, 800 bp | 20 | 44.9% / 0.95 / **0.072** / 245 | 6.0% / 1.00 / 0.007 / 38 |
| 25% divergence, 800 bp | 20 | detection recovered only 18 blocks from 1000 insertions | as `-vector` |

#### Conclusions

1. **At low divergence `-contain` was the more accurate measure in these simulations.** At 2% substitution it reached 0.99 pairwise recall against `-vector`'s 0.95, and resolved 20 clusters for the 20 planted families where `-vector` split them into 33. It was also the faster of the two on small block sets. `-contain` also groups more blocks on `MF782455` (465 vs 436); as noted above, block count alone does not establish accuracy, but the simulated ground truth is consistent with a genuine difference at low divergence.

2. **`-contain` degrades sharply as divergence rises.** At 15% substitution it clusters 6.0% of blocks (recall 0.007) against `-vector`'s 44.9% (recall 0.072) — a 7.5-fold difference in blocks clustered and an order-of-magnitude difference in recall. The mechanism is the one the two measures are built on: containment requires *exact* shared k-mers, which substitutions destroy, whereas the composition cosine still resolves similarity.

3. **`-vector` is substantially more sensitive on reverse-complement copies.** With half the insertions reverse-complemented, `-vector` clustered 76.7% of blocks against `-contain`'s 34.6%, and achieved roughly ten times the pairwise recall (0.249 vs 0.023). Canonical k-mers give `-contain` formal strand-independence, but combined with the exact-match requirement this does not survive 10% divergence.

4. **`-contain` is the more conservative measure.** Where it groups at all it rarely errs (precision 0.71–1.00 across these conditions), but it groups few blocks; `-vector` spans a wider range (0.64–0.95), trading precision for a several-fold gain in recall. Both precision figures are lower bounds (see Limitations).

5. **Both over-split once blocks are fragmented**, resolving far more clusters than the 20–30 planted families and reaching pairwise recall of only 0.007–0.25 in the fragmented conditions. This appears to be a shared ceiling of the greedy seed-and-absorb scheme rather than a property of either measure, and is the clearest target for further work.

6. **Above roughly 25% divergence the choice of measure is immaterial**, because detection fails first: k-mer masking at `kmer=19` recovered 18 blocks from 1000 planted insertions. The binding constraint is then upstream, in detection, and the remedy is a smaller `kmer=` rather than a different clustering measure.

7. **`-vector` scales better with block count.** Between the 3762-block and 7680-block runs, clustering time rose 3.3-fold for `-contain` against 1.7-fold for `-vector`. On small block sets (~1500) `-contain` is faster, because `-vector` carries a fixed chance-floor calibration cost, measured at 0.5–0.8 s across the five conditions in the runtime table, that amortises only over larger runs.

#### Choosing a measure

| Situation | Recommended |
|---|---|
| Closely related copies (low divergence), small block sets | `-contain` — higher recall at equal precision, and faster |
| Diverged copies (≳8% substitution), reverse-complement copies | `-vector` — several-fold higher recall |
| Large runs with many repeat blocks | `-vector` — markedly better scaling |
| Divergence ≳25% | Neither; reduce `kmer=` so that detection recovers the copies first |

The default (`-vector`) suits typical genome-scale runs, which combine many blocks, mixed divergence and reverse-complement copies.

**Limitations.** These are simulated sequences with uniform substitution and random spacers; real repeats also carry insertions and deletions, nesting and variable copy length, none of which is modelled here. Each condition was run once, on a single 4-core machine, and timings include JVM start-up; because each condition was run once, the reported times carry no estimate of run-to-run variance. The precision estimate is conservative, since blocks covering slightly offset regions of the same element are scored as non-homologous, so the reported values are lower bounds. Absolute recall depends on how heavily masking fragments an element, and is therefore a property of the detection stage as much as of the clustering measure.

#### Notes for both measures

- Both are deterministic and reproducible, and unaffected by the worker-thread count.
- `-homology` performs no clustering, so none of these flags have any effect in that mode.
- The output shape is the same for either measure — the same `ClusterID`s, strands, colours, per-file and combined reports, and pangenome classification — so downstream reports are unaffected by the choice of measure.
- Each option is matched as a whole token. `-vector` (the default) and `-contain` are the opt-in/opt-out pair; if both are given, `-contain` wins.

### `-lib=` — External Repeat Library

Annotate detected repeat families against a curated database. TotalRepeats does **not** ship with libraries — download them separately:

| Source | URL | Notes |
|---|---|---|
| **Dfam / FamDB** | <https://www.dfam.org/releases/current/families/FamDB> | Open access; export consensus sequences in FASTA for use with `-lib=` |
| **Repbase** | <https://www.girinst.org/> | Requires registration; curated consensus sequences |
| **Custom** | — | Build your own FASTA library for a species of interest |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.ref
```

### Comparing the comparative modes

All three comparative modes cluster repeats *jointly* across the input files, and all three now produce per-file reports, a combined report, and a pangenome report. They differ mainly in how repeats are detected and in which per-file files are written:

| Mode | Repeat detection | Combined image | Per-file outputs | Pangenome | Typical use |
|---|---|---|---|---|---|
| `-collate` | k-mer masking **per sequence**, then joint clustering | SVG | Annotation table · PNG · SVG · masked FASTA | ✅ | Compare homologous chromosomes / assemblies |
| `-joint` | sequences treated as **one** combined search space for detection | SVG | Annotation table · PNG · SVG · masked FASTA | ✅ | Treat a set of sequences as a single search space |
| `-combinemask` | mask-based detection per sequence, then joint clustering | SVG | Annotation table · PNG · SVG | ✅ | Masking-mode comparison / pangenome studies |

> All three handle a *pangenome* input larger than ~2.1 GB — see [What's New](#whats-new). The combined (cross-file) image is SVG in every mode; per-file images are PNG **and** SVG.

### `-collate` — Comparative Analysis

Performs synchronized repeat clustering across multiple sequences. All files in the input directory are processed together, enabling direct cross-genome comparison of repeat families.

**Ideal for:**
- Comparing homologous chromosomes across species or strains
- Detecting repeat polymorphisms without whole-genome alignment
- Analysing long-read assemblies (e.g., ONT, PacBio HiFi)

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/sequences/ -collate
```

**Output files produced by a `-collate` run:**

| File | Description |
|---|---|
| `<file>.gff`, `<file>.png`, `<file>.svg` | Per-file annotation and images for each input sequence |
| `<file>.msk` | Per-file soft-masked FASTA |
| `<report>.gff`, `<report>.svg` | Combined annotation and image (SVG) across all sequences |
| `<report>_pangenome.txt` | Pangenome summary (see below) |
| `<report>_pangenome.tsv` | Repeat-family presence/absence matrix |

**Pangenome report.** Using the synchronized clustering, each repeat family (cluster) — `CRP`, `UCRP` and `STR` alike, irrespective of whether it was annotated against an external library — is mapped back to the sequences it occurs in and classified as:

- **core** — present in *all* analysed sequences (shared content);
- **accessory** (shell) — present in some but not all sequences;
- **unique** (cloud) — present in a single sequence (differing content).

`<report>_pangenome.txt` reports the core/accessory/unique breakdown, the family-frequency spectrum (number of families present in exactly *k* sequences), per-sequence statistics (CRP / UCRP / STR content), the genome-specific (unique) families, and a pairwise shared-family matrix with Jaccard similarity. `<report>_pangenome.tsv` is a machine-readable presence/absence matrix — one row per family, one column per sequence — suitable for downstream analysis and plotting. The `ClusterID` column matches the GFF3 `ClusterID`, so the reports cross-reference each other.

The pangenome report is generated automatically for any multi-sequence comparative run (`-collate`, `-joint`, or `-combinemask`).

### `-combinemask` — Mask-Based Comparative Analysis

Performs pangenome-scale comparative analysis with mask-based detection. Enables synchronized clustering and annotation across assemblies, and supports comparing different masking strategies or parameters.

Like `-collate`, it writes an **individual annotation table + PNG + SVG per input file**, a **combined table + SVG**, and a **pangenome report** (`_pangenome.txt` / `_pangenome.tsv`). It does not write per-file masked FASTA. See [Comparing the comparative modes](#comparing-the-comparative-modes).

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

TotalRepeats accepts standard FASTA format or plain-text nucleotide sequences. Multi-entry FASTA files are supported. Each individual sequence (chromosome / record) may be up to 2 GB; in the comparative modes (`-collate`, `-joint`, `-combinemask`) the *pangenome* size across all files is not limited by that cap.

```
>sequence_id Optional description
ATCGATCGATCGATCGATCGATCG...
```

### Output Overview

| Output File | Format | Description |
|---|---|---|
| **Repeat annotations** | Tab-delimited table (`.gff`) | Coordinates, classification, and cluster assignments for every repeat element |
| **Masked genome** | FASTA | Soft-masked sequence (repeats in lowercase) |
| **Repeat landscape** | PNG / SVG | Vector and raster visualization of repeat distribution |
| **Summary statistics** | Text | Genome-wide repeat content, class breakdown, and family counts |
| **Per-file reports** *(comparative modes)* | Annotation table / PNG / SVG | Individual annotation and images for each input file, alongside the combined output |
| **Pangenome report** *(combine modes)* | Text + TSV | Core / accessory / unique repeat-family classification and presence/absence matrix |

### Repeat Annotation Table

Repeat annotations are written as a tab-delimited table with a header line. The layout is GFF3-*like* but is **not** [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)-conformant: the columns carry TotalRepeats-specific fields rather than GFF3's `source`, `type`, `score` and `attributes`, and no `##gff-version 3` pragma is written. Convert before use with GFF3-aware software.

| Column | Field | Description |
|---|---|---|
| 1 | **Seqid** | Source file or sequence identifier |
| 2 | **Repeat** | Repeat class: `STR` (tandem), `UCRP` (unclassified interspersed), or `CRP` (classified against an external library) |
| 3 | **ClusterID** | Repeat family cluster identifier |
| 4 | **Start** | Start position (1-based) |
| 5 | **Stop** | End position (1-based, inclusive) |
| 6 | **Length** | Repeat element length (bp) |
| 7 | **Strand** | `+` (forward) or `-` (reverse complement) |
| 8 | **Phase** | Not used |
| 9 | **Sequence** | Repeat sequence; written **only** when `-seqshow` is given |

Without `-seqshow` the file has eight columns. The file extension remains `.gff` for compatibility with existing workflows.

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
# *Homo sapiens* T2T-CHM13v2.0 assembly with Repbase annotation
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/T2T-CHM13v2.0/ \
    -lib=/data/libraries/humsub.ref -out=/report/

# *Sarcophilus harrisii* (Tasmanian devil) assembly
java -Xms16g -Xmx64g -jar TotalRepeats.jar /data/genomes/Sarcophilus_harrisii/ -out=/report/

# *Pleurodeles waltl* (Iberian ribbed newt), ~20 Gb assembly — stretched image output
java -Xms16g -Xmx64g -jar TotalRepeats.jar /data/genomes/Pleurodeles_waltl/ imgx=20
```

### Comparative Genomics

```bash
# Compare bacterial strains — per-sequence detection, then joint clustering
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -collate

# Compare bacterial strains — all sequences as one combined search space, then joint clustering
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -joint

# Compare bacterial strains — mask-based detection per strain, then joint clustering (no per-file masked FASTA)
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -combinemask

# Identify unique regions between related *Pyricularia oryzae* isolates
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

Ensure Java 26 or higher is installed and active:

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

> **`Requested string length exceeds VM limit` in a comparative run:** this specific error came from the earlier comparative code (`-collate`/`-joint`/`-combinemask`) building a single string of the whole concatenation, which Java caps at ~2.1 GB. It is resolved — these modes now concatenate virtually and use 64-bit coordinates. If you still hit a plain `OutOfMemoryError` (heap), raise `-Xmx` per the table above.

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

Clustering is multithreaded by default and already uses all available CPU cores — first make sure you have **not** passed `-normal`, which forces the slower single-threaded path:

```bash
# Multithreaded by default — do NOT add -normal
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta
```

To cap the worker-thread count (e.g. on a shared node) while keeping multithreading, add `-Djava.util.concurrent.ForkJoinPool.common.parallelism=N` before `-jar`. Multithreaded clustering is deterministic, so the thread count changes only the runtime, not the output.

---

## FAQ

**Q: Does TotalRepeats require a pre-built repeat library?**
No. TotalRepeats performs *de novo* detection using k-mer analysis — no external library is needed. However, you can optionally provide a library (Repbase, Dfam, or custom) via `-lib=` to classify and annotate the detected repeat families.

**Q: What organisms does TotalRepeats support?**
Any organism with a DNA sequence. TotalRepeats operates on nucleotide sequence alone and makes no taxon-specific assumptions; it has been applied to viral, bacterial, fungal, plant, and animal genomes.

**Q: How does TotalRepeats compare to RepeatMasker?**
RepeatMasker relies on a curated library (Repbase/Dfam) and performs homology-based detection — it can only find repeats represented in its database. TotalRepeats detects repeats *de novo*, without prior annotation, although sensitivity falls as copies diverge: under simulation, detection at `kmer=19` recovered few copies beyond ~25% substitution (see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)). The approaches are complementary: *de novo* detection for discovery, library-based search for classification. The `-maskscomp` option lets you compare their outputs directly.

**Q: Can I process multiple genomes at once?**
Yes. Point TotalRepeats at a directory containing FASTA files, and it will process all of them. Add `-collate` (or `-joint` / `-combinemask`) for synchronized cross-genome repeat clustering — each of these produces a per-file report and image for every input, a combined report and image, and a pangenome report (core / accessory / unique families plus a presence/absence matrix). See [Comparing the comparative modes](#comparing-the-comparative-modes).

**Q: What is the maximum genome size TotalRepeats can handle?**
TotalRepeats has been tested on genomes exceeding 20 Gb (e.g., *Pleurodeles waltl*). A single sequence (one FASTA record) is capped at 2 GB by Java's array/`String` limit, but in the comparative modes the *combined* size across files is not — those modes concatenate virtually and use 64-bit coordinates. The practical limiting factor is available RAM, so allocate sufficient heap via `-Xms` / `-Xmx`.

**Q: My `-collate` run crashed with `Requested string length exceeds VM limit` — is that fixed?**
Yes. That error came from concatenating all inputs into one string, which Java caps at ~2.1 GB. The comparative modes (`-collate`, `-joint`, `-combinemask`) now build the concatenation virtually and address it with 64-bit coordinates, so the combined analysis runs at pangenome scale. Any remaining memory issue would be an ordinary heap `OutOfMemoryError` — raise `-Xmx`.

**Q: Can I import results from other tools?**
Yes. Use `-readmask` to import masked FASTA files from RepeatMasker or other tools, and `-readgff` to import GFF annotations for visualization.

**Q: What does the `-normal` flag do exactly?**
By default, clustering is multithreaded — the pairwise comparison step runs across all CPU cores, reducing wall-clock clustering time on genomes with many repeat families (measured timings: see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)). Passing `-normal` turns that off and runs the clustering on a single thread: slower, but useful for debugging or when you want the single-threaded code path. If you want to keep multithreading but only limit the number of worker threads, add `-Djava.util.concurrent.ForkJoinPool.common.parallelism=N` before `-jar` instead (see [Multithreaded Clustering](#multithreaded-clustering)). The multithreaded path is deterministic: repeated runs produce identical cluster assignments and `ClusterID`s regardless of thread count.

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
  title  = {{TotalRepeats: an integrated genome-wide platform for rapid de novo identification, classification, annotation, comparative analysis, and visualization of repetitive elements}},
  url    = {https://github.com/rkalendar/TotalRepeats},
  note   = {Online version: https://primerdigital.com/tools/repeats.html}
}
```

Once a peer-reviewed description of TotalRepeats has been published, cite that article in preference to the software entry above; the reference will be added to this file on publication.

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
