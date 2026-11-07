# TotalRepeats

**an integrated genome-wide platform for rapid de novo identification, classification, annotation, comparative analysis (pangenome), and visualization of repetitive elements**

[![Java](https://img.shields.io/badge/Java-25+-orange.svg?logo=openjdk)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-blue.svg)]()
[![Online Tool](https://img.shields.io/badge/Try%20Online-TotalRepeats-green.svg)](https://digitalgens.org/repeats.html)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)

Documentation: this file · [TotalRepeats_overview.md](TotalRepeats_overview.md) (one-page summary with the workflow figure `TotalRepeats_workflow.svg`) · `TotalRepeats_Run_Guide.docx` (the same documentation as a Word document.

---

## TotalRepeats at a glance

TotalRepeats is *self-calibrated, alignment-free discovery and pangenome comparison of repeat families* in one program. One command takes a FASTA file, or a folder of genomes, through the whole repeat analysis and leaves nothing to a second tool:

- **One tool for the whole task.** Short tandem repeat detection, *de novo* repeat discovery, soft-masking, clustering of copies into families, strand assignment, optional naming against a library, publication-ready figures, and — for several genomes — a complete pangenome analysis, all from a single run. There is no pipeline of third-party programs to install, configure or keep compatible.
- **Nothing to install but Java.** A single JAR with no external dependencies, no database and no build step. The same file runs on Windows, Linux and macOS, on a laptop or an HPC node.
- **No repeat library required.** Repeats are discovered from the sequence itself by k-mer frequency, so TotalRepeats finds the repeats that are actually in the genome — including those absent from every curated database — and works for any organism. A library (Repbase, Dfam/FamDB or your own FASTA) is optional and is used only to name the families found.
- **Alignment-free, and therefore fast.** No all-against-all alignment is ever performed: detection counts k-mers in a single streaming pass over both strands, and families are built by comparing compact composition profiles. On a 4-core desktop a bacterial genome takes seconds and a 15 Mb eukaryotic chromosome about a minute end to end (see [Measured runtimes](#measured-runtimes)).
- **Self-calibrated thresholds.** A copy joins a family when its similarity to the family seed clears a chance floor that is *measured from the run's own data*, separately for each block-length class, instead of a fixed sequence-identity cutoff. The decision level therefore adapts to the length distribution and base composition of the genome being analysed. Checked against Smith–Waterman alignment on real genomes, 98–100 % of the seed–member pairs the default measure groups are confirmed homologous (see [Clustering measures](#clustering-measures)).
- **Deterministic and reproducible.** Clustering runs on all CPU cores by default, and its output is byte-identical whatever the thread count. Pangenome rarefaction curves are computed exactly, not by random permutation. Every output file carries the release identifier that produced it.
- **Pangenome scale.** The comparative modes address all input files with 64-bit coordinates over a virtual concatenation, so the combined input is limited by RAM, not by Java's 2.1 Gb single-sequence cap. The pangenome report gives core / soft-core / shell / cloud families, rarefaction with a Heaps' law openness estimate, Jaccard and Bray–Curtis genome distances with a UPGMA tree, presence patterns and copy-number variation.
- **The inverse question too.** `-homology` masks everything that is homologous within or between genomes and leaves the sequence that is *unique* to each, and `-extunique` writes those regions out — the same engine answers "what is repeated?" and "what is single-copy?".

---

## Table of Contents

- [TotalRepeats at a glance](#totalrepeats-at-a-glance)
- [Overview](#overview)
- [Key Features](#key-features)
- [How It Works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Running the test suite](#running-the-test-suite)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Basic Syntax](#basic-syntax)
  - [Memory Configuration](#memory-configuration)
  - [Common Options](#common-options)
  - [Advanced Options](#advanced-options)
- [Detailed Option Reference](#detailed-option-reference)
  - [Clustering measures](#clustering-measures)
  - [Benchmark: sensitivity and speed on simulated genomes](#benchmark-sensitivity-and-speed-on-simulated-genomes)
  - [Measured runtimes](#measured-runtimes)
  - [Choosing a measure](#choosing-a-measure)
  - [Comparing the comparative modes](#comparing-the-comparative-modes)
  - [Pangenome report](#pangenome-report)
- [Input/Output Formats](#inputoutput-formats)
- [Example Workflows](#example-workflows)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Release notes](#release-notes)
- [Author & Contact](#author--contact)
- [Citation](#citation)
- [License](#license)
- [Related Tools](#related-tools)

---

## Overview

TotalRepeats is a *de novo* tool for the genome-wide identification, classification, annotation, visualization, and pangenomic comparison of repetitive DNA elements. It detects a broad spectrum of repeat types, including:

- **Mobile genetic elements** — retrotransposons (LINEs, SINEs, LTR elements) and DNA transposons
- **Tandem arrays** — microsatellites (SSRs), minisatellites, telomeric repeats, centromeric satellite DNA
- **Gene duplication** — duplicated genes and copy-number variants
- **Large-scale structural variations** — segmental duplications, rearrangements

The tool is intended for comparative genomics, evolutionary biology, structural-variation studies, and general bioinformatics analysis. It scales from a single locus to a pangenome: the same command analyses one sequence, every record of a multi-record FASTA, every file of a folder, or a folder of assemblies compared jointly. Repeat detection is a single streaming k-mer pass over both strands; family clustering and library annotation are **multithreaded** and run across the available CPU cores on workstation and HPC hardware alike. Measured runtimes are given under [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes) and [Measured runtimes](#measured-runtimes).

In addition, the `-homology` mode inverts the search to report **unique (single-copy) sequence** — regions occurring only once across the whole genome or within an individual chromosome (see [`-homology`](#-homology--homology-masking)).

---

## Key Features

| Feature | Description |
|---|---|
| 🔍 ***De novo* detection** | k-mer-based identification of every repeat class listed in the [Overview](#overview) — no prior repeat library required, no taxon-specific assumptions |
| 🧪 **Self-calibrated clustering** | Copies are grouped into families against a chance floor measured from the run's own data per length class, not a fixed identity cutoff; 98–100 % of grouped pairs confirmed by alignment on real genomes |
| 🧬 **Comparative analysis** | Analyse multiple genomes or assemblies simultaneously with synchronized repeat clustering: the same family keeps the same `ClusterID`, colour and label in every per-file and combined report |
| 🔀 **Polymorphism detection** | Identify repeat presence/absence polymorphisms between genomes (intra- or interspecific) without whole-genome alignment |
| 🧩 **Pangenome analysis** | Core / soft-core / shell / cloud family classification, exact rarefaction curves with a Heaps' law openness estimate, Jaccard and Bray–Curtis genome distances with a UPGMA tree, presence patterns, copy-number variation, copy-number and bp matrices and a summary figure — from every comparative mode (`-collate`, `-joint`, `-combinemask`) |
| 🗄️ **Pangenome-scale inputs** | Comparative analysis is not limited by the ~2.1 GB single-sequence cap — inputs are concatenated virtually and addressed with 64-bit coordinates |
| 📚 **External annotation** | Name repeat families against Repbase, Dfam/FamDB, or a custom species-specific library (optional) |
| 🔎 **Unique sequence** | `-homology` masks all homologous sequence within and between genomes and keeps what is single-copy; `-extunique` extracts those regions |
| 📊 **Vector and raster visualization** | PNG and SVG repeat landscapes, annotated maps, and comparative views, with adjustable dimensions (`image=`, `imgx=`) |
| ⚡ **Deterministic multithreading** | Clustering and annotation run on all cores by default (`-normal` forces one thread) and give **identical, reproducible** output on every run, whatever the thread count |
| 📦 **Zero dependencies** | One JAR, Java 25+, nothing else; the test suite (297 checks) needs only a JDK |
| 🏷️ **Provenance** | The release identifier is printed by `-version`, at the start of every run and in the header of every annotation table |

---

## How It Works

```
┌──────────────────────────────────────────────────────────────────────────┐
│                       TotalRepeats Pipeline                              │
├──────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  1. INPUT                    2. MASKING                                  │
│  ┌────────────────┐          ┌─────────────────────────────┐             │
│  │ FASTA files    │────────▶ │ STR / low-complexity scan   │             │
│  │ (single or dir)│          │ k-mer counting, both strands│             │
│  └────────────────┘          │ Coverage → repeat blocks    │             │
│                              └─────────────┬───────────────┘             │
│                                            │                             │
│  3. CLUSTERING                             ▼                             │
│  ┌──────────────────────────────────────────────────────────┐            │
│  │ • Group repeat copies into families by similarity        │            │
│  │ • Classify: tandem (STR), unclustered (UCRP),            │            │
│  │   clustered family (CRP), with strand                    │            │
│  │ • Multithreaded comparison: 4-mer composition (default), │            │
│  │   or k-mer containment with -contain                     │            │
│  └─────────────────────────┬────────────────────────────────┘            │
│                            │                                             │
│  4. ANNOTATION             ▼               5. OUTPUT                     │
│  ┌──────────────────────────┐    ┌──────────────────────────┐            │
│  │ Name families against    │───▶│ Repeat annotation tables │            │
│  │ Repbase / Dfam / custom  │    │ Soft-masked FASTA        │            │
│  │ library (optional)       │    │ PNG / SVG visualizations │            │
│  └──────────────────────────┘    │ Summary statistics       │            │
│                                  │ Pangenome report         │            │
│                                  └──────────────────────────┘            │
│                                                                          │
└──────────────────────────────────────────────────────────────────────────┘
```

1. **Input** — Reads one or more FASTA files (or an entire directory of genomes). Each record of a multi-record file is analysed on its own.
2. **Masking** — Finds short tandem / low-complexity repeats with a sliding diversity score, then counts all k-mers on both strands, marks per-base coverage by k-mers that occur at least twice, and merges covered runs into repeat blocks. The result is a soft-masked sequence and the block coordinates.
3. **Clustering** — Groups repeat blocks into families by sequence similarity in both orientations and records the strand of every copy; each block ends up as a tandem repeat (`STR`), an unclustered dispersed repeat (`UCRP`) or a member of a clustered family (`CRP`).
4. **Annotation** *(optional)* — Matches each family against Repbase, Dfam/FamDB, or a custom library and labels it with the matching entry.
5. **Output** — Produces repeat annotation tables, soft-masked FASTA, PNG and SVG images, summary statistics and, for a comparative run, the pangenome report.

---

## Requirements

| Requirement | Details |
|---|---|
| **Java** | Version **25** or higher ([download](https://www.oracle.com/java/technologies/downloads/)); tested with JDK 25 and 26 |
| **OS** | Windows, Linux, or macOS — platform-independent |
| **RAM** | Depends on genome size (see [Memory Configuration](#memory-configuration)) |
| **Anything else** | Nothing: no libraries, databases, Python/R packages or external programs |

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
java -jar dist/TotalRepeats.jar genome.fasta -out=/path/to/results/
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

### Running the test suite

The repository carries a self-contained test suite covering the sequence reader, the
k-mer maps, the virtual-concatenation coordinate arithmetic, mask-block extraction,
the pangenome statistics, and command-line handling. It needs nothing but a JDK —
there is no JUnit or other external dependency to install:

```bash
./run-tests.sh          # Linux/macOS
run-tests.bat           # Windows
```

The script compiles `src/` together with `tests/` and runs `TestRunner`, which prints
each check and exits non-zero if any fails, so it can gate a build or a CI job:

```
PASSED  297 checks in 3.6 s
```

| Suite | Checks | Covers |
|---|---|---|
| `DnaTest` | 27 | complement / reverse-complement, IUPAC codes, normalisation, shared-table stability |
| `FastaReaderTest` | 28 | record parsing, IUPAC normalisation, soft-mask-preserving raw mode, text before the first header, header-less sequences (and prose that is not one), unmodifiable accessors |
| `SeqStoreTest` | 28 | segment boundaries, global-to-local coordinates, substrings spanning segments |
| `LongIntHashMapTest` | 36 | insertion, tombstones, growth, `retainAtLeast`, `mergeInto` |
| `LongCountSetTest` | 30 | saturating counts, the bit-49 packing and its k ≤ 21 bound, pruning |
| `MaskResultTest` | 14 | block extraction, length filtering, STR exclusion, letter-range boundaries |
| `PangenomeAnalysisTest` | 94 | core / accessory / unique and tier classification, exact rarefaction checked against brute-force enumeration of genome subsets, Heaps' law on open and closed pangenomes, Jaccard / Bray–Curtis, UPGMA on a hand-solved matrix, presence patterns, copy-number statistics, record-to-file grouping |
| `CliTest` | 40 | lenient integer parsing, token-based option parsing, `-pangroup=`, the directory output filter, per-record outputs of a multi-FASTA (including too-short and empty records), `-seqshow` sequences, `-extract` file naming |

A second script, `tools/golden-check.sh`, runs every file-writing mode over deterministic generated genomes and reduces the outputs to a checksum list, so that a refactor can be proven to change no output.

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

# Compare multiple genomes side by side, with the pangenome report
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -collate
```

---

## Usage

### Basic Syntax

```bash
java [JVM flags] -jar TotalRepeats.jar <input_file_or_directory> [options]
```

- **Input** can be a single FASTA file or a directory containing multiple genomes. It must be the **first** argument. The comparative modes (`-collate`, `-joint`, `-combinemask`, `-homology`) act only on a **directory**; given a single file they are silently ignored and a normal single-file analysis is performed. When the input is a directory, files that TotalRepeats itself produced are skipped rather than re-analysed, so a folder can be re-run safely after an earlier run wrote its results there. The test is mode-aware: `.msk` files are still read by `-readmask`, `-extunique` and `-combinemask`, and `.gff` files by `-readgff`; images and pangenome reports are never treated as input. The number skipped is reported at the start of the run.
- **Options** follow the input path and may be given in any order. Each argument is parsed as a single token and matched exactly: names are case-insensitive and a leading dash is optional (`joint`, `-joint`, and `--joint` are equivalent). An unrecognised option is reported on stderr and skipped. A path that happens to contain an option name (a folder called `joint_genomes`, say) is never mistaken for an option.
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
| `flanks=N` (alias `flangs=N`) | Extend each repeat sequence written with `-seqshow` by N bp on both sides (range: 0–1000; values outside it are clamped) | `0` |
| `gap=N` | Maximum gap (bp) when merging adjacent repeated k-mers into one block | `2 × kmer` |
| `-seqshow` | Include repeat sequences in the annotation table | Off |
| `-nossr` | Disable short tandem repeat (STR/SSR) detection | On (detection enabled) |
| `-maskonly` | Generate only the masked FASTA (skip clustering/annotation) | Off |
| `-normal` | Use single-threaded clustering (multithreaded is the default) | Off |
| `-vector` | Cluster by 4-mer composition profile, compared by cosine against a chance floor calibrated per length scale from the run's own data; long seeds are matched window-by-window. Requires no per-block k-mer index, and on the sequences tested groups diverged copies that share no exact k-mer. See [Clustering measures](#clustering-measures) | Default measure |
| `-contain` | Cluster by asymmetric k-mer containment: groups size-disparate homologous blocks (a short element within a longer one) and reverse-complement copies. Builds a per-block k-mer index costing ~8 bytes per distinct k-mer per block, so size `-Xmx` accordingly. No effect in `-homology`. See [Clustering measures](#clustering-measures) | Not selected |
| `-help` | Show the usage guide and exit (aliases: `--help`, `-h`, `-?`, `/?`, `/h`, `/help`, and others — see [Basic Syntax](#basic-syntax)) | — |
| `-version` | Print the release identifier (e.g. `TotalRepeats 4.1.2`) and exit (aliases: `--version`, `-v`). The same identifier is printed at the start of every run and written into the header of every annotation table | — |

`-vector` and `-contain` are mutually exclusive; if both are supplied, `-contain` takes precedence. The leading dash is optional and names are case-insensitive: `-kmer=19`, `kmer=19` and `--KMER=19` are equivalent.

### Advanced Options

| Option | Description |
|---|---|
| `-collate` 📁 | Comparative analysis — each file is masked on its own, then all files are clustered jointly. Produces per-file reports, a combined report, and a pangenome report (see [Comparing the comparative modes](#comparing-the-comparative-modes)) |
| `-joint` 📁 | Comparative analysis where all files share one k-mer table — a repeat whose copies lie in different files is detected — then clustered jointly. Produces per-file reports, a combined report and a pangenome report |
| `-combinemask` 📁 | Comparative analysis starting from previously soft-masked (`.msk`) files — synchronized clustering across files for masking-mode comparison and pangenome studies. Produces per-file reports, a combined report, and a pangenome report; no per-file masked FASTA |
| `-homology` 📁 | Mask all homologous regions within and between the files to reveal the sequence unique to each |
| `-pangroup=<file\|record>` | What one *genome* is in the pangenome report of the three comparative modes: all records of an input file together (`file`, the default — a multi-chromosome assembly stays one genome) or every FASTA record on its own (`record`). See [Pangenome report](#pangenome-report) |
| `-amask` | Mask by seed-and-extend pairwise comparison of repeated k-mer occurrences instead of k-mer coverage (slower; for small inputs) |
| `-readmask` | Import an existing soft-masked FASTA (one or more records) for clustering, annotation, and visualization without re-masking |
| `-readgff` | Import an annotation table written by TotalRepeats and redraw the figures |
| `-extract` | Split a multi-entry FASTA into individual files (one per sequence) |
| `-extunique` | Extract the unmasked (UPPERCASE) blocks of a masked sequence — the unique, non-repetitive regions left by a previous run (see [`-extunique`](#-extunique--extract-unique-regions-from-a-mask)) |
| `-maskscomp` | Compare masked outputs from different tools or pipelines |
| `-out=<path>` | Path to output folder (default: the input file's folder) |
| `-lib=<path>` (alias `-ref=<path>`) | Annotate repeat families with an external library (Repbase, Dfam, or custom FASTA) |

> 📁 **These four modes require a directory as input.** Comparative analysis is only triggered when the first argument is a folder; if it is a single file, the flag is ignored without a warning and an ordinary single-file analysis is run — even when that file holds many FASTA records. To compare the records of one multi-entry FASTA, split it first with `-extract` and point TotalRepeats at the resulting folder.

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

Clustering is **multithreaded by default** — no flag is required. The pairwise block comparison that dominates the clustering step — the 4-mer cosine under the default measure, or k-mer containment under `-contain` — and the matching of families against a library are distributed across the available CPU cores, shortening wall-clock time on genomes with many repeat families (measured timings: see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)). This applies to both [clustering measures](#clustering-measures). Repeat detection itself is a single streaming pass over the sequence and needs no threads.

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

Two different measures can group repeat blocks into families, and they answer different questions. **The 4-mer composition profile (`-vector`) is the default**; `-contain` selects asymmetric k-mer containment. Both are alignment-free, both are deterministic, and both control for chance with a floor derived from the data rather than a fixed identity cutoff.

| Flag | Measure | Groups |
|---|---|---|
| *(none)* or `-vector` | 4-mer composition, compared by cosine against a per-run calibrated chance floor | **Diverged** copies that no longer share exact k-mers; blocks of differing length via scale-matched windows, though strongly size-disparate containment may be missed |
| `-contain` | Asymmetric containment on exact k-mers | Blocks of very different length, including a short element inside a long one |

#### `-vector` — the 4-mer composition profile (the default)

Each block becomes a raw count vector over the 238 canonical 4-mers, and two blocks are compared by the **cosine** of those vectors, in both orientations. A candidate joins a seed when the cosine clears a **chance floor calibrated from the run's own data**, separately for each length scale. Three design elements make the measure self-calibrating:

- **A chance floor, per length scale.** The score that chance alone produces depends strongly on block length, and on the genome: two *unrelated* 250 bp windows of `MF782455` reach a cosine that a single flat threshold cannot tell from real homology, while at 2 kb chance never comes close. Nor is it a constant of the algorithm — at 600 bp the chance level is far higher on a compositionally homogeneous *E. coli* replicon than on `MF782455`. So the floor is **measured per run**, from windows of the input sequence screened for homology, rather than hard-coded.
- **Scale-matched windows.** When a seed is much longer than a candidate, the candidate is compared against the best **window** of the seed at its own scale, instead of against the seed's whole composition. This is what lets it see a short element that is a copy of only *part* of a longer block — real repeat homology is usually partial, and averaging it over the whole block would dilute it away.
- **Cosine over raw counts.** Comparing full count vectors keeps the information in which 4-mers are **absent** — the main signal for a short block — and gives every block, however short, something to compare.

The measure requires no per-block k-mer index, so its memory footprint is small and nearly constant, and it is the only measure that resolves homology too diverged to share exact k-mers. Each comparison is a single merge pass over two sparse vectors; the per-run calibration costs a few hundred window-vector builds, and scoring a candidate against dense seed-window vectors by direct indexing keeps the window scan cheap, so `-vector` is the faster of the two measures on large block sets (see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)).

**Validation against alignment.** Grouping more blocks is not by itself a measure of quality — a loose threshold inflates it with noise — so the families were checked against evidence independent of both measures: whether a member shares exact 25-mers with its family's seed (the probability that a given 25-mer matches at a given position under a uniform base model is 4⁻²⁵), adjudicated by Smith–Waterman local alignment (match +2, mismatch −3, gap open/extend −5/−2, both orientations) against a null built by shuffling the seed, calling a pair homologous when it scores at least twice its own shuffled null. On `MF782455` (a 1.4 Mb viral genome; 765 non-SSR blocks, median block length 190 bp; default parameters `kmer=19`, `sln=80`):

| Measure | blocks clustered | members with 25-mer evidence | pairs confirmed by Smith–Waterman (random sample of 50) |
|---|---|---|---|
| `-vector`, fixed-threshold version of release 3 | 205 | 67 | 52 % (26/50) |
| *(default)* `-vector`, self-calibrated | 436 | **354** | **98 %** (49/50) |
| `-contain` | 465 | 365 | 100 % (50/50) |

The middle column counts members sharing exact 25-mer evidence with their seed; the right-hand column is the Smith–Waterman verdict over a random sample of seed–member pairs. Against the earlier fixed-threshold version, the self-calibrated measure recovers about five times as many verified homologous members at close to `-contain`'s precision. The gain is genome-dependent: on `NC_014649` it groups 75 blocks to `-contain`'s 48, and 40 member–seed pairs drawn at random from those 75 were all confirmed homologous by alignment (40/40), although only about half of the 75 share a 25-mer with their seed — composition recognises copies diverged past sharing any exact k-mer. On a 6.6 Mb chromosome of *Pyricularia oryzae* (`ch3_MG8`) it groups 600 blocks at 88 % confirmed, against 384 at 50 % for the fixed-threshold version.

#### `-contain` — asymmetric k-mer containment

A block is absorbed into a (longer) seed when a large share of its k-mers *occur inside* that seed, measured on **canonical** k-mers so that a copy inserted in the reverse-complement orientation is still found and labelled on the minus strand:

```
C(A, B) = |kmers(A) ∩ kmers(B)| / |kmers(A)|   — how much of A occurs in B
```

Ported from [GeneDistance](https://github.com/rkalendar/GeneDistance). Because it asks "how much of the short block is contained in the long one?" rather than "how similar are their overall compositions?", it places a short element with the longer repeat that contains it.

- The k-mer length is chosen automatically from the median block length (an odd value in the range 11–23), independent of the `kmer=` used for repeat *detection*.
- The containment ratio must stand five standard deviations above its own chance level, which keeps the expected number of false positives below one across the 10⁵–10⁶ candidate tests of a large run.
- Costs roughly 8 bytes per distinct k-mer per block, so give the JVM enough heap (`-Xmx`) for large combined runs.

#### Comparing the two measures

The two measures group genuinely different blocks. `-vector` (the default) is fast and needs almost no memory; `-contain` costs ~8 bytes per distinct k-mer per block but catches size-disparate containment (a short element inside a long one) that composition can miss.

```bash
# Default clustering (4-mer composition profile) across a pangenome
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -joint

# Group size-disparate homologous repeats by k-mer containment instead
java -Xms16g -Xmx32g -jar TotalRepeats.jar /path/to/genomes/ -joint -contain
```

Measured with default parameters (`kmer=19`, `sln=80`) on `MF782455` and on 16 plastid genomes run with `-collate`:

| Mode | MF782455 | 16 plastid genomes, `-collate` |
|---|---|---|
| *(default)* `-vector` | 436 blocks / 74 clusters | 176 / 41 |
| `-contain` | 465 / 92 | 179 / 31 |

#### Benchmark: sensitivity and speed on simulated genomes

The figures above come from real genomes, where the true repeat content is not known exactly. To measure sensitivity directly, the two measures were also run on simulated genomes carrying a planted ground truth: 20–30 families of 500–2000 bp elements (see the tables below), each inserted 50–60 times at a fixed per-base substitution rate and separated by unique random spacers of 300–400 bp drawn with equal base frequencies (≈50% GC). Assembled genome sizes were 1.2–2.4 Mb. Each condition used a fixed random seed, so every run is reproducible.

Two properties make the comparison between measures controlled. First, repeat **detection is identical** for both measures — only the clustering step differs — so both are scored on precisely the same set of blocks. Second, homology is scored **by position within the element**, corrected for strand: masking fragments an element into pieces, and two non-overlapping pieces of the same element are different sequences that should not be grouped. Scoring on family membership alone would penalise both measures for a decision that is in fact correct.

Scores are pairwise: *recall* is the fraction of truly homologous block pairs placed in one cluster, and *precision* the fraction of pairs placed in one cluster that are truly homologous.

**Runtime** (seconds, wall clock, one run per condition — values are not averaged; 4 CPU cores, JDK 25, default parameters. `-maskonly` isolates the shared detection stage, so the clustering columns are derived as total − `-maskonly` and are estimates):

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

#### Measured runtimes

End-to-end wall-clock times of complete default runs — STR detection, k-mer masking, clustering, the annotation table and the PNG and SVG figures — on an ordinary 4-core desktop (Intel Core i5-4570, 3.2 GHz, a 2013 CPU; JDK 26; `-Xmx8g`; one run each, JVM start-up included). The last column isolates detection and masking with `-maskonly`. Blocks and families are counted from the annotation table (`STR` and `UCRP` rows included in the block count; families are `ClusterID` ≥ 3):

| Input | Sequence | Repeat coverage | Blocks / families | Complete run | `-maskonly` |
|---|---|---|---|---|---|
| *Bodo saltans* virus NG1 (`MF782455`) | 1.39 Mb | 23.1 % | 1808 / 76 | 6.5 s | — |
| *Shigella flexneri* STLE1 chromosome (`CP058826`) | 4.51 Mb | 4.2 % | 376 / 68 | 9.3 s | — |
| *Pyricularia oryzae* 70-15 chromosome 3 (`NC_017849`, `ch3_MG8`) | 6.61 Mb | 12.4 % | 2162 / 106 | 23.1 s | 4.2 s |
| *Nomada hirtipes* chromosome 10 (`OX637925`) | 15.36 Mb | 60.4 % | 8692 / 932 | 64.2 s | 7.3 s |
| 16 plastid genomes, `-collate` (per-file masking, joint clustering, pangenome report) | 2.64 Mb in total | — | 1096 / 41 | 22.1 s | — |
| 16 plastid genomes, `-joint` (one shared k-mer table over all files) | 2.64 Mb in total | — | 3107 / 383 | 59.0 s | — |

Detection and masking take a few seconds even for the 15 Mb chromosome; the rest of a run is clustering, whose cost grows with the number of repeat blocks rather than with sequence length (the repeat-rich insect chromosome yields 8700 blocks), plus figure rendering. `-joint` reports many more blocks than `-collate` on the same plastid set because a sequence present once in each genome is a repeat across the set, which is exactly what that mode is for. All runs completed within the 8 GB heap.

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

Annotate detected repeat families against a curated database. `-ref=<path>` is accepted as an alias; if both are given, `-lib=` wins. TotalRepeats does **not** ship with libraries — download them separately:

| Source | URL | Notes |
|---|---|---|
| **Dfam / FamDB** | <https://www.dfam.org/releases/current/families/FamDB> | Open access; export consensus sequences in FASTA for use with `-lib=` |
| **Repbase** | <https://www.girinst.org/> | Requires registration; curated consensus sequences |
| **Custom** | — | Build your own FASTA library for a species of interest |

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.ref
```

A family that matches a library entry is labelled `<n>:<name>` in the `Repeat` column of the annotation table and in the figure legend, where *n* is the 1-based position of the entry in the library and *name* its FASTA identifier. Families without a match keep the label `CRP`.

### Comparing the comparative modes

All three comparative modes cluster repeats *jointly* across the input files, and all three produce per-file reports, a combined report, and a pangenome report. All three — like `-homology` — require the input path to be a **directory**; passed a single file, they are ignored without a warning. They differ mainly in how repeats are detected and in which per-file files are written:

| Mode | Repeat detection | Combined image | Per-file outputs | Pangenome | Typical use |
|---|---|---|---|---|---|
| `-collate` | k-mer masking **per file**, then joint clustering | SVG | Annotation table · PNG · SVG · masked FASTA | ✅ | Compare homologous chromosomes / assemblies |
| `-joint` | **one shared k-mer table** over all files — a repeat whose copies lie in different files is detected — then joint clustering | SVG | Annotation table · PNG · SVG · masked FASTA | ✅ | Treat a set of sequences as a single search space |
| `-combinemask` | repeats read from existing soft-masked (`.msk`) files, then joint clustering | SVG | Annotation table · PNG · SVG | ✅ | Masking-mode comparison / pangenome studies |

> All three handle a *pangenome* input larger than ~2.1 GB — see [Release notes](#release-notes). The combined (cross-file) image is SVG in every mode; per-file images are PNG **and** SVG.

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
| `<report>_pangenome.txt` | Pangenome summary (see [Pangenome report](#pangenome-report)) |
| `<report>_pangenome.tsv` | Repeat-family copy-number matrix (families × genomes) with category, tier, presence pattern and copy-number statistics |
| `<report>_pangenome_bp.tsv` | The same matrix in covered bp |
| `<report>_pangenome_distance.tsv` | Square Jaccard-distance matrix between genomes |
| `<report>_pangenome.nwk` | UPGMA tree of the genomes (Newick, with branch lengths) |
| `<report>_pangenome.svg` | Four-panel pangenome figure |

### Pangenome report

Every multi-genome comparative run (`-collate`, `-joint`, `-combinemask`) ends with a pangenome analysis of the jointly clustered repeat families. A **genome** is, by default, one input file: all records of a multi-record FASTA — chromosomes, contigs, plasmids — are counted together, so a folder of assemblies compares assemblies, not chromosomes. Pass `-pangroup=record` to treat every FASTA record as a genome instead (a folder holding a single multi-record file falls back to this automatically, since one file is no pangenome). Each family (`CRP` cluster, with its `ClusterID` as in the `.gff`) is mapped back to the genomes it occurs in, giving a copy-number profile over the genomes; `UCRP` and `STR` content is reported per genome as aggregates. From the profiles the report derives:

- **Classification.** Each family is **core** (present in *all* genomes), **accessory** (in some) or **unique** (in exactly one). Because these three classes get coarse once the genome count grows, every family also carries the four-tier label used by bacterial pangenome tools: **core** ≥ 99 % of genomes, **soft-core** 95–99 %, **shell** 15–95 %, **cloud** < 15 %.
- **Family frequency spectrum** — the number (and bp) of families present in exactly *k* genomes: the classic U-shaped pangenome curve.
- **Rarefaction and openness.** The expected pangenome size and core size when *n* genomes are sampled, for *n* = 1 … *N*. These are computed **exactly** from the frequency spectrum — the average over all *N*! orderings of the genomes — rather than by random permutations, so the curves are reproducible to the last digit. The number of new families contributed by the *n*-th genome is fitted to Heaps' law (Tettelin *et al.* 2008), *n*(*N*) = κ·*N*^−α: **α ≤ 1 marks an open pangenome** (new families keep appearing), **α > 1 a closed one**. The fit needs at least three genomes and is reported as an estimate below five.
- **Genome distances and tree.** For every genome pair: shared and union family counts, **Jaccard** similarity over all families and over the accessory families alone (core families are shared by definition and only compress the differences), and **Bray–Curtis** dissimilarity over copy numbers, which is sensitive to family expansions that presence/absence cannot see. A **UPGMA** tree built on the Jaccard distance is written in Newick format for any tree viewer, and the square distance matrix as TSV.
- **Presence patterns.** Families are grouped by the exact set of genomes they occur in (e.g. `1101` over four genomes), with the number of families and bp per pattern — the input of an UpSet plot, and the quickest way to see which genome subsets share lineage-specific repeats.
- **Copy-number variation.** The core families whose copy number varies most between genomes (coefficient of variation, min/max, fold range), and per-genome **expansions**: families with at least twice the mean copy number of the other genomes and at least three extra copies.
- **Per-genome summary** — genome length, families present by class, bp of core / accessory / unique content, and `UCRP` / `STR` counts and bp with the percentage of the genome they cover.

`<report>_pangenome.txt` is the human-readable summary of all of the above. `<report>_pangenome.tsv` is a machine-readable matrix with one row per family (plus an `STR` and a `UCRP` aggregate row) and one column per genome holding the copy number; it keeps the column layout of earlier releases and appends `Tier`, `Pattern`, `MinCopies`, `MaxCopies`, `MeanCopies` and `CV`. `<report>_pangenome_bp.tsv` holds the same matrix in covered bp. The figure `<report>_pangenome.svg` shows the rarefaction curves with the Heaps' exponent, the frequency spectrum, the per-genome composition, and a families × genomes presence/absence map under the UPGMA tree (when there are more families than drawable rows, consecutive families are binned and cell opacity shows the fraction present). All numbers use `.` as the decimal separator regardless of the system locale.

> **What counts as a family.** A repeat is detected only where it has at least two copies in the search space — per record in `-collate` and `-combinemask`, across the whole set in `-joint`. Single-copy sequence that two genomes merely share therefore does not form a family and is not part of this pangenome; use `-homology` to mask shared sequence directly.

### `-combinemask` — Mask-Based Comparative Analysis

Performs pangenome-scale comparative analysis starting from existing soft-masked files. Enables synchronized clustering and annotation across assemblies, and supports comparing different masking strategies or parameters — including masks produced by other tools.

Like `-collate`, it writes an **individual annotation table + PNG + SVG per input file**, a **combined table + SVG**, and the full **pangenome report** (see [Pangenome report](#pangenome-report)). It does not write per-file masked FASTA. See [Comparing the comparative modes](#comparing-the-comparative-modes).

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

The output is one soft-masked FASTA per input file; the unmasked (upper-case) stretches are the unique sequence, and `-extunique` writes them out as FASTA records.

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

Imports pre-computed masking files (e.g., from RepeatMasker or a previous TotalRepeats run) for clustering, annotation, and visualization without re-running the masking step. Lower-case letters are read as repeat; a multi-record file is processed record by record.

```bash
java -Xms16g -Xmx32g -jar TotalRepeats.jar masked_genome.fasta -readmask
```

### `-readgff` — Import an Annotation Table for Visualization

Imports an annotation table written by TotalRepeats and regenerates the figures directly, bypassing repeat detection (useful after changing `image=` or `imgx=`).

```bash
java -jar TotalRepeats.jar annotation.gff -readgff imgx=20
```

### `-extract` — Split Multi-Entry FASTA

Splits a multi-entry FASTA file into separate files (one per sequence record). Useful for preparing input when you need per-chromosome processing.

```bash
java -jar TotalRepeats.jar multi_entry.fasta -extract
```

Each record is written to `<id>.fasta`, where `<id>` is the first word of its header. Characters that are not allowed in file names (such as the `|` of NCBI identifiers) become `_`, and a repeated identifier — or one that would name the input file itself — gets `_2`, `_3`, …, so no record overwrites another; the console lists every record whose file name differs from its identifier. `report.txt` gives each record's length and soft-masked share; a folder run writes one report for all its files.

### `-extunique` — Extract Unique Regions from a Mask

The inverse of the mask: takes a soft-masked FASTA (a `.msk` file from an earlier run, or any sequence where repeats are lowercase) and writes out the **UPPERCASE** stretches — the unique, non-repetitive sequence.

Neighbouring uppercase runs separated by no more than `gap=` bases are merged, and only runs longer than `sln=` are kept, so both parameters control the result:

```bash
java -jar TotalRepeats.jar genome.fasta.msk -extunique sln=100 gap=40
```

The output is a FASTA-like file named `<input>.report` (`<input>_N.report` per record for a multi-entry input). Each extracted region carries its coordinates in the header, 1-based and inclusive:

```
>12345-67890
ATCGATCGATCG...
```

Unlike `-homology`, which computes the masking itself, `-extunique` only reads an existing mask.

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

TotalRepeats accepts standard FASTA format or plain-text nucleotide sequences. A file with no header line is read as one sequence named after the file, provided at least 90% of its letters are nucleotide codes (so a stray text file in an input folder is not analysed as DNA); text before the first header of a FASTA file is ignored. IUPAC ambiguity codes are accepted, `U` is read as `T` and `I` as `C`. Multi-entry FASTA files are supported: each record is analysed on its own and gets its own outputs, numbered by its position in the file — `<input>_N.gff`, `<input>_N.msk`, `<input>_N.png`, `<input>_N.svg` for record *N* (in the comparative modes, `<file>_N_.gff` and so on). A record not longer than `-sln` is skipped with a message on the console, so its number is absent from the outputs. Each individual sequence (chromosome / record) may be up to 2 GB; in the comparative modes (`-collate`, `-joint`, `-combinemask`) the *pangenome* size across all files is not limited by that cap.

```
>sequence_id Optional description
ATCGATCGATCGATCGATCGATCG...
```

### Output Overview

| Output File | Format | Description |
|---|---|---|
| **Repeat annotations** | Tab-delimited table (`.gff`) | Coordinates, classification, strand and family assignment for every repeat block |
| **Masked genome** | FASTA (`.msk`) | Soft-masked sequence (repeats in lowercase), 70 columns per line; the header records the % repeat coverage |
| **Repeat landscape** | PNG / SVG | Vector and raster visualization of repeat distribution: the all-repeat track, the STR and UCRP tracks, and one row per family with forward copies in blue and reverse-complement copies in red |
| **Summary statistics** | Console and table headers | Genome-wide repeat content, STR content, gap (N) content, masking and total run time, release identifier |
| **Per-file reports** *(comparative modes)* | Annotation table / PNG / SVG | Individual annotation and images for each input file, alongside the combined output |
| **Pangenome report** *(comparative modes)* | Text + TSV + Newick + SVG | Core / soft-core / shell / cloud classification, rarefaction and Heaps' law openness, genome distances with a UPGMA tree, presence patterns, copy-number variation, copy-number and bp matrices, and a four-panel figure (see [Pangenome report](#pangenome-report)) |

### Repeat Annotation Table

Repeat annotations are written as a tab-delimited table with a header line. The layout is GFF3-*like* but is **not** [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)-conformant: the columns carry TotalRepeats-specific fields rather than GFF3's `source`, `type`, `score` and `attributes`, and no `##gff-version 3` pragma is written. Convert before use with GFF3-aware software. The file extension remains `.gff` for compatibility with existing workflows.

Comment lines (`#`) at the top record the release, `kmer`, the minimum block length, the sequence length, the repeat / STR / gap percentages and the run time, followed by the column header line.

| Column | Field | Description |
|---|---|---|
| 1 | **Seqid** | Sequence identifier (the first word of the FASTA header) |
| 2 | **Repeat** | Repeat class: `STR` (short tandem / low-complexity), `UCRP` (unclustered dispersed repeat), `CRP` (member of a clustered family), or `<n>:<name>` for a family matched to entry *n* of the `-lib=` library |
| 3 | **ClusterID** | Family identifier: `1` is the STR row, `2` the UCRP row, `3` and above are families |
| 4 | **Start** | Start position (1-based) |
| 5 | **Stop** | End position (1-based, inclusive) |
| 6 | **Length** | Repeat block length (bp) |
| 7 | **Strand** | `+` (forward) or `-` (reverse complement of the family seed); always `+` for STR and UCRP |
| 8 | **Phase** | Always empty (kept for the GFF-like column layout) |
| 9 | **Sequence** | Only with `-seqshow`: the block sequence from Start to Stop, reverse-complemented for `-` strand copies, with `flanks=` bases of flanking sequence in upper case on either side (flanks stop at the record boundary) |

Rows have eight tab-separated fields without `-seqshow` and nine with it; the sequence, when present, is always the last field.

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
# Compare bacterial strains — per-file detection, then joint clustering and the pangenome report
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -collate

# Compare bacterial strains — one shared k-mer table over all strains, then joint clustering
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/Shigella/ -joint

# Compare bacterial strains — start from existing masks, then joint clustering (no per-file masked FASTA)
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

The recognised help flags are listed under [Basic Syntax](#basic-syntax).

### Java version errors

Ensure Java 25 or higher is installed and active:

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

> **`Requested string length exceeds VM limit` in a comparative run:** this specific error came from the comparative code of releases before 4.0 (`-collate`/`-joint`/`-combinemask`) building a single string of the whole concatenation, which Java caps at ~2.1 GB. It is resolved — these modes now concatenate virtually and use 64-bit coordinates. If you still hit a plain `OutOfMemoryError` (heap), raise `-Xmx` per the table above.

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
No. TotalRepeats performs *de novo* detection using k-mer analysis — no external library is needed. However, you can optionally provide a library (Repbase, Dfam, or custom) via `-lib=` to name the detected repeat families.

**Q: What organisms does TotalRepeats support?**
Any organism with a DNA sequence. TotalRepeats operates on nucleotide sequence alone and makes no taxon-specific assumptions; it has been applied to viral, bacterial, fungal, plant, and animal genomes.

**Q: How does TotalRepeats compare to RepeatMasker?**
RepeatMasker relies on a curated library (Repbase/Dfam) and performs homology-based detection — it can only find repeats represented in its database. TotalRepeats detects repeats *de novo*, without prior annotation, although sensitivity falls as copies diverge: under simulation, detection at `kmer=19` recovered few copies beyond ~25% substitution (see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)). The approaches are complementary: *de novo* detection for discovery, library-based search for classification. The `-maskscomp` option lets you compare their outputs directly, and `-readmask` lets TotalRepeats cluster and draw a mask produced by RepeatMasker.

**Q: What is the difference between TotalRepeats and a de novo pipeline such as RepeatModeler?**
A *de novo* pipeline chains several programs (k-mer or alignment-based discovery, consensus building, classification, masking) each with its own installation and parameters. TotalRepeats performs discovery, STR detection, family clustering, strand assignment, masking, optional naming, figures and the pangenome comparison in one dependency-free program, without any alignment step, and its clustering threshold is calibrated from the data rather than fixed.

**Q: Can I process multiple genomes at once?**
Yes. Point TotalRepeats at a directory containing FASTA files, and it will process all of them. Add `-collate` (or `-joint` / `-combinemask`) for synchronized cross-genome repeat clustering — each of these produces a per-file report and image for every input, a combined report and image, and a pangenome report (core / accessory / unique and core / soft-core / shell / cloud families, rarefaction curves with Heaps' law, genome distances and a UPGMA tree, presence patterns, copy-number variation, matrices and a figure). Each input file counts as one genome; `-pangroup=record` compares FASTA records instead. See [Comparing the comparative modes](#comparing-the-comparative-modes) and [Pangenome report](#pangenome-report).

**Q: What is the maximum genome size TotalRepeats can handle?**
TotalRepeats has been tested on genomes exceeding 20 Gb (e.g., *Pleurodeles waltl*). A single sequence (one FASTA record) is capped at 2 GB by Java's array/`String` limit, but in the comparative modes the *combined* size across files is not — those modes concatenate virtually and use 64-bit coordinates. The practical limiting factor is available RAM, so allocate sufficient heap via `-Xms` / `-Xmx`.

**Q: Are the results reproducible?**
Yes. Given the same input and parameters, every mode writes byte-identical output on every run and on any number of threads: the parallel clustering is deterministic by construction, pangenome rarefaction is computed exactly rather than by random permutation, and no step uses a random seed. The release identifier is written into every annotation table, so the version behind any result can be recovered from the output itself.

**Q: Can I import results from other tools?**
Yes. Use `-readmask` to import masked FASTA files from RepeatMasker or other tools, and `-readgff` to redraw figures from an annotation table written by TotalRepeats.

**Q: What does the `-normal` flag do exactly?**
By default, clustering is multithreaded — the pairwise comparison step runs across all CPU cores, reducing wall-clock clustering time on genomes with many repeat families (measured timings: see [Benchmark](#benchmark-sensitivity-and-speed-on-simulated-genomes)). Passing `-normal` turns that off and runs the clustering on a single thread: slower, but useful for debugging or when you want the single-threaded code path. If you want to keep multithreading but only limit the number of worker threads, add `-Djava.util.concurrent.ForkJoinPool.common.parallelism=N` before `-jar` instead (see [Multithreaded Clustering](#multithreaded-clustering)). The multithreaded path is deterministic: repeated runs produce identical cluster assignments and `ClusterID`s regardless of thread count.

---

## Release notes

**4.1.2** — Multi-record FASTA files, record by record. Records shorter than the k-mer, empty records and all-gap records no longer crash the comparative modes or report `NaN`; text before the first header is ignored; a header-less nucleotide sequence is read as one record; `-extract` writes valid, distinct file names and never overwrites its input; `-joint` / `-homology` name their masks after the record, as `-collate` does. The annotation table written with `-seqshow` now carries the sequence in its own ninth column, after the (empty) `Phase` column, so the data rows match the header. Repeat detection, masking and clustering are unchanged from 4.1.0.

**4.1.1** — Every record of a multi-entry FASTA gets its own outputs — `<input>_N.gff`, `_N.msk`, `_N.png`, `_N.svg` — with its own name, statistics and sequences; `-seqshow` sequences run exactly from Start to Stop with flanks that stop at the record boundary.

**4.1.0** — Pangenome analysis rebuilt. The report of `-collate`, `-joint` and `-combinemask` treats each input file as one genome (`-pangroup=record` restores one genome per record) and adds the four-tier core / soft-core / shell / cloud classification, exact rarefaction curves with a Heaps' law fit, Jaccard and Bray–Curtis distances with a UPGMA tree in Newick format, presence patterns, copy-number variation and a four-panel SVG figure. New files: `_pangenome_bp.tsv`, `_pangenome_distance.tsv`, `_pangenome.nwk`, `_pangenome.svg`. Numbers are written with a `.` decimal separator whatever the system locale.

**4.0** —
- *Self-calibrated clustering.* The 4-mer composition measure (`-vector`, now the default) compares raw 4-mer counts by cosine against a chance floor calibrated per length scale from the run's own data and matches long seeds window by window; on `MF782455` this raised the share of grouped pairs confirmed by Smith–Waterman from 52 % to 98 % while grouping five times as many verified homologous members (see [Clustering measures](#clustering-measures)). `-contain` selects exact-k-mer containment instead.
- *Pangenome-scale comparative runs.* `-collate`, `-joint` and `-combinemask` concatenate their inputs virtually and address them with 64-bit coordinates, so the combined input is no longer capped at ~2.1 GB (`Requested string length exceeds VM limit`).
- *Per-file reports in every comparative mode*, as exact slices of the joint clustering, so a family keeps the same `ClusterID`, colour and label in every view; the combined image is written as SVG.
- *Multithreaded clustering by default*, deterministic and independent of the thread count; `-normal` selects the single-threaded path.
- *Token-based command-line parsing*: options are matched as whole words, case-insensitively, with an optional leading dash; a path containing an option name no longer switches that option on, and unrecognised options are reported.
- *Performance*: k-mer counting moved to a value-less open-addressing set (`LongCountSet`), and the STR detector and clustering inner loops were optimised.

---

## Author & Contact

**Ruslan Kalendar**
📧 ruslan.kalendar@helsinki.fi

🌐 **Online version:** <https://digitalgens.org/repeats.html>
💻 **Source code:** <https://github.com/rkalendar/TotalRepeats>

---

## Citation

If you use TotalRepeats in your research, please cite the software release you used. The release identifier is printed by `-version` and written into the header of every annotation table, so the version behind any result can always be recovered from the output itself. Machine-readable metadata are in [`CITATION.cff`](CITATION.cff) and [`.zenodo.json`](.zenodo.json).

> Kalendar, R. (2026). *TotalRepeats: self-calibrated, alignment-free discovery and pangenome comparison of repeat families* (version 4.1.2) [Computer software]. University of Helsinki. <https://github.com/rkalendar/TotalRepeats>. Archived at Zenodo: <https://doi.org/10.5281/zenodo.13735264>

```bibtex
@software{Kalendar_TotalRepeats_2026,
  author  = {Kalendar, Ruslan},
  title   = {{TotalRepeats: self-calibrated, alignment-free discovery and pangenome comparison of repeat families}},
  version = {4.1.2},
  year    = {2026},
  url     = {https://github.com/rkalendar/TotalRepeats},
  doi     = {10.5281/zenodo.13735264},
  note    = {Online version: https://digitalgens.org/repeats.html}
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
| **[GeneDistance](https://github.com/rkalendar/GeneDistance)** | Alignment-free sequence comparison; the source of the k-mer containment measure used by `-contain` |
| **[Repeater](https://github.com/rkalendar/Repeater)** | Pairwise sequence alignment tool |
| **[NCBI RefSeq GenBank Downloader](https://github.com/rkalendar/genbanktools)** | Batch-download GenBank records by gene or accession |

---

<p align="center">
  <em>TotalRepeats — Genome-wide repeat discovery, from laptop to supercomputer.</em>
</p>
