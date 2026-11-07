# TotalRepeats — general information and workflow

*Reference summary derived from the application source code (release 4.1.2), for the manuscript, users and reviewers. The companion figure is `TotalRepeats_workflow.svg`.*

## Overview

TotalRepeats is a stand-alone Java command-line application for the **de novo discovery, classification, annotation, masking, visualization and pangenome comparison of repetitive DNA**. It is one program for the whole task: a single JAR with no external dependencies, no database and no pipeline of third-party tools. Repeats are found directly from the sequence by k-mer frequency, so no repeat library is needed and the same command applies to a single locus, a chromosome, a whole genome (multi-record FASTA), a folder of genomes, or a multi-assembly pangenome comparison. A run reports the soft-masked sequence, a tabular repeat annotation, genome-scale figures and — when several genomes are analysed together — a full pangenome analysis of the repeat families.

Three properties distinguish it from library-based annotators and from alignment-based *de novo* tools:

- **Alignment-free throughout.** Detection counts k-mers; family clustering compares 4-mer composition profiles (or exact k-mer containment); library annotation matches families to library sequences by the same measures. No all-against-all alignment is ever performed, which is what makes genome-scale runs take seconds to minutes.
- **Self-calibrated.** A copy joins a family when its similarity to the family seed clears a chance floor that is measured from the run's own data, separately for each block-length class, instead of a fixed sequence-identity cutoff. The threshold therefore adapts to the length distribution and base composition of the genome being analysed.
- **Deterministic and scalable.** Multithreaded clustering gives byte-identical output whatever the thread count. Comparative modes address all inputs with 64-bit coordinates over a virtual concatenation, so the combined input is limited by RAM, not by Java's 2.1 Gb single-sequence cap.

The analysis distinguishes three repeat classes: short tandem and low-complexity repeats (**STR**), dispersed repeat blocks that do not group into a family in the analysed material (**UCRP**, unclustered repeats), and **CRP** repeat families (clusters of homologous repeat copies), each optionally named after a matching entry of a user-supplied library of known repeats or genes.

## What the program does (tasks)

The core run performs four sequential tasks on each input sequence, followed by two optional steps:

1. **Sequence parsing and normalization.** FASTA records are read and nucleotides normalized (IUPAC codes kept, U→T, I→C); case is preserved where soft-masking must be respected (`-readmask`, `-extunique`, `-combinemask`). Each record of a multi-record file is analysed on its own.
2. **STR / low-complexity detection.** A sliding 17-bp window scores local di- and tri-nucleotide diversity; runs that fall below the diversity threshold are masked, short gaps between them are bridged, and surviving blocks (≥ 50 bp) form the STR track. STR positions are excluded from the repeat coverage of the next step.
3. **De novo repeat masking.** All gap-free k-mers (default *k* = 19) are counted on the forward strand and its reverse complement in an open-addressing table of base-5 codes (`LongCountSet`), which stores the only count that matters — "seen at least twice" — in one spare bit of the key. Unique k-mers are dropped, per-base coverage by repeated k-mers is accumulated, and sufficiently covered runs are merged (gaps up to `-gap`, default 2*k*) and length-filtered (`-sln`, default 80 bp) into repeat blocks. The result is a soft-masked sequence plus block coordinates.
4. **Clustering into families.** Each repeat block becomes a raw count vector over the 238 canonical 4-mers. Blocks are compared by the cosine of those vectors, in both orientations, against a chance floor calibrated per run and per length class from windows of the input itself; a long seed is matched window by window at the candidate's own scale, so a short element that is a copy of part of a longer block is still recognised. Linked blocks form CRP families; blocks that link to nothing become UCRP. `-contain` switches to asymmetric exact-k-mer containment (k chosen from the median block length, odd, 11–23). Clustering is multithreaded and deterministic.
5. **Reference annotation (optional).** With a library (`-lib=` / `-ref=`), each family is matched to the library sequences and labelled `<n>:<name>` in the annotation table, where *n* is the 1-based index of the matching library entry.
6. **Pangenome analysis (comparative modes).** Family copies are mapped back to their genome (by default one input file = one genome; `-pangroup=record` makes every FASTA record a genome). The report classifies families as core / accessory / unique and as core / soft-core / shell / cloud, gives the family frequency spectrum, exact rarefaction curves with a Heaps' law openness estimate, Jaccard and Bray–Curtis genome distances with a UPGMA tree, presence patterns, copy-number variation and per-genome summaries.

## Input data

| Input | Flag | Description |
|---|---|---|
| Target sequence(s) | *(first argument)* | A FASTA file (single or multi-record) **or** a folder of FASTA files. Nucleotide sequences only. Each record may be up to 2 GB; the combined input of a comparative run is limited only by RAM (`-Xmx`). |
| Reference library | `-lib=` / `-ref=` | Optional FASTA of known repeats or genes used to name CRP families. |
| Pre-masked FASTA | `-readmask` | Optional: a previously soft-masked FASTA, re-used for clustering, annotation and plotting. |
| Annotation table | `-readgff` | Optional: a table written by an earlier run, re-used for plotting. |
| Output folder | `-out=` | Optional output directory (default: the input's folder). |

## Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `-kmer=` | 19 | K-mer size for repeat detection (9–21; smaller for short sequences and small genomes). |
| `-sln=` | 80 | Minimum reported repeat-block length, bp. |
| `-gap=` | 2 × kmer | Largest gap bridged when merging repeated k-mers into one block. |
| `-flanks=` | 0 | Flank length added to each repeat sequence written with `-seqshow` (alias `-flangs=`). |
| `-seqshow` | off | Write the repeat (and flank) sequence into the annotation table. |
| `-nossr` | off | Disable STR / low-complexity detection (on by default). |
| `-contain` | off | Cluster by exact-k-mer containment instead of the default 4-mer composition profile. |
| `-normal` | off | Force single-threaded clustering (multithreaded by default). |
| `-image=` / `-imgx=` | auto / 10 | Figure dimensions / width scale (1–30). |

## Expected results (outputs)

| Output | File | Contents |
|---|---|---|
| Soft-masked sequence | `*.msk` | Input sequence with repeats in lower case, 70 columns per line; header records the % repeat coverage. |
| Repeat annotation | `*.gff` | One row per repeat block: *Seqid · Repeat (STR / UCRP / CRP or `n:library name`) · ClusterID · Start · Stop · Length · Strand · Phase (empty) · [Sequence ± flanks]*. Header carries the release, run parameters and genome-level metrics. |
| Repeat-landscape figures | `*.png`, `*.svg` | Genome-scale map showing the all-repeat track and per-class / per-strand tracks (STR, UCRP, CRP forward / reverse). |
| Pangenome report | `*_pangenome.txt`, `.tsv`, `_bp.tsv`, `_distance.tsv`, `.nwk`, `.svg` | Human-readable summary; copy-number and bp matrices (families × genomes) with tier, pattern and copy-number statistics; Jaccard distance matrix; UPGMA tree in Newick; four-panel figure. |
| Genome metrics | console + file headers | % sequence covered by repeats, % covered by STR, % gap (N) content, masking and total run time. |
| Auxiliary reports | `report.txt`, `result.txt`, `*.report` | Produced by `-extract`, `-maskscomp` and `-extunique`. |

A multi-record FASTA gets one set of outputs per record, numbered by position (`<input>_N.gff`, `_N.msk`, `_N.png`, `_N.svg`).

## Repeat classes and figure colour key

- **STR** — short tandem / low-complexity repeats (dark green in figures).
- **UCRP** — unclustered repeats, i.e. dispersed repeat blocks that do not form a family in the analysed material (grey).
- **CRP** — clustered repeat families of homologous copies; drawn **blue** on the forward strand and **red** on the reverse-complement strand. With a library, a matched family carries the library name.
- **Σ all-repeat track** — the projection of all repeats onto the sequence (brown).

## Operational modes

The same engine is driven into different workflows by command-line flags:

| Mode | Flag | Purpose |
|---|---|---|
| Standard (default) | *(none)* | Mask → cluster → annotate → visualize each record of a file, or each file of a folder. |
| Mask only | `-maskonly` | Produce the soft-masked FASTA only (skip clustering, annotation, figures). |
| Pairwise masking | `-amask` | Masking by seed-and-extend pairwise comparison of repeated k-mer occurrences (small inputs). |
| Comparative | `-collate` / `-joint` / `-combinemask` | Analyse a folder of genomes jointly: per-file reports, a combined report and the pangenome report. `-collate` masks each file on its own; `-joint` uses one shared k-mer table so repeats shared between files are detected; `-combinemask` starts from existing `.msk` files. |
| Comparative homology | `-homology` | Mask every region homologous within or between the files; what remains is unique (single-copy) sequence. |
| Re-use precomputed | `-readmask` / `-readgff` | Re-cluster and redraw a soft-masked FASTA, or redraw an annotation table. |
| Split | `-extract` | Split a multi-record FASTA into one file per record. |
| Unique regions | `-extunique` | Write the unmasked stretches of a soft-masked FASTA. |
| Mask comparison | `-maskscomp` | Compare masks produced by different software or parameters. |

## Representative commands

```text
# Standard run with repeat sequences and 100 bp flanks in the table
java -Xms16g -Xmx32g -jar TotalRepeats.jar <input> -seqshow -flanks=100

# Custom k-mer and block length
java -Xms16g -Xmx32g -jar TotalRepeats.jar <input> -kmer=18 -sln=100

# A folder of genomes with a reference library and an output folder
java -Xms16g -Xmx32g -jar TotalRepeats.jar /data/genomes/ -lib=/data/known.fasta -out=/data/report/

# Pangenome comparison of a folder of assemblies
java -Xms16g -Xmx64g -jar TotalRepeats.jar /data/genomes/<species>/ -collate -out=/data/report/
```

## Implementation notes (for reviewers)

- **No reference dependency.** Repeats are identified purely from k-mer frequencies in the target itself; a library is optional and used only to name families.
- **Both strands are always considered**, during masking (forward + reverse-complement k-mer counting), clustering (forward + reverse-complement comparison; the strand of each copy is reported) and annotation.
- **Self-calibrated thresholds.** The composition measure compares against a chance floor measured per run and per length class from windows of the input; the containment measure stands a 5-σ margin above its own chance level. Neither uses a fixed identity cutoff.
- **Scalability.** K-mers are compact base-5 codes in a custom open-addressing table; coordinates are 64-bit offsets over a segmented sequence store, so a comparative run may exceed the ~2.1 Gb single-string limit. Memory is controlled through the JVM heap (`-Xmx`).
- **Determinism.** Multithreaded clustering writes only to private per-candidate slots and applies assignments in index order, so the output is byte-identical for any thread count; `-normal` selects a separate single-threaded path with the same result. Pangenome rarefaction is computed exactly, not by random permutation.
- **Provenance.** Every run prints the release identifier and writes it into the header of every annotation table.
- **Tests.** A self-contained suite (297 checks, no external dependency) covers the reader, the k-mer tables, the coordinate arithmetic, mask extraction, the pangenome statistics and the command line: `./run-tests.sh`.

---
*Source: github.com/rkalendar/TotalRepeats. This document was compiled from the application source to accompany the figure and the manuscript.*
