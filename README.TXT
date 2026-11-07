TotalRepeats 4.1.2
Self-calibrated, alignment-free discovery and pangenome comparison of repeat families.
An integrated genome-wide platform for rapid de novo identification, classification,
annotation, comparative (pangenome) analysis and visualization of repetitive elements.

Author:  Ruslan Kalendar (ruslan.kalendar@helsinki.fi)
Source:  https://github.com/rkalendar/TotalRepeats
Online:  https://digitalgens.org/repeats.html
License: GNU GPL v3.0 (LICENSE.txt)

Full documentation: README.md (this folder) and TotalRepeats_Run_Guide.docx.


REQUIREMENTS
  Java 25 or later (https://www.oracle.com/java/technologies/downloads/).
  Nothing else: TotalRepeats is a single JAR with no external dependencies,
  no database and no build step. Windows, Linux and macOS.


BASIC USAGE
  java -jar TotalRepeats.jar <input_file_or_folder> [options]

  The input (a FASTA file, or a folder of FASTA files) is always the first
  argument. Options follow in any order; names are case-insensitive and the
  leading dash is optional (kmer=19, -kmer=19 and --KMER=19 are equivalent).
  Run with no arguments or with -help for the full usage guide.

  Large inputs: give the JVM enough heap, e.g.
  java -Xms16g -Xmx32g -jar TotalRepeats.jar genome.fasta


EXAMPLES
  java -jar TotalRepeats.jar genome.fasta
  java -jar TotalRepeats.jar genome.fasta -seqshow -flanks=100
  java -jar TotalRepeats.jar genome.fasta -kmer=18 -sln=100
  java -jar TotalRepeats.jar genome.fasta -lib=/path/to/library.fasta
  java -jar TotalRepeats.jar /path/to/genomes/ -collate -out=/path/to/report/
  java -jar TotalRepeats.jar /path/to/genomes/ -homology


COMMON OPTIONS                                            default
  -kmer=N        k-mer size for repeat detection (9-21)   19
  -sln=N         minimum repeat block length, bp          80
  -gap=N         largest gap bridged inside a block       2 x kmer
  -flanks=N      flank length added to -seqshow output    0
  -seqshow       write repeat sequences into the table    off
  -nossr         disable short tandem repeat detection    detection on
  -maskonly      soft-masked FASTA only                   off
  -lib=<path>    annotate families against a FASTA library (Repbase, Dfam, custom)
  -out=<path>    output folder                            folder of the input
  -image=WxH     figure size in pixels                    automatic
  -imgx=N        figure width scale (1-30)                10
  -contain       cluster by k-mer containment instead of the default
                 4-mer composition profile
  -normal        single-threaded clustering               multithreaded

COMPARATIVE MODES (the input must be a folder)
  -collate       mask each file on its own, cluster all files jointly
  -joint         one shared k-mer table over all files, then joint clustering
  -combinemask   joint clustering of previously masked (.msk) files
  -homology      mask everything homologous within or between the files;
                 what remains is unique to each sequence
  Every comparative run writes per-file reports, a combined report and a
  pangenome report (core / accessory / unique families, rarefaction and
  Heaps' law, Jaccard and Bray-Curtis distances with a UPGMA tree, presence
  patterns, copy-number variation, matrices and an SVG figure).

OTHER UTILITIES
  -readmask      re-cluster and redraw an existing soft-masked FASTA
  -readgff       redraw figures from an existing annotation table
  -extract       split a multi-record FASTA into one file per record
  -extunique     write the unmasked (unique) stretches of a masked FASTA
  -maskscomp     compare masks produced by different tools
  -version       print the release identifier


OUTPUT
  <input>.gff    tab-delimited repeat annotation table (one row per block)
  <input>.msk    soft-masked FASTA (repeats in lower case)
  <input>.png    repeat landscape, raster
  <input>.svg    repeat landscape, vector
  report_pangenome.txt / .tsv / _bp.tsv / _distance.tsv / .nwk / .svg
                 pangenome report of a comparative run
  A multi-record FASTA gets one set of outputs per record (<input>_N.*).
