## TotalRepeats
## A universal genome-wide tool for rapid *de novo* identification, classification, annotation, comparative analysis, and visualization of repetitive elements

## Overview

TotalRepeats is a universal, *de novo* tool for genome-wide identification, classification, annotation, visualization, and comparison of repetitive DNA elements.
It efficiently detects a wide spectrum of repeats, including: mobile genetic elements (transposons, retrotransposons), tandem arrays (microsatellites, minisatellites, telomers, satellite DNA),
Large-scale structural variations (duplications, rearrangements). 
Unlike conventional tools, TotalRepeats does not require prior sequence knowledge. 
It supports:
Rapid identification and classification of repeats;
Comparative analysis across multiple genomes or assemblies;
Interspecific polymorphism detection without whole-genome alignment;
Annotation of repeats using external libraries (Repbase, Dfam/FamDB, or custom datasets).
The tool is particularly well-suited for comparative genomics, evolutionary biology, structural variation studies, and bioinformatics research. Masking, clustering, and annotation are fully multithreaded, scaling efficiently from personal computers to HPC/supercomputing clusters.

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

## TotalRepeats online: 
https://primerdigital.com/tools/repeats.html

## Availability and requirements:

Operating system(s): Cross-platform (Windows, Linux, macOS)

Programming language: Java 24 or higher

Java Downloads: https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable: https://www.java.com/en/download/help/path.html


## Running TotalRepeats

The main executable is ```TotalRepeats.jar``` (located in the ```dist``` directory). Copy it to any folder and run:

```java -jar <TotalRepeatsPath>/TotalRepeats.jar <input_file_or_directory> [options]```

Input can be a single sequence file (FASTA/plain text) or a directory with multiple genomes.
No extra dependencies are required.

## Memory Considerations

For small genomes/files (<100 MB): No JVM memory flags needed.
For medium genomes (>500 MB): Use ≥64 GB RAM and specify -Xms / -Xmx.
For very large genomes (>2 Gb): >128 GB RAM recommended.

### Examples

```java -Xms32g -Xmx128g -jar TotalRepeats.jar E:/Genomes/Sarcophilus_harrisii/```

### Usage Examples
### Windows
```
# Basic run on a FASTA file
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\Genomes\NC_014637.fasta  

# Custom k-mer and string length
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\Genomes\ kmer=16 sln=30  

# Extract repeats with 100-nt flanks
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\Genomes\NC_014637.fasta -seqshow -flanks=100  

# Large human genome with image compression
java -Xms16g -Xmx32g -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ -imgx=20  

# Combine multiple chromosomes for synchronized classification
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Shigella\ -combine

# This is the complete process, including de novo identification, classification and visualisation of repeats, and annotation according to the user file specified in '-lib='
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ -lib=C:\TotalRepeats\test\humsub.ref
```
### Linux
```
java -Xms32g -Xmx128g -jar /data/user/dist/TotalRepeats.jar /data/genomes/Sarcophilus_harrisii/  

java -Xms64g -Xmx256g -jar /data/user/dist/TotalRepeats.jar /data/genomes/Pleurodeles_waltl/ -imgx=20

java -jar -Xms16g -Xmx32g /data/user/dist/TotalRepeats.jar /data/genomes/Shigella/ -combine

java -jar -Xms16g -Xmx32g /data/user/dist/TotalRepeats.jar /data/genomes/T2T-CHM13v2.0/ -lib=/data/user/test/humsub.ref
```

**Java Memory Management Parameters**

**-Xms (Initial Heap Size)**:
This parameter sets the initial amount of memory allocated to the Java Virtual Machine (JVM) at startup. For example, -Xms32g allocates 32 GB of heap memory immediately. This is particularly important for analyzing large genomes (e.g., >300 Mb), as it helps avoid memory allocation delays and reduces the risk of **OutOfMemoryError** during runtime.

**-Xmx (Maximum Heap Size) – Optional**:
This parameter sets the maximum amount of heap memory the JVM can use. For instance, -Xmx64g allows the application to grow the heap up to 64 GB if needed. If -Xmx is not specified, the JVM relies on a platform-dependent default, which may be too low for memory-intensive tasks such as genome-wide analysis.

**Recommendation**:
On systems with 64–128 GB of RAM, explicitly setting -Xms32g -Xmx64g is not always necessary, as the JVM typically manages memory dynamically. However, for large datasets or if you encounter **OutOfMemoryError**, specifying these parameters in your command-line execution is strongly advised. 
You need to increase the amount of memory used by the program. Restart application with "-Xms32g -Xmx64g" flags: 

```
java -Xms32g -Xmx128g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Sarcophilus_harrisii\
```

## Common options
| Command        | Description                                                               |
| -------------- | ------------------------------------------------------------------------- |
| `kmer=18`      | k-mer size (9–21; default = 18)                                           |
| `sln=90`       | Minimum repeat block length (default = 80)                                |
| `image=W×H`    | Output image size (default = auto)                                        |
| `imgx=5`       | Image width scaling (1 = compressed, 20 = stretched; default = 5)         |
| `flanks=100`   | Extend repeats by N bases (default = 0)                                   |
| `-seqshow`     | Extract repeat sequences                                                  |
| `-maskonly`    | This option generates solely the masked output file                       |
| `-readmask`    | Import masking files for clustering/annotation/visualization              |
| `-readgff`     | This function imports one or more GFF files for direct visualization of annotated repeats|

## ⚙️ Advanced Usage
| Option           | Purpose                                                                                                                                                                       |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **-combine**     | This function conducts pangenome comparative analyses utilizing multiple files for analysis, classification of repeats, and synchronized visualization. Useful for detecting repeat polymorphisms across species or assemblies. |
| **-homology**    | This involves the comparative masking of homologous regions to emphasize unique sequences, effectively masking homologous regions between sequences to highlight distinct regions. Masks homologous regions between sequences to highlight **unique regions**. Ideal for comparing closely related chromosomes (e.g., human X vs. Y).                            |
| **combinemask**  | This function performs pangenome comparative analyses using multiple **masking files**, supporting synchronized repeat clustering, annotation, visualization, and cross-tool benchmarking by comparing outputs from different assemblies or algorithms. |
| **-amask**       | Masking is executed through pairwise sequence alignment. |
| **-extract**     | Split multi-entry FASTA into separate files                               |
| **-maskscomp**   | This function compares masked files produced by different software tools or algorithms, facilitating cross-tool benchmarking.|
| **-lib=path**    | Use external repeat libraries (e.g., Repbase, Dfam/FamDB, or a custom library) for annotation.                                                                                |

## Example Commands
```
# Comparative repeat analysis across multiple bacterial strains
java -Xms16g -Xmx32g -jar TotalRepeats.jar E:/Genomes/Pyricularia_oryzae/ -combine  

# Comparative homology masking (unique vs. shared repeats)
java -Xms16g -Xmx32g -jar TotalRepeats.jar E:/Genomes/Pyricularia_oryzae/ -homology  

# Benchmark different repeat-masking outputs
java -Xms16g -Xmx32g -jar TotalRepeats.jar E:/Genomes/ -combinemask  

# Annotate repeats using Repbase/Dfam
java -Xms16g -Xmx32g -jar TotalRepeats.jar E:/Genomes/T2T-CHM13v2.0/ -lib=E:/Lib/humsub.ref  
```

## kmer=
This is the minimum value for repeat sequence masking (9-21). It can be as low as 9 for very short repeats, such as CRISPR repeats. However, a value of 19 is recommended for eukaryotic chromosomes. It is also possible to use values of 18 for short chromosomes.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ kmer=21  
```

## sln=
The minimum length of the sequence that is next to be used in the analysis. Some repeats are about 100 nucleotides long, such as Alu, so this value can be either above 100 to ignore short repeats or equal to 60-80 to detect short repeats like Alu. 
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ sln=100 
```

## lib=reference_file_path 
The TE Libraries are not included with the software; however, they can be used with user-supplied libraries by selecting the '-lib=' option. TE libraries in FamDB can be downloaded from Dfam at https://www.dfam.org/releases/current/families/FamDB.
The latest Repbase library can also be obtained at [https://www.girinst.org/](https://www.girinst.org/). 
If the program is directed to the FASTA file containing the database of existing repeats, each cluster will be annotated in accordance with this database.
The user can compile a local database containing specific elements for a specific species, or use files from the Repbase database:

```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ -lib=C:\TotalRepeats\test\humsub.ref
```

## combine
This involves the comparative masking of homologous regions to emphasize unique sequences, effectively masking homologous regions between sequences to highlight distinct regions. 
This option is employed in genome-wide comparative analyses to analyze homologous sequences or chromosomes from the same or different species, or long sequences following ONT sequencing. In the context of studying homologous sequences, this approach enables the detection of heterogeneity and polymorphism without the necessity of a full-genome alignment. In this scenario, all repetitive sequences are detected individually for each sequence. However, repeat clustering is performed for all sequences simultaneously, resulting in the identification of repeat polymorphism. In addition, this analysis is recommended for chromosomes of the same (or different) species but different strains, e.g., bacteria or fungi, or for chromosomes from different assemblies of eukaryotes or prokaryotes. Furthermore, the application of this option is not limited to chromosome fragments from different genomic assemblies. The range of application of this analysis is not limited. It is necessary that all target sequences are collected in one directory, and the full path to this directory must be provided to the application. It is acceptable to analyse a file containing more than one FASTA format record. 
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Pyricularia_oryzae\ -combine
```

## homology
Comparative Homology Masking: This option performs a comparative analysis of homologous regions and generates an individual mask for each sequence. All homologous regions both within a sequence and between different sequences are treated as repeats. It is particularly useful for analyzing homologous sequences or chromosomes from the same or different species, as well as long-read assemblies (e.g., ONT sequencing). By focusing on homologous segments, this method enables the detection of unique regions, sequence heterogeneity, and polymorphisms without requiring a full-genome alignment. Importantly, all repetitive elements including exons and introns are classified as repeats, which facilitates polymorphism discovery.  It is best applied to chromosomes from different strains (e.g., bacteria or fungi), or to assemblies from different versions of eukaryotic or prokaryotic genomes. For example, it can be used to identify unique regions distinguishing human chromosomes X and Y: shared regions are masked, while unique sequences remain in uppercase for clarity. This option is not restricted to chromosomal fragments; it can be applied to any collection of target sequences. All sequences must be placed in a single directory, and the full directory path must be provided. Multi-entry FASTA files are fully supported.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Pyricularia_oryzae\ -homology
```

## combinemask
This function performs pangenome comparative analyses using multiple masking files, supporting synchronized repeat clustering, annotation, visualization, and cross-tool benchmarking by comparing outputs from different assemblies or algorithms.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Pyricularia_oryzae\ -combinemask
```

## imgx=
The resulting image can be stretched in width if more detailed analysis is required (default imgx=5). This value determines the width of the image, the higher the value, the longer the width. The minimum value of imgx=1 (maximum compression), and a value of imgx=20 for the most stretched image width. Two files with the *.png and *.svg extensions are used to save graphic results.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ imgx=1	
```

## readmask             
This function imports one or more masking files for use in repeat clustering, annotation, and visualization. The file may contain a single FASTA entry or multiple entries.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\NC_134482.1.fasta -readmask

java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ -readmask
```

## amask           
Masking is executed through pairwise sequence alignment (https://github.com/rkalendar/Repeater).  
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ -amask
```

## readgff             
This function imports one or more GFF files for direct visualization of annotated repeats.
```
java -jar -Xms16g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\NC_134482.1.fasta.gff -readgff

java -jar -Xms16g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\ -readgff
```

## extract
This function divides a multi-entry FASTA file into separate files, producing one file per entry.
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\7\GCF_000002495.2_MG8_genomic.fna -extract
```

## maskscomp
This function compares masked files produced by different software tools or algorithms, facilitating cross-tool benchmarking.
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Test\runfiles.txt -maskscomp
```

Specify the file in which the related files for analysis are listed in two columns, separated by tabs. Example file: runfiles.txt
```
NC_017844.1.fasta	NC_017844.1.fasta.msk
NC_017849.1.fasta	NC_017849.1.fasta.msk
NC_017850.1.fasta	NC_017850.1.fasta.msk
NC_017851.1.fasta	NC_017851.1.fasta.msk
NC_017852.1.fasta	NC_017852.1.fasta.msk
NC_017853.1.fasta	NC_017853.1.fasta.msk
NC_017854.1.fasta	NC_017854.1.fasta.msk
```
 
## Sequence Entry:
Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.


## FASTA format description:
A sequence in FASTA format consists of the following:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence. 
The software supports large FASTA files (2 GB).
 

## The output is saved in a GFF3 file: a nine-column, tab-delimited, plain text file. 
 
GFF format General Feature Format describes genes and other features associated with DNA, RNA, and Protein sequences. GFF lines have nine tab-separated fields:
Generic Feature Format Version 3 (GFF3) 
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
1. Seqid - The file name.
2. Repeat - Short tandem repeat (STR) sequence (perfect and imperfect microsatellite repeats, and any short or long tandem repeat); UCRP - sequences with no identified classification; CRP - sequences with a specific classification.
3. ClusterID - a cluster ID.
4. Start - The initial position of the repeat element in the sequence. The first base has the number 1.
5. Stop - The ending position of the repeat element (inclusive).
6. Length - length of the repeat element.
7. Strand - Valid entries include '+' (forward direction), '-' (complement direction).
8. Phase -  is not analysed for the presence of a reading frame; therefore, this value is '.'
9. Sequence (if parameter used: seqshow=true)
