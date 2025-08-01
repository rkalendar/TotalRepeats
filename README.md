## TotalRepeats
## A universal genome-wide tool for rapid *de novo* identification, classification, annotation, comparative analysis, and visualization of repetitive elements

Novel and universal tool for *de novo* identification, classification, visualization, and comparison of DNA profiles of repetitive elements at the genomic scale. This tool efficiently detects a wide range of repeats, including mobile genetic elements, tandem arrays, and large-scale genomic rearrangements, without reliance on prior sequence knowledge. The software is not limited only to the purpose of identifying any repetitive sequences and their rapid classification, but also to the comparative analysis of any target sequences, including interspecific analysis for the detection of polymorphisms without using genomic alignment. The tool uses a database of known repeats to enable annotation of repetitive sequences.
TotalRepeats is ideal for applications in comparative genomics, evolutionary biology, structural variation analysis, and general bioinformatics research.

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

## TotalRepeats online: 
https://primerdigital.com/tools/repeats.html

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 24 or higher

Java Downloads: https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable: https://www.java.com/en/download/help/path.html

To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```TotalRepeats.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar <TotalRepeatsPath>\TotalRepeats.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <TotalRepeatsPath>\TotalRepeats.jar <target_file_path> optional_commands```

This command launches the TotalRepeats application, which is packaged as a standalone .jar file. The input can be either a single genomic file or a directory containing multiple sequences. Since the tool is implemented as a standard Java application, no additional software or dependencies are required. 

Masking, clustering, and annotation are all key stages of the software that use multithreading functions for parallel computing. This significantly improves the software's performance. 

When working with large chromosome files, it is necessary to use computers with at least 64 GB of free RAM and specify the maximum amount of available RAM with the following flags: -Xms and -Xmx. You do not need to use these flags for files smaller than 100 MB. However, for files larger than 500 MB, you need more than 64 GB of free RAM and must specify the maximum values of these flags. For example, use the flags -Xms16g and -Xmx64g.

### Examples:
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\NC_014637.fasta -nomask -nogff

java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ kmer=12 sln=30 

java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\NC_014637.fasta -seqshow -flanks=100

java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\

java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Shigella\ -combine -nomask -nogff 

```

### Analysing all files in the folder:

```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ 

java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Shigella\ 
```

### Analysing large files (>200 MB, minimal RAM 64 GB):

```
java -jar -Xms16g -Xmx64g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Sarcophilus_harrisii\
```

**Java Memory Management Parameters**

**-Xms (Initial Heap Size)**:
This parameter sets the initial amount of memory allocated to the Java Virtual Machine (JVM) at startup. For example, -Xms32g allocates 32 GB of heap memory immediately. This is particularly important for analyzing large genomes (e.g., >300 Mb), as it helps avoid memory allocation delays and reduces the risk of **OutOfMemoryError** during runtime.

**-Xmx (Maximum Heap Size) – Optional**:
This parameter sets the maximum amount of heap memory the JVM can use. For instance, -Xmx64g allows the application to grow the heap up to 64 GB if needed. If -Xmx is not specified, the JVM relies on a platform-dependent default, which may be too low for memory-intensive tasks such as genome-wide analysis.

**Recommendation**:
In case of **Exception in thread ‘main’ java.lang.OutOfMemoryError: Java heap space** during runtime: you need to increase the amount of memory used by the program. On systems with 64–128 GB of RAM, explicitly setting -Xms32g -Xmx64g is not always necessary, as the JVM typically manages memory dynamically. However, for large datasets or if you encounter **OutOfMemoryError**, specifying these parameters in your command-line execution is strongly advised. 
Restart application with "-Xms32g -Xmx64g" flags: 

```
java -Xms32g -Xmx64g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Sarcophilus_harrisii\
```

**Common options:**
| Command               | Description                                   |
| ----------------------| --------------------------------------------- |
| kmer=	                | kmer size: 9 to 21 (default: k-mer = 19)|
| sln=	                | repeat block length (default sln=90) |
| nsize=                | speed and sensitivity of sequence clustering: nsize=0 - ignores clustering; nsize=1 - very fast clustering without sequence chain direction detection; nsize >1 - efficient clustering (default nsize=12) |
| image=                | the dimensionality of the image (by default, the dimensionality of the image is automatically determined), example: image=10000x300 |
| imgx=                 | figure width compression (default imgx=5), the minimum value of imgx=1 (maximum compression), and a value of imgx=20 for the most extended figure width |
| flanks=               | extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0) |
| nomask                | generate a new file with masking repeats (default performed) |
| nogff                 | generate a GFF file (default performed)|
| seqshow               | extract repeat sequences (default not performed) |
| combine               | multiple sequences can be analysed as one entire sequence (default not performed)|
| ref=file_path         | uses a database of known repeats to enable annotation of repeats (default not performed)|

	
## kmer=
This is the minimum value for repeat sequence masking (9-21). It can be as low as 9 for very short repeats, such as CRISPR repeats. However, a value of 19 is recommended for eukaryotic chromosomes. It is also possible to use values of 18 for short chromosomes.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ kmer=20 
```

## sln=
The minimum length of the sequence that is next to be used in the analysis. Some repeats are about 100 nucleotides long, such as Alu, so this value can be either above 100 to ignore short repeats or equal to 60-80 to detect short repeats like Alu. 
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ sln=100 
```

## ref=reference_file_path 
If the program is directed to the FASTA file containing the database of existing repeats, each cluster will be annotated in accordance with this database.
The user can compile a local database containing specific elements for a specific species, or use files from the Repbase database:
[https://www.girinst.org/](https://www.girinst.org/)

```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ ref=C:\TotalRepeats\test\humsub.ref
```

## combine
This option is employed in genome-wide comparative analyses to analyze homologous sequences or chromosomes from the same or different species, or long sequences following ONT sequencing. In the context of studying homologous sequences, this approach enables the detection of heterogeneity and polymorphism without the necessity of a full-genome alignment. In this scenario, all repetitive sequences are detected individually for each sequence. However, repeat clustering is performed for all sequences simultaneously, resulting in the identification of repeat polymorphism. In addition, this analysis is recommended for chromosomes of the same (or different) species but different strains, e.g., bacteria or fungi, or for chromosomes from different assemblies of eukaryotes or prokaryotes. Furthermore, the application of this option is not limited to chromosome fragments from different genomic assemblies. The range of application of this analysis is not limited. It is necessary that all target sequences are collected in one directory, and the full path to this directory must be provided to the application. It is acceptable to analyse a file containing more than one FASTA format record. 
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Pyricularia_oryzae\ -combine
```

## combine2
This option is employed in genome-wide comparative analyses to analyze homologous sequences or chromosomes from the same or different species, or long sequences following ONT sequencing. In the context of studying homologous sequences, this approach enables the detection of heterogeneity and polymorphism without the necessity of a full-genome alignment. In this scenario, all repetitive sequences, inclusive of exons and introns, are regarded as repeats and subsequently clustered, leading to the identification of polymorphism. This approach is not recommended for comparative analysis of nearly identical sequences (genomes). In addition, this analysis is recommended for chromosomes of the same (or different) species but different strains, e.g., bacteria or fungi, or for chromosomes from different assemblies of eukaryotes or prokaryotes. Furthermore, the application of this option is not limited to chromosome fragments from different genomic assemblies. The range of application of this analysis is not limited. It is necessary that all target sequences are collected in one directory, and the full path to this directory must be provided to the application. It is acceptable to analyse a file containing more than one FASTA format record.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Pyricularia_oryzae\ -combine2
```

## nsize=
A rather important parameter for the classification of sequences. Value nsize=0 is used to ignore classification; nsize=1 - very fast clustering without determining the direction of sequences. Values 2 and higher (up to 230) are used for classification. The higher the value, the slower the algorithm will run, but the sequences will be effectively classified. The maximum value (nsize=230) can only be used when absolute stringency in classification is required. Absolute classification of sequences in practice, for certain sequences, is not reachable. Therefore, already with the parameter nsize=36, the maximum in the classification of sequences will already be reached. Therefore, the default values are 7 to 12, which in reality is ideal in terms of efficiency and speed. It is recommended to use the parameter nsize=1, for maximal fast analysis of both small and large genomes and in comparative analysis (with parameter: -combine), when it is necessary to get the overall structure of repeats in genomes as quickly as possible: 
```
java -jar -Xms32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Sorghum_bicolor\ nsize=1
```
In the TotalRepeats online version (https://primerdigital.com/tools/repeats.html), the default setting is to classify sequences with a parameter (nsize=1).

## imgx=
The resulting image can be stretched in width if more detailed analysis is required (default imgx=5). This value determines the width of the image, the higher the value, the longer the width. The minimum value of imgx=1 (maximum compression), and a value of imgx=20 for the most stretched image width.
```
java -jar -Xms16g -Xmx32g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ imgx=1	
```

## Sequence Entry:
Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of the following:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.


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
