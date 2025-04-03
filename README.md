## TotalRepeats
## Genome-wide tool for rapid de novo identification, classification, comparative analysis, and visualization of repeats

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

## TotalRepeats online: 
https://primerdigital.com/tools/repeats.html

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 23 or higher

Java Downloads: https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable: https://www.java.com/en/download/help/path.html

To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```TotalRepeats.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar <TotalRepeatsPath>\TotalRepeats.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <TotalRepeatsPath>\TotalRepeats.jar <target_file_path> optional_commands```


### Examples:
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\4.txt -nomask -nogff nsize=1

java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ kmer=12 sln=30 

java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\2.txt -seqshow flanks=100

java -jar C:\TotalRepeats\dist\TotalRepeats.jar D:\Genomes\Hydra_vulgaris\ 

```

### Large genome usage:
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ nsize=1 kmer=21 sln=90
```

Analysing all files in the folder:

```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ 

java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\
```


**Common options:**

```
kmer=19	kmer size (default kmer=19)

sln=60	repeat block length (default sln=60)

nsize=12 Speed and sensitivity of sequence clustering: nsize=0 - ignores clustering; nsize=1 - very fast clustering without sequence chain direction detection; nsize >2 - efficient clustering.

image=	the dimensionality of the image (by default, the dimensionality of the image is automatically determined), example: image=10000x300

imgx=3	figure width compression, the minimum value of imgx=1 (maximum compression), and a value of imgx=30 for the most extended figure width

flanks=100	extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

-nomask	generate a new file with masking repeats (default performed)

-nogff generate a GFF file (default performed)

-seqshow	extract repeat sequences (default not performed)

-combine multiple sequences can be analysed as one entire sequence (default not performed)

```
## kmer=
minimum value for the repeat ‘growth’ initialisation sequence. The value can be as low as 9, but for short repeats, this value can be used, but for chromosomes a value of 19 is recommended as a minimum. But it is also possible to use this value of 18.  
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ kmer=21 
```

## sln=
The minimum length of the sequence that is next to be used in the analysis. Some repeats are about 100 nucleotides long, such as Alu, so this value can be either above 100 to ignore short repeats or equal to 60-80 to detect short repeats like Alu. 
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ sln=90 
```

## -combine
This option is employed in genome-wide comparative analyses to analyze homologous sequences or chromosomes from the same or different species, or long sequences following ONT sequencing. In the context of studying homologous sequences, this approach enables the detection of heterogeneity and polymorphism without the necessity of a full-genome alignment. In this scenario, all repetitive sequences, inclusive of exons and introns, are regarded as repeats and subsequently clustered, leading to the identification of polymorphism. Additionally, the analysis of chromosomes from the same species but different strains, such as bacteria or fungi, is possible. Furthermore, the application of this option is not limited to chromosome fragments from different genomic assemblies. It is imperative that all target sequences are collected under one directory, and the full path to this directory must be provided to the application. It is allowed to analyze file containing more than one record in FASTA format. 
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ -combine
```

## nsize=
A rather important parameter for the classification of sequences. Value nsize=0 is used to ignore classification; nsize=1 - very fast clustering without determining the direction of sequences. Values 2 and higher (up to 230) are used for classification. The higher the value, the slower the algorithm will run, but the sequences will be effectively classified. The maximum value (nsize=230) can only be used when absolute stringency in classification is required. Absolute classification of sequences in practice, for certain regions, is not reachable. Therefore, already with the parameter nsize=36, the maximum in the classification of sequences will already be reached. Therefore, the default values are 7 to 12, which in reality is ideal in terms of efficiency and speed. It is recommended to use the parameter nsize=1, for maximal fast analysis of both small and large genomes and in comparative analysis (with parameter: -combine), when it is necessary to get the overall structure of repeats in genomes as quickly as possible: 
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ nsize=1
```
In the TotalRepeats online version (https://primerdigital.com/tools/repeats.html), the default setting is to classify sequences with a parameter (nsize=1).

## imgx=
The resulting image can be stretched in width if more detailed analysis is required. This value determines the width of the image, the higher the value, the longer the width. The minimum value of imgx=1 (maximum compression), and a value of imgx=30 for the most stretched image width.
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\ imgx=30	
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
2. Repeat - SSR (perfect and imperfect microsatellite repeats, and any short or long tandem repeat); UCRP - a group of repeats that did not form clusters; CRP - a group of repeats that formed interspersed clusters.
3. ClusterID - a cluster ID.
4. Start - The initial position of the repeat element in the sequence. The first base has the number 1.
5. Stop - The ending position of the repeat element (inclusive).
6. Length - length of the repeat element.
7. Strand - Valid entries include '+' (forward direction), '-' (complement direction).
8. Phase -  is not analysed for the presence of a reading frame; therefore, this value is '.'
9. Sequence (if parameter used: seqshow=true)
