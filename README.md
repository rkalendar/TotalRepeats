## TotalRepeats
## Genome-wide tool for quick de novo identification, classification, and visualisation of repeats

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

[Web](http://primerdigital.com/tools/repeats.html)

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 23 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)


How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)


To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```TotalRepeats.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar <TotalRepeatsPath>\TotalRepeats.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <TotalRepeatsPath>\TotalRepeats.jar <target_file_path> optional_commands```


### Examples:
```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\4.txt mask=false gff=false

java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ kmer=18 sln=60 image=5000x300 mask=false seqshow=true

java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\2.txt ssr=true seqshow=true flanks=100

java -jar C:\TotalRepeats\dist\TotalRepeats.jar D:\Genomes\Hydra_vulgaris\ kmer=20 sln=100 image=10000x300

```

### Large genome usage:
```
java -jar -Xms16g -Xmx128g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\Cycas_panzhihuaensis\ kmer=21 sln=90 nsize=1 imgx=10 
```

Analysing all files in the folder:

```
java -jar C:\TotalRepeats\dist\TotalRepeats.jar C:\TotalRepeats\test\ 

java -jar -Xms16g -Xmx64g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\GRCh38.p14\
```


**Common options:**

```
kmer=	kmer size (default kmer=19)

sln=	repeat block length (default sln=60)

nsize=	speed and sensitivity of sequence clustering: nsize=0 - very fast clustering without chain direction detection; nsize=1 - used when ignoring clustering; nsize=2 - complete clustering

image=	the dimensionality of the image (by default, the dimensionality of the image is automatically determined), example: image=10000x300

imgx=	figure width compression, the minimum value of imgx=1 (maximum compression), and a value of imgx=30 for the most extended figure width

flanks=	extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

mask=true/false	generate a new file with masking repeats (default mask=true)

gff=true/false	generate a GFF file (default gff=true)

seqshow=true/false	extract repeat sequences (default seqshow=false)

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
