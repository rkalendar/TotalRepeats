To distribute this project, zip up the dist folder (including the lib folder) and distribute the ZIP file.

To run the project from the command line. Command-line options, separated by spaces. 

The executive file is in the dist directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified.

To run the project from the command line, go to the dist folder and type the following:

java -jar TotalRepeats.jar <target_file_path/Folder_path>

For large chromosomes exceeding 100 Mbp in size, the software requires a minimum of 32 GB of system memory to ensure optimal performance. 
Users must allocate the desired amount of memory by specifying the initial Java heap size using the -Xms flag. 
For example, to allocate 32 GB of memory, the following command should be used: 

java -jar -Xms32g TotalRepeats.jar <inputfile/folder>

java -jar -Xms16g C:\TotalRepeats\dist\TotalRepeats.jar E:\Genomes\T2T-CHM13v2.0\

This parameter sets the initial heap size allocated to the Java Virtual Machine (JVM). In this example, -Xms32g tells the JVM to allocate 32 gigabytes of memory at startup. This is especially important for processing large genomes, as it helps prevent memory-related errors and ensures efficient performance during analysis. 
Optional: -Xmx flag. The -Xmx parameter defines the maximum heap size that the JVM is allowed to use. For example, -Xmx64g would allow the program to use up to 64 GB of RAM if needed. If not specified, the JVM uses a platform-dependent default maximum, which may be insufficient for large-scale analyses:

java -Xms32g -Xmx64g TotalRepeats.jar <inputfile/folder>

The comparative analysis mode, activated via the -combine flag, is specifically designed for genome-wide comparison of two or more input sequences. Example usage: 

java -jar TotalRepeats.jar <inputfile/folder> -combine 

Additionally, repetitive sequences can be annotated using a user-defined reference database. Users may import established repositories such as GIRI (https://www.girinst.org/) or construct custom databases tailored to specific research needs. To enable reference-based annotation, the -ref=<file_path> flag should be specified: 

java -jar TotalRepeats.jar <inputfile/folder> -ref=humsub.ref
