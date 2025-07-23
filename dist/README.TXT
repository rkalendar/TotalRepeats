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

The comparative analysis mode, activated via the -combine flag, is specifically designed for genome-wide comparison of two or more input sequences. Example usage: 

java -jar TotalRepeats.jar <inputfile/folder> -combine 

Additionally, repetitive sequences can be annotated using a user-defined reference database. Users may import established repositories such as GIRI (https://www.girinst.org/) or construct custom databases tailored to specific research needs. To enable reference-based annotation, the -ref=<file_path> flag should be specified: 

java -jar TotalRepeats.jar <inputfile/folder> -ref=humsub.ref
