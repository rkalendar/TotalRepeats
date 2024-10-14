TotalRepeats "Genome-wide tool for rapid de novo identification, classification and visualisation of interspersed and tandem repeats"
by Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi
web: https://primerdigital.com/tools/


Availability and requirements:
Operating system(s): Platform independent
Programming language: Java
Other requirements: Java 23 or higher
Java Downloads:
https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable:
https://www.java.com/en/download/help/path.html


To run the project from the command line, go to the target folder and type the following; an individual file or a file folder can be specified:

java -jar Repeater3.jar <inputfile>/<inputfolderpath> <optional_commands> 

Basic usage:

java -jar \Repeater3\dist\Repeater3.jar <inputfilepath> <optional_commands>

Examples:

java -jar \Repeater3\dist\Repeater3.jar \Repeater3\test\ef067844.txt kmer=21 seqshow=true flanks=100

Large genome usage:

java -jar -Xms4g -Xmx32g \Repeater3\dist\Repeater3.jar <inputfolderpath> kmer=21 min=100 sln=550


Common options:
kmer=19	 minimal kmer=12 (default kmer=19)

sln=100	repeat block length (default sln=90)

image=5000x3000 (by default, the dimensionality of the image is automatically determined)

flanks=100	 extend the flanks of the repeat with an appropriate length (100 nt) (default flanks=0)

mask=true/false	 generate a new file with masking repeats (default mask=true)

seqshow=true/false	 extract repeat sequences (default seqshow=false)

ssr=true	 analyzing only the SSR/telomers loci (default ssr=false) 


 
