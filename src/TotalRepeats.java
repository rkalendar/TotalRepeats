import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class TotalRepeats {

    public static void main(String[] args) throws IOException {
        if (args.length > 0) {
            String infile = args[0]; // file path or Folder
            String reffile = "";
            String s = String.join(" ", args).toLowerCase() + " ";
            String c = String.join(" ", args) + " ";
            int kmer = 19;   //optimal rules: kmer=19-21 seqlen=30...100, gap=kmer kmer=12-18 for short seqyences
            int seqlen = 90;
            int gap = kmer + kmer;
            int width = 0;
            int hight = 0;
            int imaged = 5;    //1...20
            int flanksshow = 0;
            int nkmer = 12;     //0,1,2-230: nsize=0: Used when ignoring clustering; Use size=1 for very fast clustering without chain direction detection; nsize >1: Used for clustering.        
            int combine = 0;
            boolean maskshow = true;
            boolean seqshow = false;
            boolean gffshow = true;
            boolean sensitivity = true;
            boolean readmask = false;
            boolean readgff = false;

            boolean extract = false;   // Split a single FASTA file into multiple FASTA files.
            boolean maskfiles = false; // A comparison analysis of masked files obtained from different software or algorithms.

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            if (s.contains("ref=")) {
                int j = s.toLowerCase().indexOf("ref=");
                int x = s.indexOf(" ", j);
                reffile = c.substring(j + 4, x);
            }
            if (s.contains("readgff")) {
                readgff = true;
            }
            if (s.contains("readmask")) {
                readmask = true;
            }
            if (s.contains("smask")) {
                sensitivity = true;
            }
            if (s.contains("maskscomp")) {
                maskfiles = true;
            }
            if (s.contains("extract")) {
                extract = true;
            }
            if (s.contains("combine")) {
                combine = 1;
            }
            if (s.contains("combine2")) {
                combine = 2;
            }
            if (s.contains("nomask")) {
                maskshow = false;
            }
            if (s.contains("nogff")) {
                gffshow = false;
            }
            if (s.contains("seqshow")) {
                seqshow = true;
            }
            if (s.contains("imgx=")) {
                int j = s.indexOf("imgx=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    imaged = StrToInt(s.substring(j + 5, x));
                    if (imaged < 1) {
                        imaged = 1;
                    }
                    if (imaged > 20) {
                        imaged = 20;
                    }
                }
            }
            if (s.contains("nsize=")) {
                int j = s.indexOf("nsize=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    nkmer = StrToInt(s.substring(j + 6, x));
                }
            }
            if (s.contains("flanks=")) {
                flanksshow = 0;
                int j = s.indexOf("flanks=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    flanksshow = StrToInt(s.substring(j + 7, x));
                }
                if (flanksshow < 0) {
                    flanksshow = 0;
                }
                if (flanksshow > 1000) {
                    flanksshow = 1000;
                }
            }
            if (s.contains("kmer=")) {
                int j = s.indexOf("kmer=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    kmer = StrToInt(s.substring(j + 5, x));
                }
            }
            if (s.contains("sln=")) {
                int j = s.indexOf("sln=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    seqlen = StrToInt(s.substring(j + 4, x));
                }
            }
            if (s.contains("gap=")) {
                int j = s.indexOf("gap=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    gap = StrToInt(s.substring(j + 4, x));
                }
            }
            if (s.contains("image=")) { // image=40000x5000
                int j = s.indexOf("image=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    String[] d = s.substring(j + 6, x).split("x");
                    if (d.length > 0) {
                        width = StrToInt(d[0]);
                        hight = StrToInt(d[1]);
                        if (hight == 0 | width == 0) {
                            width = 0;
                            hight = 0;
                        }
                    }
                }
            }

            File folder = new File(infile);
            if (folder.exists() && (folder.isDirectory() || folder.isFile())) {
                if (folder.isDirectory()) {
                    File[] files = folder.listFiles();
                    if (files.length == 0) {
                        System.err.println("No files: " + folder.toString());
                        return;
                    }
                    int k = -1;
                    String[] filelist = new String[files.length];
                    for (File file : files) {
                        if (file.isFile()) {
                            filelist[++k] = file.getAbsolutePath();
                        }
                    }
                    if (combine > 0) {
                        SaveResult2(combine, reffile, filelist, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, sensitivity);
                    } else {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    SaveResult(nfile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, sensitivity);
                                } catch (Exception e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                    }

                } else {
                    if (maskfiles) {
                        ComparisionMaskFiles(infile);
                        return;
                    }
                    if (readmask) {
                        ReadingMaskFile(infile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, seqshow, width, hight);
                        return;
                    }
                    if (extract) {
                        ExtractFiles(infile);
                        return;
                    }
                    SaveResult(infile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, sensitivity);
                }
            }
        } else {
            System.out.println("TotalRepeats (2024-2025) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/TotalRepeats\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("kmer=19\t kmer=9-21 (default kmer=19)");
            System.out.println("sln=90\trepeat block length (default sln=90), it can be equal to 'kmer'");
            System.out.println("nsize=12\tspeed and sensitivity of sequence clustering: nsize=0 - ignoring clustering;  nsize=1 - very fast clustering without chain direction detection; nsize=12 - default for complete clustering.");
            System.out.println("-smask\ta more sensitive method for identifying repetitive sequences (default not performed)");
            System.out.println("flangs=100\textend the flanks of the repeat with an appropriate length (100 nt) (default flangs=0)");
            System.out.println("image=10000x300\t (by default, the dimensionality of the image is automatically determined)");
            System.out.println("imgx=5\t (figure width compression, minimum value of imgx=1 (maximum compression), and a value of imgx=20 for the longest figure length)");
            System.out.println("-nomask\tquick generation a new file with masking repeats (default performed)");
            System.out.println("-nogff\tgenerate a GFF file (default performed)");
            System.out.println("-maskpic\tgenerate a image file with masking repeats (default not performed)");
            System.out.println("-seqshow\textract repeat sequences (default not performed)");
            System.out.println("-combine\tthis option is employed in genome-wide comparative analyses (each sequence is analyzed for repeats individually) (default not performed)");
            System.out.println("-combine2\tthis option is employed in genome-wide comparative analyses (all sequences are analyzed together) (default not performed)");
            System.out.println("-ref=target_file_path\tthe application enables annotation of repeats using a database of known repeats/genes (default not performed)");
            System.out.println("-readmask\ttransfer the masking file to the software, which will then be used for clustering repeats and visualisation. The file contains only one FASTA entry.");
            System.out.println("-readgff\ttransfer the GFF file to the software, which will then be used for visualisation.");
            System.out.println("-extract\tSplit a single FASTA file into multiple FASTA files.");
            System.out.println("-maskscomp\tA comparison analysis of masked files obtained from different software or algorithms.");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> kmer=18 sln=100 mask=false seqshow=true flanks=100\n");
            System.out.println("java -jar -Xms16g -Xmx64g \\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\T2T-CHM13v2.0\\ -ref=C:\\TotalRepeats\\test\\humsub.ref\n");
            System.out.println("Large chromosome usage (1 GB): you will need to show the program to use more RAM, up to 128-256 GB of memory:\n");
            System.out.println("java -jar -Xms32g -Xmx128g C:\\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\Pleurodeles_waltl\\ \n");
            System.out.println("Analysing all files in the directory:");
            System.out.println("java -jar -Xms16g -Xmx64g C:\\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\T2T-CHM13v2.0\\\n");
            System.out.println("java -jar C:\\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\Shigella\\ -combine -nomask -nogff \n");
        }
    }

    public static int StrToInt(String str) {
        StringBuilder r = new StringBuilder();
        int z = 0;
        r.append(0);
        for (int i = 0; i < str.length(); i++) {
            char chr = str.charAt(i);
            if (chr > 47 && chr < 58) {
                r.append(chr);
                z++;
                if (z > 10) {
                    break;
                }
            }
            if (chr == '.' || chr == ',') {
                break;
            }
        }
        return (Integer.parseInt(r.toString()));
    }

    private static void SaveResult2(int combine, String reffile, String[] filelist, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean maskshow, boolean seqshow, int width, int hight, boolean sensitivity) throws IOException {
        long startTime = System.nanoTime();
        List<String> seqs = new ArrayList<>();
        List<String> names = new ArrayList<>();
        List<String> fnames = new ArrayList<>();

        for (String nfile : filelist) {
            if (nfile != null) {
                try {
                    byte[] binaryArray = Files.readAllBytes(Paths.get(nfile));
                    ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);

                    if (rf.getNseq() > 0) {
                        System.out.println("Target FASTA sequences = " + rf.getNseq());
                        names.addAll(Arrays.asList(rf.getNames()));
                        seqs.addAll(Arrays.asList(rf.getSequences()));
                        fnames.add(nfile);
                    }
                    if (rf.getNseq() == 0) {
                        System.out.println("There is no sequence(s).");
                        System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                        System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                    } else {
                        System.out.println("Target file: " + nfile);
                    }
                } catch (IOException e) {
                    System.err.println("Failed to open file: " + nfile);
                }
            }
        }

        System.out.println("\nRunning...");
        System.out.println("kmer=" + kmer);
        System.out.println("Classification index (0-230)=" + nkmer);
        System.out.println("Repeat block length =" + seqlen);

        String[] nms = names.toArray(String[]::new);
        String[] sqs = seqs.toArray(String[]::new);
        String[] fnms = fnames.toArray(String[]::new);

        TotalRepeatsSearching s2 = new TotalRepeatsSearching();
        s2.SetSequences(sqs, nms);
        s2.SetRepeatLen(kmer, seqlen, gap);
        s2.SetShowSeq(seqshow);
        s2.SetFlanks(flanksshow);
        s2.SetMasked(maskshow);
        s2.SetGFF(gffshow);
        s2.SetFileNames(fnms);

        Path path = Paths.get(fnms[0]);
        Path parentDir = path.getParent();
        if (parentDir != null) {
            String fileName = "combined";
            Path combinedPath = parentDir.resolve(fileName);
            fileName = combinedPath.toString();
            System.out.println("Combined path: " + fileName);
            s2.SetFileName(fileName);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (reffile.length() > 0) {
                OpenSeqFiles fastafile = new OpenSeqFiles(reffile);
                System.out.println("Reference file=" + reffile);
                s2.SetRefSequences(fastafile.getSeqs(), fastafile.getNames());
            }
            if (combine == 1) {
                s2.RunCombine(imgx, nkmer, sensitivity);
            }
            if (combine == 2) {
                s2.RunCombine2(imgx, nkmer, sensitivity);
            }

        }
        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Total duration: " + duration + " seconds\n");
    }

    private static void SaveResult(String infile, String reffile, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean maskshow, boolean seqshow, int width, int hight, boolean sensitivity) {
        try {
            long startTime = System.nanoTime();
            byte[] binaryArray = Files.readAllBytes(Paths.get(infile));
            ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return;
            }
            System.out.println("\nRunning...");
            System.out.println("kmer=" + kmer);
            System.out.println("Classification index (0-230)=" + nkmer);
            System.out.println("Repeat block length =" + seqlen);
            if (sensitivity) {
                System.out.println("Sensitive detection enabled.");
            }

            System.out.println("Target file: " + infile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }
            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetRepeatLen(kmer, seqlen, gap);
            s2.SetShowSeq(seqshow);
            s2.SetFlanks(flanksshow);
            s2.SetMasked(maskshow);
            s2.SetGFF(gffshow);
            s2.SetFileName(infile);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (reffile.length() > 0) {
                OpenSeqFiles fastafile = new OpenSeqFiles(reffile);
                System.out.println("Reference file=" + reffile);
                s2.SetRefSequences(fastafile.getSeqs(), fastafile.getNames());
            }

            s2.Run(imgx, nkmer, sensitivity);
            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Total duration: " + duration + " seconds\n");
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    public static long countLower(String text) {
        long lowercase = 0;
        for (char c : text.toCharArray()) {
            if (Character.isLowerCase(c)) {
                lowercase++;
            }
        }
        return lowercase;
    }

    private static void ExtractFiles(String inputFile) {
        File input = new File(inputFile);
        String parentDir = input.getParent();
        if (parentDir == null) {
            parentDir = ".";
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
            System.out.println("\nRunning...");
            long l = 0;
            long r = 0;
            String line;
            BufferedWriter writer = null;
            BufferedWriter report = new BufferedWriter(new FileWriter(parentDir + File.separator + "report.txt"));

            while ((line = reader.readLine()) != null) {
                if (line.contains(">")) {
                    if (writer != null) {
                        report.write("Total length (nt)= " + l);
                        report.newLine();
                        report.write("Total masked length (nt)= " + r);
                        report.newLine();
                        report.write("Masked= " + String.format("%.2f", (float) ((r * 100) / l)) + "%\n");
                        writer.close();
                    }
                    String[] s = line.split(">");
                    if (s.length > 1) {
                        String[] s2 = s[1].trim().split(" ");
                        String outputFile = parentDir + File.separator + s2[0] + ".fasta";
                        writer = new BufferedWriter(new FileWriter(outputFile));
                        writer.write(line);
                        writer.newLine();
                        report.write(line);
                        report.newLine();
                        l = 0;
                        r = 0;
                    }
                } else {
                    if (writer != null) {
                        writer.write(line);
                        writer.newLine();
                        l += line.length();
                        r += countLower(line);
                    }
                }
            }
            if (writer != null) {
                writer.close();
                report.write("Total length= " + l);
                report.newLine();
                report.write("Total masked length= " + r);
                report.newLine();
                report.write("Masked= " + String.format("%.2f", (float) ((r * 100) / l)) + "%\n");
                report.close();
            }
            System.out.println("File(s) processed successfully.");
        } catch (IOException e) {

        }
    }

    private static void ComparisionMaskFiles(String inputFile) {
        File input = new File(inputFile);
        String parentDir = input.getParent();
        if (parentDir == null) {
            parentDir = ".";
        }
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFile))) {
            System.out.println("\nRunning...");
            String line;
            StringBuilder sb = new StringBuilder();
            String saveFile = parentDir + File.separator + "result.txt";
            while ((line = reader.readLine()) != null) {
                String[] s = line.trim().split("[\t ]+");
                if (s.length > 1) {
                    String inputFile1 = parentDir + File.separator + s[0].trim();
                    String inputFile2 = parentDir + File.separator + s[1].trim();

                    FileCaseComparator fc = new FileCaseComparator();
                    sb.append(fc.AnalysisFiles(inputFile1, inputFile2)).append("\n\n");
                }
            }
            try (BufferedWriter report = new BufferedWriter(new FileWriter(saveFile))) {
                report.write(sb.toString());
            }
            System.out.println("Files processed successfully.");
        } catch (IOException e) {
        }
    }

    private static void ReadingMaskFile(String inputFile, String reffile, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean seqshow, int width, int hight) {
        try {
            long startTime = System.nanoTime();
            byte[] binaryArray = Files.readAllBytes(Paths.get(inputFile));
            ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return;
            }
            System.out.println("\nRunning...");
            System.out.println("Target file: " + inputFile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }
            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetRepeatLen(kmer, seqlen, gap);
            s2.SetShowSeq(seqshow);
            s2.SetFlanks(flanksshow);
            s2.SetGFF(gffshow);
            s2.SetFileName(inputFile);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (reffile.length() > 0) {
                OpenSeqFiles fastafile = new OpenSeqFiles(reffile);
                System.out.println("Reference file=" + reffile);
                s2.SetRefSequences(fastafile.getSeqs(), fastafile.getNames());
            }

            s2.RunThroughMask(imgx, nkmer, inputFile);
            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Total duration: " + duration + " seconds\n");
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }
}
