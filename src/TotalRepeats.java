import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
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
            int kmer = 19;   // optimal rules: kmer=18-19 seqlen=30...90, kmer=12-18 for short sequences
            int seqlen = 80;
            int gap = kmer + kmer;
            int width = 0;
            int hight = 0;
            int imaged = 10;           // 1...20
            int flanksshow = 0;
            int nkmer = 49;           // 0,1,2-230: nsize=0: Used when ignoring clustering; Use size=1 for very fast clustering without chain direction detection; nsize >2: Used for clustering.        
            int combine = 0;          // Two or more files analysed simultaneously.
            boolean amask = false;    // masking is performed by Pairwise Sequence Alignment: Repeater2 based 
            boolean maskonly = false; // The process of masking is performed without the use of clustering or annotation.
            boolean seqshow = false;
            boolean readmask = false;  // reading masked FASTA for clustering and annotation
            boolean extupmask = false; // extraction of the UPPER blocks of the masked chromosome contain unique sequences.
            boolean readgff = false;
            boolean extract = false;   // Split a larger FASTA file into multiple FASTA files, works with single file or folder
            boolean maskfiles = false; // A comparison analysis of masked files obtained from different software or algorithms.

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            if (s.contains("-lib=")) {
                int j = s.toLowerCase().indexOf("lib=");
                int x = s.indexOf(" ", j);
                if (x == -1) {
                    x = s.length();
                }
                reffile = s.substring(j + 4, x).trim();
                File ref = new File(reffile);
                if (!ref.exists() || !ref.isFile()) {
                    System.err.println("Reference file does not exist: " + reffile);
                    reffile = "";
                } else {
                    System.out.println("Reference file found: " + ref.getAbsolutePath());
                }
            }
            if (s.contains("amask")) { // Pairwise Sequence Alignment: Repeater2 based 
                amask = true;
            }
            if (s.contains("extunique")) { // -seqlen>100
                extupmask = true;
            }
            if (s.contains("maskonly")) { // -nsize=0
                maskonly = true;
            }
            if (s.contains("readgff")) {
                readgff = true;
            }
            if (s.contains("readmask")) {
                readmask = true;
            }
            if (s.contains("maskscomp")) {
                maskfiles = true;
            }
            if (s.contains("extract")) {
                extract = true;
            }
            if (s.contains("combine")) {
                combine = 1;
                if (s.contains("combine2")) {
                    combine = 2;
                }
                if (s.contains("combinemask")) {
                    combine = 3;
                }
            }
            if (s.contains("homology")) {
                combine = 4;
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
                        TotalRepeatsCombinedResult(combine, reffile, filelist, nkmer, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight);
                        return;
                    }

                    if (extract) {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    ExtractFiles(nfile);
                                } catch (Exception e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                        return;
                    }

                    if (readmask) {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    ReadingMaskFile(nfile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight);
                                } catch (Exception e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                        return;
                    }

                    if (extupmask) {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    ReadingUpperMaskFile(nfile, seqlen, gap);
                                } catch (IOException e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                        return;
                    }

                    if (readgff) {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    ReadingGffFile(nfile, imaged, width, hight);
                                } catch (Exception e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                        return;
                    }

                    for (String nfile : filelist) {
                        if (nfile != null) {
                            try {
                                TotalRepeatsResult(nfile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, maskonly, amask);
                            } catch (Exception e) {
                                System.err.println("Failed to open file: " + nfile);
                            }
                        }
                    }

                } else {
                    if (maskfiles) {
                        ComparisionMaskFiles(infile);
                        return;
                    }
                    if (readmask) {
                        ReadingMaskFile(infile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight);
                        return;
                    }
                    if (readgff) {
                        ReadingGffFile(infile, imaged, width, hight);
                        return;
                    }
                    if (extract) {
                        ExtractFiles(infile);
                        return;
                    }
                    if (extupmask) {
                        ReadingUpperMaskFile(infile, seqlen, gap);
                        return;
                    }

                    TotalRepeatsResult(infile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, maskonly, amask);
                }
            }
        } else {
            System.out.println("TotalRepeats (2024-2025) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/TotalRepeats\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar /data/user/dist/TotalRepeats.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("kmer=19\t kmer=9-21 (default kmer=19)");
            System.out.println("sln=90\trepeat block length (default sln=90), it can be equal to 'kmer'");
            System.out.println("nsize=12\tspeed and sensitivity of sequence clustering: nsize=0 - ignoring clustering;  nsize=1 - very fast clustering without chain direction detection; nsize=12 - default for complete clustering.");
            System.out.println("-flangs=100\textend the flanks of the repeat with an appropriate length (100 nt) (default flangs=0)");
            System.out.println("-image=10000x300\t (by default, the dimensionality of the image is automatically determined)");
            System.out.println("-imgx=5\t (figure width compression, minimum value of imgx=1 (maximum compression), and a value of imgx=20 for the longest figure length)");
            System.out.println("-maskonly\\tit only generates the mask file; classification, annotation and visualisation are not performed");
            System.out.println("-amask\\tmasking is performed by pairwise sequence alignment (Repeater2 based)");
            System.out.println("-seqshow\textract repeat sequences (not default)");
            System.out.println("-combine\tthis option is employed in genome-wide comparative analyses (each sequence is analyzed for repeats individually)");
            System.out.println("-combine2\tthis option is employed in genome-wide comparative analyses (all sequences are analyzed together)");
            System.out.println("-homology\tComparative Homology Masking. It performs comparative analysis of individual sequences (chromosomes) using multiple files to analyse homologous regions (and repeats) between target sequences.");
            System.out.println("-combinemask\tthis option is used for genome-wide comparative analyses, for which masking files serve as the input data");
            System.out.println("-readmask\ttransfer the masking file to the software, which will then be used for clustering repeats and visualisation. The file contains only one FASTA entry.");
            System.out.println("-readgff\ttransfer the GFF file to the software, which will then be used for visualisation.");
            System.out.println("-extract\tSplit a single FASTA file into multiple FASTA files.");
            System.out.println("-maskscomp\tA comparison analysis of masked files obtained from different software or algorithms.");
            System.out.println("-lib=target_file_path\tthe application enables annotation of repeats using a database of known repeats/genes");
            System.out.println("java -jar -Xms16g -Xmx64g /data/user/dist/TotalRepeats.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("java -jar -Xms16g -Xmx64g /data/user/dist/TotalRepeats.jar <inputfile> kmer=18 sln=100 mask=false seqshow=true flanks=100\n");
            System.out.println("java -jar -Xms16g -Xmx64g /data/user/dist/TotalRepeats.jar /data/user/genomes/T2T-CHM13v2.0/ -ref=/data/user/test/humsub.ref\n");
            System.out.println("Large chromosome usage (>2 GB): you will need to show the program to use more RAM, up to 256 GB of memory:\n");
            System.out.println("java -jar -Xms64g -Xmx256g /data/user/dist/TotalRepeats.jar /data/user/genomes/Viscum_album/ \n");
            System.out.println("Analysing all files in the directory:");
            System.out.println("java -jar -Xms64g -Xmx128g /data/user/dist/TotalRepeats.jar /data/user/genomes/Aegilops_tauschii/ \n");
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

    private static void TotalRepeatsCombinedResult(int combine, String reffile, String[] filelist, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean seqshow, int width, int hight) throws IOException {
        long startTime = System.nanoTime();
        List<String> seqs = new ArrayList<>();
        List<String> names = new ArrayList<>();
        List<String> fnames = new ArrayList<>();

        for (String nfile : filelist) {
            if (nfile != null) {
                ReadingSequencesFiles rf;
                try {
                    if (combine == 3) {
                        rf = ReadingSequencesFiles.readMasking(Paths.get(nfile));
                    } else {
                        rf = new ReadingSequencesFiles(Paths.get(nfile));
                    }

                    if (rf.getNseq() > 0) {
                        System.out.println("The target file contains " + rf.getNseq() + " sequence(s)");
                        names.addAll(Arrays.asList(rf.getNames()));
                        seqs.addAll(Arrays.asList(rf.getSequences()));
                        for (int i = 0; i < rf.getNseq(); i++) {
                            fnames.add(nfile + "_1");
                        }

                    }
                    if (rf.getNseq() == 0) {
                        fnames.add(nfile);
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
        s2.SetFileNames(fnms);

        Path path = Paths.get(fnms[0]);
        Path parentDir = path.getParent();
        if (parentDir != null) {

            if (combine == 3) {
                String fileName = "combined";
                Path combinedPath = parentDir.resolve(fileName);
                fileName = combinedPath.toString();
                System.out.println("Combined path: " + fileName);
                s2.SetFileName(fileName);
            }

            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (reffile.length() > 0) {
                ReadingSequencesFiles fastafile = new ReadingSequencesFiles(Paths.get(reffile));
                s2.SetRefSequences(fastafile.getSequences(), fastafile.getNames());
                System.out.println("Reference file=" + reffile);
            }
            if (combine == 1) {
                s2.RunCombine(imgx, nkmer);
            }
            if (combine == 4) {
                s2.RunHomologyMasking(imgx, nkmer);
            }
            if (combine == 3) {
                s2.RunCombineMask(imgx, nkmer);
            }
            if (combine == 2) {
                s2.RunCombine2(imgx, nkmer);
            }
        }
        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Total duration: " + duration + " seconds\n");
    }

    private static void TotalRepeatsResult(String infile, String reffile, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean seqshow, int width, int hight, boolean maskonly, boolean amask) {
        try {
            long startTime = System.nanoTime();
            ReadingSequencesFiles rf = new ReadingSequencesFiles(Paths.get(infile));
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
            System.out.println("Target file: " + infile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }
            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetRepeatLen(kmer, seqlen, gap);
            s2.SetMaskGenerate(maskonly);
            s2.SetShowSeq(seqshow);
            s2.SetFlanks(flanksshow);
            s2.SetFileName(infile);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (reffile.length() > 0) {
                ReadingSequencesFiles fastafile = new ReadingSequencesFiles(Paths.get(reffile));
                s2.SetRefSequences(fastafile.getSequences(), fastafile.getNames());
                System.out.println("Reference file=" + reffile);
            }
            if (amask) {
                s2.RunAlignmentMask(imgx, nkmer);
            } else {
                s2.Run(imgx, nkmer);
            }
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

    private static void ReadingUpperMaskFile(String inputFile, int seqlen, int gap) throws IOException {
        try {
            ReadingSequencesFiles rf = ReadingSequencesFiles.readMasking(Paths.get(inputFile));
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
            s2.SetFileName(inputFile);
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.RunUniquesMaskSaving(gap, seqlen);

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }

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

                    FileMasksComparator fc = new FileMasksComparator();
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

    private static void ReadingMaskFile(String inputFile, String reffile, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean seqshow, int width, int hight) {
        try {
            ReadingSequencesFiles rf = ReadingSequencesFiles.readMasking(Paths.get(inputFile));
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
            s2.SetFileName(inputFile);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            if (reffile.length() > 0) {
                ReadingSequencesFiles fastafile = new ReadingSequencesFiles(Paths.get(reffile));
                s2.SetRefSequences(fastafile.getSequences(), fastafile.getNames());
                System.out.println("Reference file=" + reffile);
            }
            s2.RunThroughMask(imgx, nkmer);
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    private static void ReadingGffFile(String inputGFFfile, int imgx, int width, int hight) {
        try {
            System.out.println("\nRunning...");
            System.out.println("GFF file: " + inputGFFfile);
            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetFileName(inputGFFfile);
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            s2.RunThroughGFF(inputGFFfile, imgx);
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }
}
