import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class TotalRepeats {

    public static void main(String[] args) throws IOException {
        if (args.length > 0) {
            if (isHelpRequest(args)) { // -help, --help, -h, -?, /?, /h, help, usage, ...
                printHelp();
                return;
            }
            String infile = args[0]; // file path or Folder
            String reffile = "";
            String argsRaw = String.join(" ", args) + " ";   // case-preserving (for paths)
            String s = argsRaw.toLowerCase();                // case-insensitive (for flags/keys)
            int kmer = 19;   // optimal rules: kmer=18-19 seqlen=30...90, kmer=12-18 for short sequences
            int seqlen = 80;
            int gap = kmer + kmer;
            int width = 0;
            int hight = 0;
            int imaged = 10;               // 1...30
            int flanksshow = 0;
            int combine = 0;               // Two or more files analysed simultaneously.
            boolean amask = false;         // masking is performed by Pairwise Sequence Alignment: Repeater2 based 
            boolean maskonly = false;      // The process of masking is performed without the use of clustering or annotation.
            boolean seqshow = false;
            boolean fastclustering = true; // Multithreading clustering 
            boolean readmask = false;       // reading masked FASTA for clustering and annotation
            boolean extupmask = false;      // extraction of the UPPER blocks of the masked chromosome contain unique sequences.
            boolean readgff = false;
            boolean extract = false;        // Split a larger FASTA file into multiple FASTA files, works with single file or folder
            boolean maskfiles = false;      // A comparison analysis of masked files obtained from different software or algorithms.
            boolean ssrdetect = true;
            String outdir = "";             // output folder; default = folder of input file

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            if (s.contains("-out=")) {
                outdir = strArg(argsRaw, s, "-out=");
                File outFolder = new File(outdir);
                if (!outFolder.exists()) {
                    if (outFolder.mkdirs()) {
                        System.out.println("Output folder created: " + outFolder.getAbsolutePath());
                    } else {
                        System.err.println("Cannot create output folder: " + outdir + ". Using default.");
                        outdir = "";
                    }
                } else if (!outFolder.isDirectory()) {
                    System.err.println("Output path is not a folder: " + outdir + ". Using default.");
                    outdir = "";
                } else {
                    System.out.println("Output folder: " + outFolder.getAbsolutePath());
                }
            }

            if (s.contains("-lib=") || s.contains("-ref=")) {
                reffile = s.contains("-lib=")
                        ? strArg(argsRaw, s, "-lib=")
                        : strArg(argsRaw, s, "-ref=");
                File ref = new File(reffile);
                if (!ref.exists() || !ref.isFile()) {
                    System.err.println("Reference file does not exist: " + reffile);
                    reffile = "";
                } else {
                    System.out.println("Reference file found: " + ref.getAbsolutePath());
                }
            }
            if (s.contains("-amask")) { // Pairwise Sequence Alignment: Repeater2 based 
                amask = true;
            }
            if (s.contains("-nossr")) {
                ssrdetect = false;
            }
            if (s.contains("-normal")) {// (Prevent Multithreading clustering) Using fully multithreading significantly accelerates the classification of sequences into individual clusters.
                fastclustering = false;
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
            imaged = clamp(intArg(s, imaged, "imgx="), 5, 30);
            flanksshow = clamp(intArg(s, flanksshow, "flanks=", "flangs="), 0, 1000);
            kmer = intArg(s, kmer, "kmer=");
            seqlen = intArg(s, seqlen, "sln=");
            gap = intArg(s, gap, "gap=");
            if (s.contains("image=")) { // image=40000x5000
                int j = s.indexOf("image=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    String[] d = s.substring(j + 6, x).split("x");
                    if (d.length > 1) {
                        width = StrToInt(d[0]);
                        hight = StrToInt(d[1]);
                        if (hight == 0 || width == 0) {
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
                    if (files == null) {
                        System.err.println("Cannot read directory: " + folder.toString());
                        return;
                    }
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
                        TotalRepeatsCombinedResult(combine, reffile, filelist, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, fastclustering, ssrdetect, outdir);
                        return;
                    }

                    if (extract) {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    ExtractFiles(nfile, outdir);
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
                                    ReadingMaskFile(nfile, reffile, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, fastclustering, ssrdetect, outdir);
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
                                    ReadingUpperMaskFile(nfile, seqlen, gap, outdir);
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
                                    ReadingGffFile(nfile, imaged, width, hight, outdir);
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
                                TotalRepeatsResult(nfile, reffile, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, maskonly, amask, fastclustering, ssrdetect, outdir);
                            } catch (Exception e) {
                                System.err.println("Failed to open file: " + nfile);
                            }
                        }
                    }

                } else {
                    if (maskfiles) {
                        ComparisionMaskFiles(infile, outdir);
                        return;
                    }
                    if (readmask) {
                        ReadingMaskFile(infile, reffile, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, fastclustering, ssrdetect, outdir);
                        return;
                    }
                    if (readgff) {
                        ReadingGffFile(infile, imaged, width, hight, outdir);
                        return;
                    }
                    if (extract) {
                        ExtractFiles(infile, outdir);
                        return;
                    }
                    if (extupmask) {
                        ReadingUpperMaskFile(infile, seqlen, gap, outdir);
                        return;
                    }

                    TotalRepeatsResult(infile, reffile, kmer, seqlen, gap, flanksshow, imaged, seqshow, width, hight, maskonly, amask, fastclustering, ssrdetect, outdir);
                }
            } else {
                // No file or folder to analyse: report the problem and show usage.
                System.err.println();
                System.err.println("ERROR: input file or folder not found: " + infile);
                System.err.println("There is no sequence file to analyse. Provide the input path as the first argument.\n");
                printHelp();
            }
        } else {
            printHelp();
        }
    }

    /**
     * Returns {@code true} if the command line is an explicit request for help
     * or usage information. Recognised case-insensitively, anywhere in the
     * argument list: {@code -help}, {@code --help}, {@code -h}, {@code --h},
     * {@code -?}, {@code --?}, {@code /?}, {@code /h}, {@code /help},
     * {@code help}, {@code ?}, {@code -usage}, {@code --usage}, {@code /usage}.
     */
    private static boolean isHelpRequest(String[] args) {
        if (args == null) {
            return false;
        }
        for (String arg : args) {
            if (arg == null) {
                continue;
            }
            switch (arg.trim().toLowerCase()) {
                case "-help", "--help", "-h", "--h", "-?", "--?", "/?", "/h", "/help", "help", "?", "-usage", "--usage", "/usage" -> {
                    return true;
                }
            }
        }
        return false;
    }

    private static void printHelp() {
        String jar = "java -jar /path/to/TotalRepeats.jar";
        String jarMem = "java -jar -Xms16g -Xmx32g /path/to/TotalRepeats.jar";

        String[] lines = {
            "TotalRepeats (2024-2026) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)",
            "https://github.com/rkalendar/TotalRepeats",
            "",
            "BASIC USAGE:",
            "  " + jar + " <input_file_or_folder> [options]",
            "",
            "HELP:",
            "  -help, --help, -h, -?, /?, /h    Show this help message and exit",
            "                                  (also shown when no input file is given)",
            "",
            "CORE OPTIONS:",
            "  -kmer=<9-21>         K-mer size for repeat detection (default: 19)",
            "  -sln=<int>           Minimum repeat block length in bp (default: 80; can equal kmer)",
            "  -flangs=<int>        Extend repeat flanks by N nucleotides (default: 0)",
            "",
            "IMAGE OPTIONS:",
            "  -image=<WxH>         Output image dimensions, e.g. -image=10000x300",
            "                       (default: auto-determined)",
            "  -imgx=<1-30>         Figure width compression; 1 = max compression, 20 = longest figure",
            "",
            "ANALYSIS MODES:",
            "  -maskonly            Generate repeat mask file only; skip classification, annotation,",
            "                       and visualisation",
            "  -amask               Mask using pairwise sequence alignment (Repeater2-based)",
            "  -seqshow             Extract and output repeat sequences (default: off)",
            "  -nossr               Disable SSR (Simple Sequence Repeat) detection (default: on)",
            "  -quick               Accelerate repeat classification using multithreading (default: off)",
            "",
            "COMPARATIVE / GENOME-WIDE OPTIONS:",
            "  -combine             Genome-wide analysis: each sequence analyzed individually",
            "  -combine2            Genome-wide analysis: all sequences analyzed together",
            "  -combinemask         Genome-wide comparative analysis using masking files as input",
            "  -homology            Comparative Homology Masking: cross-file analysis of homologous",
            "                       regions and repeats between target sequences (e.g. chromosomes)",            
            "",
            "INPUT / OUTPUT OPTIONS:",
            "  -readmask            Load a masking file (single FASTA entry) for clustering and",
            "                       visualisation",
            "  -readgff             Load a GFF file for visualisation",
            "  -extract             Split a single multi-entry FASTA file into individual FASTA files",
            "  -maskscomp           Compare masking files produced by different software or algorithms",
            "  -lib=<path>          Annotate repeats using a database of known repeats or genes",
            "  -out=<path>          Path to output folder (default: current folder)",
            "",
            "EXAMPLES:",
            "  # Standard run with SSR and flanking sequences:",
            "  " + jarMem + " <input> -flangs=100 -seqshow",
            "",
            "  # Custom k-mer, block length, no masking, with flanks:",
            "  " + jarMem + " <input> -kmer=18 -sln=100 -seqshow -flangs=100",
            "",
            "  # Directory input with a reference file:",
            "  " + jarMem + " /data/genomes/T2T-CHM13v2.0/ -ref=/data/test/humsub.ref -out=/data/report/",
            "",
            "  # Large chromosome (>2 GB) — allocate up to 64 GB RAM:",
            "  java -jar -Xms16g -Xmx64g /path/to/TotalRepeats.jar /data/genomes/Viscum_album/ -out=/data/report/",
            "",
            "  # Analyse all FASTA files in a directory:",
            "  java -jar -Xms16g -Xmx32g /path/to/TotalRepeats.jar /data/genomes/Aegilops_tauschii/",
            "",
            "NOTE: For sequences larger than 2 GB, increase heap memory with -Xmx (e.g. -Xmx64g).",};

        for (String line : lines) {
            System.out.println(line);
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

    /**
     * Returns the integer value following the first matching {@code key} (e.g.
     * "kmer=") in the lowercased argument string {@code s}, or {@code def} if
     * none of the keys is present.
     */
    private static int intArg(String s, int def, String... keys) {
        for (String key : keys) {
            int j = s.indexOf(key);
            if (j >= 0) {
                int x = s.indexOf(' ', j);
                if (x > j) {
                    return StrToInt(s.substring(j + key.length(), x));
                }
            }
        }
        return def;
    }

    /**
     * Extracts the case-preserving value following {@code key} from the raw
     * argument string, or {@code null} if the key is absent. The key is matched
     * case-insensitively via {@code rawLower}; the value is taken from
     * {@code raw}.
     */
    private static String strArg(String raw, String rawLower, String key) {
        int j = rawLower.indexOf(key);
        if (j < 0) {
            return null;
        }
        int x = raw.indexOf(' ', j);
        if (x < 0) {
            x = raw.length();
        }
        return raw.substring(j + key.length(), x).trim();
    }

    private static int clamp(int v, int lo, int hi) {
        return (v < lo) ? lo : (v > hi) ? hi : v;
    }

    private static void printFastaFormatHelp() {
        System.out.println("There is no sequence(s).");
        System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
        System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
    }

    private static void TotalRepeatsCombinedResult(int combine, String reffile, String[] filelist,
            int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean seqshow,
            int width, int hight, boolean fst, boolean ssr, String outdir) throws IOException {

        long startTime = System.nanoTime();

        List<String> seqs = new ArrayList<>();
        List<String> names = new ArrayList<>();
        List<String> fnames = new ArrayList<>();

        String combinedFile = getReportFile(filelist, outdir);
        if (combinedFile != null) {
            System.out.println("Combined report path: " + combinedFile);
        } else {
            System.out.println("File list is empty or invalid.");
        }

        for (String nfile : filelist) {
            if (nfile == null) {
                continue;
            }

            try {
                FastaReader rf = (combine == 3)
                        ? FastaReader.fromPathRaw(Paths.get(nfile)) // without normalisation
                        : FastaReader.fromPath(Paths.get(nfile));   // with normalisation

                if (!rf.isEmpty()) {
                    System.out.println("The target file contains " + rf.getSequenceCount() + " sequence(s)");

                    names.addAll(rf.getNames());
                    seqs.addAll(rf.getSequences());

                    if (rf.getSequenceCount() > 1) {
                        for (int i = 0; i < rf.getSequenceCount(); i++) {
                            fnames.add(nfile + "_" + (i + 1) + "_");
                        }
                    } else {
                        fnames.add(nfile);
                    }

                    System.out.println("Target file: " + nfile);
                } else {
                    fnames.add(nfile);
                    printFastaFormatHelp();
                }
            } catch (IOException e) {
                System.err.println("Failed to open file: " + nfile);
            }
        }

        System.out.println("\nRunning...");
        System.out.println("kmer=" + kmer);
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
        s2.SetReportFile(combinedFile);
        s2.SetSSRdetection(ssr);

        Path path = Paths.get(fnms[0]);
        Path parentDir = path.getParent();
        Path effectiveDir = (outdir != null && !outdir.isEmpty()) ? Paths.get(outdir) : parentDir;
        if (effectiveDir == null) {
            System.err.println("Cannot determine output directory for combined result.");
            return;
        }

        if (combine == 3) {
            String fileName = effectiveDir.resolve("combined").toString();
            System.out.println("Combined path: " + fileName);
            s2.SetFileName(fileName);
        }

        if (width > 0 && hight > 0) {
            s2.SetImage(width, hight);
        }

        if (!reffile.isEmpty()) {
            FastaReader ref = FastaReader.fromPath(Paths.get(reffile));
            s2.SetRefSequences(
                    ref.getSequences().toArray(String[]::new),
                    ref.getNames().toArray(String[]::new)
            );
            System.out.println("Reference file=" + reffile);
        }

        switch (combine) {
            case 1 ->
                s2.RunCombine(imgx, fst);
            case 2 ->
                s2.RunCombine2(imgx, fst);
            case 3 ->
                s2.RunCombineMask(imgx, fst);
            case 4 ->
                s2.RunHomologyMasking(imgx);
        }

        long duration = (System.nanoTime() - startTime) / 1_000_000_000;
        System.out.println("Total duration: " + duration + " seconds\n");
    }

    private static void TotalRepeatsResult(String infile, String reffile, int kmer, int seqlen,
            int gap, int flanksshow, int imgx, boolean seqshow, int width, int hight,
            boolean maskonly, boolean amask, boolean fst, boolean ssr, String outdir) {
        try {
            long startTime = System.nanoTime();

            FastaReader rf = FastaReader.fromPath(Paths.get(infile));

            if (rf.isEmpty()) {
                printFastaFormatHelp();
                return;
            }

            System.out.println("\nRunning...");
            System.out.println("kmer=" + kmer);
            System.out.println("Repeat block length =" + seqlen);
            System.out.println("Target file: " + infile);
            if (rf.getSequenceCount() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getSequenceCount());
            }

            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(
                    rf.getSequences().toArray(String[]::new),
                    rf.getNames().toArray(String[]::new)
            );
            s2.SetRepeatLen(kmer, seqlen, gap);
            s2.SetMaskGenerate(maskonly);
            s2.SetShowSeq(seqshow);
            s2.SetSSRdetection(ssr);
            s2.SetFlanks(flanksshow);
            s2.SetFileName(resolveOutputDir(infile, outdir) + File.separator + new File(infile).getName());

            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }

            if (!reffile.isEmpty()) {
                FastaReader ref = FastaReader.fromPath(Paths.get(reffile));
                s2.SetRefSequences(
                        ref.getSequences().toArray(String[]::new),
                        ref.getNames().toArray(String[]::new)
                );
                System.out.println("Reference file=" + reffile);
            }

            if (amask) {
                s2.RunAlignmentMask(imgx, fst);
            } else {
                s2.Run(imgx, fst);
            }

            long duration = (System.nanoTime() - startTime) / 1_000_000_000;
            System.out.println("Total duration: " + duration + " seconds\n");

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    public static long countLower(String text) {
        long lowercase = 0;
        for (int i = 0, n = text.length(); i < n; i++) {
            if (Character.isLowerCase(text.charAt(i))) {
                lowercase++;
            }
        }
        return lowercase;
    }

    private static void ReadingUpperMaskFile(String inputFile, int seqlen, int gap, String outdir) throws IOException {
        try {
            FastaReader rf = FastaReader.fromPathRaw(Paths.get(inputFile));

            if (rf.isEmpty()) {
                printFastaFormatHelp();
                return;
            }

            System.out.println("\nRunning...");
            System.out.println("Target file: " + inputFile);
            if (rf.getSequenceCount() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getSequenceCount());
            }

            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetFileName(resolveOutputDir(inputFile, outdir) + File.separator + new File(inputFile).getName());
            s2.SetSequences(
                    rf.getSequences().toArray(String[]::new),
                    rf.getNames().toArray(String[]::new)
            );
            s2.RunUniquesMaskSaving(gap, seqlen);

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    private static void ExtractFiles(String inputFile, String outdir) {
        String parentDir = resolveOutputDir(inputFile, outdir);

        try (BufferedReader reader = new BufferedReader(new FileReader(inputFile)); BufferedWriter report = new BufferedWriter(new FileWriter(parentDir + File.separator + "report.txt"))) {
            System.out.println("\nRunning...");
            long l = 0;
            long r = 0;
            String line;
            BufferedWriter writer = null;

            while ((line = reader.readLine()) != null) {
                if (line.contains(">")) {
                    if (writer != null) {
                        report.write("Total length (nt)= " + l);
                        report.newLine();
                        report.write("Total masked length (nt)= " + r);
                        report.newLine();
                        report.write("Masked= " + String.format("%.2f", l > 0 ? (100.0 * r) / l : 0.0) + "%\n");
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
                report.write("Masked= " + String.format("%.2f", l > 0 ? (100.0 * r) / l : 0.0) + "%\n");
            }
            System.out.println("File(s) processed successfully.");
        } catch (IOException e) {
            System.err.println("Error processing file: " + inputFile + " — " + e.getMessage());
        }
    }

    private static void ComparisionMaskFiles(String inputFile, String outdir) {
        String parentDir = resolveOutputDir(inputFile, outdir);
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
            System.err.println("Error processing mask comparison: " + e.getMessage());
        }
    }

    private static void ReadingMaskFile(String inputFile, String reffile, int kmer, int seqlen,
            int gap, int flanksshow, int imgx, boolean seqshow, int width, int hight,
            boolean fst, boolean ssr, String outdir) {
        try {
            FastaReader rf = FastaReader.fromPathRaw(Paths.get(inputFile));

            if (rf.isEmpty()) {
                printFastaFormatHelp();
                return;
            }

            System.out.println("\nRunning...");
            System.out.println("Target file: " + inputFile);
            if (rf.getSequenceCount() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getSequenceCount());
            }

            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(
                    rf.getSequences().toArray(String[]::new),
                    rf.getNames().toArray(String[]::new)
            );
            s2.SetRepeatLen(kmer, seqlen, gap);
            s2.SetShowSeq(seqshow);
            s2.SetSSRdetection(ssr);
            s2.SetFlanks(flanksshow);
            s2.SetFileName(resolveOutputDir(inputFile, outdir) + File.separator + new File(inputFile).getName());

            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }

            if (!reffile.isEmpty()) {
                FastaReader ref = FastaReader.fromPath(Paths.get(reffile));
                s2.SetRefSequences(
                        ref.getSequences().toArray(String[]::new),
                        ref.getNames().toArray(String[]::new)
                );
                System.out.println("Reference file=" + reffile);
            }

            s2.RunThroughMask(imgx, fst);

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    private static void ReadingGffFile(String inputGFFfile, int imgx, int width, int hight, String outdir) {
        try {
            System.out.println("\nRunning...");
            System.out.println("GFF file: " + inputGFFfile);
            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetFileName(resolveOutputDir(inputGFFfile, outdir) + File.separator + new File(inputGFFfile).getName());
            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            s2.RunThroughGFF(inputGFFfile, imgx);
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }

    /**
     * Returns the effective output directory for a given input file. Uses
     * outdir if specified; otherwise falls back to the input file's parent
     * folder.
     */
    private static String resolveOutputDir(String inputFile, String outdir) {
        if (outdir != null && !outdir.isEmpty()) {
            return outdir;
        }
        String parent = new File(inputFile).getParent();
        return (parent != null) ? parent : ".";
    }

    private static String getReportFile(String[] filelist, String outdir) {
        if (filelist == null || filelist.length == 0) {
            return "";
        }
        String dir = (outdir != null && !outdir.isEmpty()) ? outdir : resolveOutputDir(filelist[0], "");
        return dir + File.separator + "report";
    }

}
