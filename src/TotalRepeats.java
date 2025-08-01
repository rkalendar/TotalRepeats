import java.io.IOException;
import java.io.File;
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
            int kmer = 19;   //optimal rules: kmer=19-21 seqlen=30...100, gap=kmer kmer=12-18 for short seqyences
            int seqlen = 90;
            int gap = kmer;
            int width = 0;
            int hight = 0;
            int imaged = 5;    //1...20
            int flanksshow = 0;
            int nkmer = 12;     //0,1,2-230: nsize=0: Used when ignoring clustering; Use size=1 for very fast clustering without chain direction detection; nsize >1: Used for clustering.        
            int combine = 0;
            boolean maskshow = true;
            boolean seqshow = false;
            boolean gffshow = true;
            boolean maskpicture = false;
            boolean sensitivity = false;
            boolean quickmask = false; // a little faster, less RAM used and a slight deviation

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            if (s.contains("ref=")) {
                int j = s.indexOf("ref=");
                int x = s.indexOf(" ", j);
                reffile = s.substring(j + 4, x);
            }
            if (s.contains("-quick")) {
                quickmask = true;
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
            if (s.contains("maskpic")) {
                maskpicture = true;
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
                        SaveResult2(combine, reffile, filelist, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, maskpicture, sensitivity, quickmask);
                    } else {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    SaveResult(nfile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, maskpicture, sensitivity, quickmask);
                                } catch (Exception e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                    }

                } else {
                    SaveResult(infile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, maskpicture, sensitivity, quickmask);
                }
            }
        } else {
            System.out.println("TotalRepeats (2024-2025) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/TotalRepeats\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("kmer=19\tminimal kmer=12 (default kmer=19)");
            System.out.println("sln=90\trepeat block length (default sln=90), it can be equal to 'kmer'");
            System.out.println("nsize=12\tspeed and sensitivity of sequence clustering: nsize=0 - ignoring clustering;  nsize=1 - very fast clustering without chain direction detection; nsize=12 - default for complete clustering.");
            System.out.println("flangs=100\textend the flanks of the repeat with an appropriate length (100 nt) (default flangs=0)");
            System.out.println("image=10000x300\t (by default, the dimensionality of the image is automatically determined)");
            System.out.println("imgx=5\t (figure width compression, minimum value of imgx=1 (maximum compression), and a value of imgx=20 for the longest figure length)");
            System.out.println("-nomask\tquick generation a new file with masking repeats (default performed)");
            System.out.println("-nogff\tgenerate a GFF file (default performed)");
            System.out.println("-maskpic\tgenerate a image file with masking repeats (default not performed)");
            System.out.println("-seqshow\textract repeat sequences (default not performed)");
            System.out.println("-combine\tthis option is employed in genome-wide comparative analyses (each sequence is analyzed for repeats individually) (default not performed)");
            System.out.println("-combine2\tthis option is employed in genome-wide comparative analyses (all sequences are analyzed together) (default not performed)");
            System.out.println("-quick\tthis flag accelerates masking slightly and uses less RAM memory on a computer. It is recommended for large chromosomes (default not performed)");
            System.out.println("-ref=target_file_path\tthe application enables annotation of repeats using a database of known repeats/genes (default not performed)");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> kmer=18 sln=100 mask=false seqshow=true flanks=100\n");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\T2T-CHM13v2.0\\ -ref=C:\\TotalRepeats\\test\\humsub.ref\n");
            System.out.println("Large chromosome usage (1 GB): you will need to show the program to use more RAM, up to 64–128 GB of memory:\n");
            System.out.println("java -jar -Xms32g C:\\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\T2T-CHM13v2.0\\ \n");
            System.out.println("Analysing all files in the directory:");
            System.out.println("java -jar C:\\TotalRepeats\\dist\\TotalRepeats.jar E:\\Genomes\\T2T-CHM13v2.0\\\n");
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

    private static void SaveResult2(int combine, String reffile, String[] filelist, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean maskshow, boolean seqshow, int width, int hight, boolean maskpic, boolean sensitivity, boolean quickmask) throws IOException {
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
        s2.SetQuickSearch(quickmask);
        s2.SetFlanks(flanksshow);
        s2.SetMasked(maskshow);
        s2.SetMaskedPicture(maskpic);
        s2.SetGFF(gffshow);
        s2.SetFileNames(fnms);

        Path path = Paths.get(fnms[0]);
        Path parentDir = path.getParent();
        if (parentDir != null) {
            String fileName = "result";
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
                s2.Run2(imgx, nkmer, sensitivity);
            }
            if (combine == 2) {
                s2.Run3(imgx, nkmer, sensitivity);
            }

        }
        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Time taken: " + duration + " seconds\n");
    }

    private static void SaveResult(String infile, String reffile, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean maskshow, boolean seqshow, int width, int hight, boolean maskpic, boolean sensitivity, boolean quickmask) {
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
            System.out.println("Target file: " + infile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }

            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetQuickSearch(quickmask);
            s2.SetRepeatLen(kmer, seqlen, gap);
            s2.SetShowSeq(seqshow);
            s2.SetFlanks(flanksshow);
            s2.SetMasked(maskshow);
            s2.SetMaskedPicture(maskpic);
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
            System.out.println("Time taken: " + duration + " seconds\n");
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }
}
