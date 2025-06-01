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
            int kmer = 19;   //quick search: optimal rules: kmer=19-21 seqlen=30...100, gap=kmer kmer=12-18 for short seqyences
            int seqlen = 90;
            int gap = kmer;
            int width = 0;
            int hight = 0;
            int imaged = 20;    //1...40
            int flanksshow = 0;
            int nkmer = 20;     //0..2-230  nsize=1 - very fast clustering without chain direction detection; nsize=0 - used when ignoring clustering; nsize >2 - complete clustering        
            int combine = 0;
            boolean maskshow = true;
            boolean seqshow = false;
            boolean gffshow = true;
            boolean maskpicture = false;
            boolean sensitivity = false;

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            if (s.contains("ref=")) {
                int j = s.indexOf("ref=");
                int x = s.indexOf(" ", j);
                reffile = s.substring(j + 4, x);
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
                    if (imaged > 40) {
                        imaged = 40;
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
                        SaveResult2(combine, reffile, filelist, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, maskpicture, sensitivity);
                    } else {
                        for (String nfile : filelist) {
                            if (nfile != null) {
                                try {
                                    SaveResult(nfile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, maskpicture, sensitivity);
                                } catch (Exception e) {
                                    System.err.println("Failed to open file: " + nfile);
                                }
                            }
                        }
                    }

                } else {
                    SaveResult(infile, reffile, nkmer, kmer, seqlen, gap, flanksshow, imaged, gffshow, maskshow, seqshow, width, hight, maskpicture, sensitivity);
                }
            }
        } else {
            System.out.println("TotalRepeats (2024-2025) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/TotalRepeats\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("kmer=19\tminimal kmer=12 (default kmer=19)");
            System.out.println("sln=60\trepeat block length (default sln=60), it can be equal to 'kmer'");
            System.out.println("nsize=12\tspeed and sensitivity of sequence clustering: nsize=0 - ignoring clustering;  nsize=1 - very fast clustering without chain direction detection; nsize=12 - default for complete clustering.");
            System.out.println("flangs=100\textend the flanks of the repeat with an appropriate length (100 nt) (default flangs=0)");
            System.out.println("image=10000x300\t (by default, the dimensionality of the image is automatically determined)");
            System.out.println("imgx=3\t (figure width compression, minimum value of imgx=1 (maximum compression), and a value of imgx=10 for the longest figure length)");
            System.out.println("-nomask\tquick generation a new file with masking repeats (default performed)");
            System.out.println("-nogff\tgenerate a GFF file (default performed)");
            System.out.println("-maskpic\tgenerate a image file with masking repeats (default not performed)");
            System.out.println("-seqshow\textract repeat sequences (default not performed)");
            System.out.println("-combine\tthis option is employed in genome-wide comparative analyses (each sequence is analyzed for repeats individually) (default not performed)");
            System.out.println("-combine2\tthis option is employed in genome-wide comparative analyses (all sequences are analyzed together) (default not performed)");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> kmer=18 sln=100 mask=false seqshow=true flanks=100\n");
            System.out.println("Large genome settings:");
            System.out.println("java -jar -Xms16g -Xmx64g \\TotalRepeats\\dist\\TotalRepeats.jar <inputfile> kmer=21 sln=90 \n");
            System.out.println("Analysing all files in the directory:");
            System.out.println("java -jar \\TotalRepeats\\dist\\TotalRepeats.jar \\TotalRepeats\\test\\\n");
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

    private static void SaveResult2(int combine, String reffile, String[] filelist, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean maskshow, boolean seqshow, int width, int hight, boolean maskpic, boolean sensitivity) throws IOException {
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

        System.out.println("Running...");
        System.out.println("kmer=" + kmer);
        System.out.println("Classification index (0-230)=" + nkmer);
        System.out.println("Repeat block length =" + seqlen);
        System.out.println("Shown repeated sequence is " + seqshow);
        if (flanksshow > 0) {
            System.out.println("Flanks around sequence is " + flanksshow);
        }

        String[] nms = names.toArray(String[]::new);
        String[] sqs = seqs.toArray(String[]::new);
        String[] fnms = fnames.toArray(String[]::new);

        TotalRepeatsSearching s2 = new TotalRepeatsSearching();
        s2.SetSequences(sqs, nms);
        s2.SetRepeatLen(kmer, seqlen, gap);
        s2.SetShowSeq(seqshow);
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

    private static void SaveResult(String infile, String reffile, int nkmer, int kmer, int seqlen, int gap, int flanksshow, int imgx, boolean gffshow, boolean maskshow, boolean seqshow, int width, int hight, boolean maskpic, boolean sensitivity) {
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

            System.out.println("Running...");
            System.out.println("kmer=" + kmer);
            System.out.println("Classification index (0-230)=" + nkmer);
            System.out.println("Repeat block length =" + seqlen);
            System.out.println("Target file: " + infile);
            if (rf.getNseq() > 1) {
                System.out.println("Target FASTA sequences = " + rf.getNseq());
            }
            System.out.println("Shown repeated sequence is " + seqshow);
            if (flanksshow > 0) {
                System.out.println("Flanks around sequence is " + flanksshow);
            }

            TotalRepeatsSearching s2 = new TotalRepeatsSearching();
            s2.SetSequences(rf.getSequences(), rf.getNames());

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
                OpenSeqFiles ffiles = new OpenSeqFiles(reffile);
                System.out.println("Reference file=" + reffile);
                s2.SetRefSequences(ffiles.getSeqs(), ffiles.getNames());
            }

            s2.Run(imgx, nkmer, sensitivity);

            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds\n");
        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
    }
}
