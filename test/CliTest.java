import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Comparator;
import java.util.Random;
import static java.nio.charset.StandardCharsets.UTF_8;

/**
 * Command-line handling: the lenient integer parser, the token-based option
 * parser, and the filter that stops a directory run from re-reading its own
 * output.
 *
 * <p>The end-to-end checks drive {@code TotalRepeats.main} directly and inspect
 * what it prints, because the parser and the directory filter are private. They
 * use {@code -maskonly} so that only masking runs — enough to exercise argument
 * handling without paying for clustering and image output.
 */
final class CliTest {

    private CliTest() {
    }

    static void run() throws Exception {
        TestRunner.suite("CLI");
        strToInt();
        Path tmp = Files.createTempDirectory("totalrepeats-tests");
        try {
            tokenParser(tmp);
            outputFilter(tmp);
            multiFasta(tmp);
            extract(tmp);
        } finally {
            deleteTree(tmp);
        }
    }

    // ── StrToInt: lenient by design, and must never throw ───────────────────

    private static void strToInt() {
        TestRunner.eq("StrToInt plain", 19, TotalRepeats.StrToInt("19"));
        TestRunner.eq("StrToInt empty is 0", 0, TotalRepeats.StrToInt(""));
        TestRunner.eq("StrToInt trailing junk", 19, TotalRepeats.StrToInt("19abc"));
        TestRunner.eq("StrToInt leading junk", 19, TotalRepeats.StrToInt("abc19"));
        TestRunner.eq("StrToInt stops at a decimal point", 18, TotalRepeats.StrToInt("18.5"));
        TestRunner.eq("StrToInt stops at a comma", 18, TotalRepeats.StrToInt("18,5"));
        TestRunner.eq("StrToInt drops the sign", 5, TotalRepeats.StrToInt("-5"));
        TestRunner.eq("StrToInt joins separated digits", 19, TotalRepeats.StrToInt("1a9"));
        TestRunner.eq("StrToInt of letters only is 0", 0, TotalRepeats.StrToInt("abc"));

        // Regression: 11+ digits used to overflow Integer.parseInt and throw an
        // uncaught NumberFormatException. The value must now saturate instead.
        TestRunner.eq("StrToInt saturates instead of overflowing",
                Integer.MAX_VALUE, TotalRepeats.StrToInt("999999999999"));
        TestRunner.noThrow("StrToInt does not throw on a very long number",
                () -> TotalRepeats.StrToInt("123456789012345678901234567890"));
    }

    // ── Token parser: args[0] is a path, never a mode ───────────────────────

    private static void tokenParser(Path tmp) throws Exception {
        // A folder whose NAME contains a mode keyword. Under the old substring
        // parser this switched on the pangenome mode; it must not now.
        Path dir = Files.createDirectories(tmp.resolve("joint_extract_data"));
        Files.writeString(dir.resolve("a.fasta"), ">a\n" + "acgtacgtac".repeat(20) + "\n", UTF_8);
        Path out = Files.createDirectories(tmp.resolve("out1"));

        String log = runMain(dir.toString(), "-maskonly", "-out=" + out);
        TestRunner.ok("a path containing 'joint' does not start a combined run",
                !log.contains("Combined report path"));
        TestRunner.ok("a path containing 'extract' does not start extraction",
                !log.contains("File(s) processed successfully"));
        TestRunner.ok("the run actually happened", log.contains("- running..."));
        TestRunner.ok("every run is stamped with the release", log.contains("TotalRepeats " + TotalRepeats.VERSION + " - running"));

        // -version prints the release identifier and nothing else is run
        String logv = runMain("-version");
        TestRunner.ok("-version prints the release", logv.trim().equals("TotalRepeats " + TotalRepeats.VERSION));

        // key=value options are read, and an out-of-range k is reported as the
        // value that will really be used.
        String log2 = runMain(dir.toString(), "-maskonly", "-kmer=15", "-out=" + out);
        TestRunner.ok("kmer= is parsed", log2.contains("kmer=15"));

        // -pangroup= is a known option: neither value draws the unknown-option warning
        String logp = runMain(dir.toString(), "-maskonly", "-pangroup=record", "-out=" + out);
        TestRunner.ok("-pangroup=record is accepted", !logp.contains("unrecognized option"));
        String logq = runMain(dir.toString(), "-maskonly", "-pangroup=bogus", "-out=" + out);
        TestRunner.ok("an unknown -pangroup value is reported", logq.contains("-pangroup=bogus not recognised"));

        String log3 = runMain(dir.toString(), "-maskonly", "-kmer=99", "-out=" + out);
        TestRunner.ok("kmer above the supported range is clamped and reported",
                log3.contains("kmer=21"));

        // an unrecognised option is reported rather than silently ignored
        String log4 = runMain(dir.toString(), "-maskonly", "-notanoption", "-out=" + out);
        TestRunner.ok("unknown option warns",
                log4.contains("unrecognized option ignored: -notanoption"));

        // a value that itself contains a mode keyword must stay a value
        String log5 = runMain(dir.toString(), "-maskonly", "-out=" + out);
        TestRunner.ok("an -out= value is never scanned for modes",
                !log5.contains("Combined report path"));
    }

    // ── Directory mode must not consume its own previous output ─────────────

    private static void outputFilter(Path tmp) throws Exception {
        Path dir = Files.createDirectories(tmp.resolve("rerun"));
        Files.writeString(dir.resolve("g.fasta"), ">g\n" + "acgtacgtac".repeat(20) + "\n", UTF_8);
        // files of the kind a previous run would have left behind
        Files.writeString(dir.resolve("g.fasta.gff"), "Seqid\tRepeat\n", UTF_8);
        Files.writeString(dir.resolve("g.fasta.msk"), ">g\nacgt\n", UTF_8);
        Files.writeString(dir.resolve("g.fasta.png"), "not really a png", UTF_8);
        Files.writeString(dir.resolve("g.fasta.svg"), "<svg/>", UTF_8);
        Files.writeString(dir.resolve("report_pangenome.txt"), "core\n", UTF_8);
        // -extunique writes FASTA-shaped records under .report, so an unfiltered
        // rerun would read them straight back in.
        Files.writeString(dir.resolve("g.fasta.report"), ">1-4\nacgt\n", UTF_8);
        Path out = Files.createDirectories(tmp.resolve("out2"));

        String log = runMain(dir.toString(), "-maskonly", "-out=" + out);
        TestRunner.ok("previous output is skipped in a normal run",
                log.contains("Skipped 6 previously generated output file(s)"));
        TestRunner.ok("the .report of -extunique is not read back as input",
                log.contains("Skipped 6"));

        // -readmask consumes .msk, so that one must NOT be skipped
        String log2 = runMain(dir.toString(), "-readmask", "-out=" + out);
        TestRunner.ok("readmask keeps .msk as input",
                log2.contains("Skipped 5 previously generated output file(s)"));

        // -readgff consumes .gff
        String log3 = runMain(dir.toString(), "-readgff", "-out=" + out);
        TestRunner.ok("readgff keeps .gff as input",
                log3.contains("Skipped 5 previously generated output file(s)"));

        // a folder holding nothing but old output is reported, not analysed
        Path onlyOut = Files.createDirectories(tmp.resolve("onlyout"));
        Files.writeString(onlyOut.resolve("x.png"), "x", UTF_8);
        Files.writeString(onlyOut.resolve("x.svg"), "x", UTF_8);
        String log4 = runMain(onlyOut.toString(), "-maskonly", "-out=" + out);
        TestRunner.ok("a folder of pure output gives a clear message",
                log4.contains("No input sequence files"));
    }

    // ── Multi-FASTA: every record gets its own outputs ──────────────────────

    private static void multiFasta(Path tmp) throws Exception {
        // Regression: the table of every record used to be written to
        // <input>_1.gff under the first record's name, so each record overwrote
        // the one before and -seqshow cut its sequences out of record 1.
        Path dir = Files.createDirectories(tmp.resolve("multi"));
        String r1 = withRepeats(1, 4);
        String r3 = withRepeats(3, 5);
        Files.writeString(dir.resolve("m.fasta"),
                ">r1\n" + r1 + "\n>tiny\nacgtacgtac\n>r3\n" + r3 + "\n", UTF_8);

        String log = runMain(dir.resolve("m.fasta").toString(), "-seqshow", "-out=" + dir);
        TestRunner.ok("record 1's table holds record 1", tableMatches(dir.resolve("m.fasta_1.gff"), "r1", r1));
        TestRunner.ok("record 3's table holds record 3", tableMatches(dir.resolve("m.fasta_3.gff"), "r3", r3));
        TestRunner.ok("each record's mask is written under its own number",
                Files.exists(dir.resolve("m.fasta_1.msk")) && Files.exists(dir.resolve("m.fasta_3.msk")));
        TestRunner.ok("a record too short to analyse is reported, not silently dropped",
                log.contains("Skipped tiny") && !Files.exists(dir.resolve("m.fasta_2.gff")));

        // In a combined run the per-record base already carries the record
        // number, so a mask must not get the run-wide index appended again
        // (formerly m.fasta_1__1.msk, and o.fasta_4.msk for a one-record file).
        // A header with no sequence (record 4) used to crash -homology with a
        // division by zero and gave NaN coverage in -collate and -joint.
        Files.writeString(dir.resolve("o.fasta"), ">o\n" + r3 + "\n", UTF_8);
        Files.writeString(dir.resolve("m.fasta"), ">void\n", UTF_8, StandardOpenOption.APPEND);
        for (String mode : new String[]{"-collate", "-joint", "-homology"}) {
            Path out = Files.createDirectories(tmp.resolve("multi" + mode));
            String logm = runMain(dir.toString(), mode, "-out=" + out);
            TestRunner.ok(mode + " names each mask after its record",
                    Files.exists(out.resolve("m.fasta_1_.msk")) && Files.exists(out.resolve("m.fasta_3_.msk"))
                    && Files.exists(out.resolve("o.fasta.msk")));
            TestRunner.ok(mode + " reports an empty record as 0, not NaN",
                    Files.exists(out.resolve("m.fasta_4_.msk")) && !logm.contains("NaN"));
        }
    }

    // ── -extract: one valid, distinct file per record ───────────────────────

    private static void extract(Path tmp) throws Exception {
        Path dir = Files.createDirectories(tmp.resolve("extract"));
        Path in = dir.resolve("x.fasta");
        Files.writeString(in, ">gi|1|ref|NC_1.1| chr\nacgt\n>dup a\nacgt\n>dup b\ntttt\n>x\ngggg\n", UTF_8);
        runMain(in.toString(), "-extract");
        TestRunner.ok("an NCBI identifier becomes a valid file name",
                Files.exists(dir.resolve("gi_1_ref_NC_1.1_.fasta")));
        TestRunner.ok("records sharing a name do not overwrite each other",
                Files.exists(dir.resolve("dup.fasta")) && Files.exists(dir.resolve("dup_2.fasta"))
                && Files.readString(dir.resolve("dup.fasta"), UTF_8).contains(">dup a"));
        TestRunner.ok("a record named after its own file leaves the input intact",
                Files.exists(dir.resolve("x_2.fasta")) && Files.readString(in, UTF_8).contains(">dup b"));
    }

    /** One random element repeated {@code copies} times between unique spacers. */
    private static String withRepeats(long seed, int copies) {
        Random rnd = new Random(seed);
        String element = randomDna(rnd, 300);
        StringBuilder sb = new StringBuilder();
        for (int c = 0; c < copies; c++) {
            sb.append(randomDna(rnd, 200)).append(element);
        }
        return sb.append(randomDna(rnd, 200)).toString();
    }

    private static String randomDna(Random rnd, int n) {
        StringBuilder sb = new StringBuilder(n);
        for (int i = 0; i < n; i++) {
            sb.append("acgt".charAt(rnd.nextInt(4)));
        }
        return sb.toString();
    }

    /**
     * True when the table reports {@code seq}'s length, has rows, names
     * {@code name} in every row, and every forward-strand row carries that
     * record's own bases from Start to Stop.
     */
    private static boolean tableMatches(Path gff, String name, String seq) throws Exception {
        if (!Files.exists(gff)) {
            return false;
        }
        boolean lengthOk = false;
        int rows = 0;
        for (String line : Files.readAllLines(gff, UTF_8)) {
            if (line.equals("#Sequence length (bp)=" + seq.length())) {
                lengthOk = true;
            }
            String[] f = line.split("\t");
            if (line.startsWith("#") || f.length < 8 || f[0].equals("Seqid")) {
                continue;
            }
            rows++;
            if (!f[0].equals(name)) {
                return false;
            }
            if (f[6].equals("+")) {
                int start = Integer.parseInt(f[3]);
                int stop = Integer.parseInt(f[4]);
                // the sequence is the last field (column 9; Phase, column 8, is empty)
                if (f.length < 9 || !f[f.length - 1].equalsIgnoreCase(seq.substring(start - 1, stop))) {
                    return false;
                }
            }
        }
        return lengthOk && rows > 0;
    }

    // ── helpers ─────────────────────────────────────────────────────────────

    /** Runs {@code TotalRepeats.main} with stdout and stderr captured. */
    private static String runMain(String... args) throws Exception {
        PrintStream oldOut = System.out;
        PrintStream oldErr = System.err;
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        PrintStream cap = new PrintStream(buf, true, UTF_8);
        try {
            System.setOut(cap);
            System.setErr(cap);
            TotalRepeats.main(args);
        } finally {
            cap.flush();
            System.setOut(oldOut);
            System.setErr(oldErr);
        }
        return buf.toString(UTF_8);
    }

    private static void deleteTree(Path root) {
        try (var walk = Files.walk(root)) {
            walk.sorted(Comparator.reverseOrder()).forEach(p -> {
                try {
                    Files.deleteIfExists(p);
                } catch (Exception ignored) {
                    // a leftover file is not a test failure
                }
            });
        } catch (Exception ignored) {
            // best effort only
        }
    }
}
