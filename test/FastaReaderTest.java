import java.util.List;
import static java.nio.charset.StandardCharsets.US_ASCII;

/** FASTA parsing, IUPAC normalisation, and the soft-mask-preserving raw mode. */
final class FastaReaderTest {

    private FastaReaderTest() {
    }

    static void run() {
        TestRunner.suite("FastaReader");

        // ── single record, normalised ──
        FastaReader r = FastaReader.fromString(">seq1\nACGTacgt\n");
        TestRunner.eq("one record", 1, r.getSequenceCount());
        TestRunner.eq("name without '>'", "seq1", r.getNames().get(0));
        TestRunner.eq("normalised to lowercase", "acgtacgt", r.getSequences().get(0));
        TestRunner.eq("residue count", 8, r.getTotalResidues());
        TestRunner.ok("not empty", !r.isEmpty());

        // ── IUPAC equivalences applied by the normalising path ──
        FastaReader iupac = FastaReader.fromString(">s\nUuIi\n");
        TestRunner.eq("U/u -> t and I/i -> c", "ttcc", iupac.getSequences().get(0));

        FastaReader amb = FastaReader.fromString(">s\nNRYKMSWBVDH\n");
        TestRunner.eq("ambiguity codes kept, lowercased", "nrykmswbvdh",
                amb.getSequences().get(0));

        // ── characters outside the alphabet are dropped, sequence lines joined ──
        FastaReader dirty = FastaReader.fromString(">s\nac 12 gt\n\nac-gt\n");
        TestRunner.eq("invalid characters filtered out", "acgtacgt",
                dirty.getSequences().get(0));

        // ── multiple records ──
        FastaReader multi = FastaReader.fromString(">a desc here\nAC\nGT\n>b\nTTTT\n");
        TestRunner.eq("two records", 2, multi.getSequenceCount());
        TestRunner.eq("header keeps description", "a desc here", multi.getNames().get(0));
        TestRunner.eq("multi-line sequence joined", "acgt", multi.getSequences().get(0));
        TestRunner.eq("second record", "tttt", multi.getSequences().get(1));

        // ── raw mode preserves case, which is what carries soft-masking ──
        FastaReader raw = FastaReader.fromBytesRaw(">s\nACGTacgt\n".getBytes(US_ASCII));
        TestRunner.eq("raw mode preserves case", "ACGTacgt", raw.getSequences().get(0));

        FastaReader norm = FastaReader.fromBytes(">s\nACGTacgt\n".getBytes(US_ASCII));
        TestRunner.eq("byte mode normalises", "acgtacgt", norm.getSequences().get(0));

        // ── empty / headerless input ──
        FastaReader empty = FastaReader.fromString("");
        TestRunner.ok("empty input is empty", empty.isEmpty());
        TestRunner.eq("empty input has no records", 0, empty.getSequenceCount());

        FastaReader headerOnly = FastaReader.fromString(">lonely\n");
        TestRunner.eq("header with no sequence still yields a record", 1,
                headerOnly.getSequenceCount());
        TestRunner.eq("that record is empty", "", headerOnly.getSequences().get(0));

        // Regression: text before the first header used to be prepended to the
        // first record (here the "rtdrmt" among its letters).
        FastaReader preamble = FastaReader.fromString("Exported from tool, 2026\n>a\nACGT\n>b\nTT\n");
        TestRunner.eq("text before the first header is no record", 2, preamble.getSequenceCount());
        TestRunner.eq("nor part of the first one", "acgt", preamble.getSequences().get(0));
        TestRunner.eq("nor counted as residues", 6, preamble.getTotalResidues());

        // a plain sequence with no header at all is one record
        FastaReader plain = FastaReader.fromString("ACGTN\nacgt 60\n");
        TestRunner.eq("a header-less sequence is one record", 1, plain.getSequenceCount());
        TestRunner.eq("named by default", FastaReader.PLAIN_NAME, plain.getNames().get(0));
        TestRunner.eq("its residues", "acgtnacgt", plain.getSequences().get(0));

        // ... but prose is not, although many of its letters are IUPAC codes
        FastaReader prose = FastaReader.fromString("Run these commands on the cluster tonight\n");
        TestRunner.ok("prose is not read as a sequence", prose.isEmpty());
        TestRunner.eq("and has no residues", 0, prose.getTotalResidues());

        // ── accessors are unmodifiable, so callers cannot corrupt the reader ──
        List<String> seqs = r.getSequences();
        TestRunner.throwsEx("getSequences() is unmodifiable",
                UnsupportedOperationException.class, () -> seqs.add("x"));
        List<String> names = r.getNames();
        TestRunner.throwsEx("getNames() is unmodifiable",
                UnsupportedOperationException.class, () -> names.add("x"));
    }
}
