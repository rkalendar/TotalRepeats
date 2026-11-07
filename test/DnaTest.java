import static java.nio.charset.StandardCharsets.US_ASCII;

/**
 * Complement / reverse-complement / normalisation in {@link Dna}.
 *
 * <p>These lock down the refactor that hoisted the per-call 128-byte lookup
 * tables into static fields: a table that were accidentally mutated between
 * calls would show up as a second call returning something different from the
 * first, which {@code repeatedCallsAreStable} checks explicitly.
 */
final class DnaTest {

    private DnaTest() {
    }

    static void run() {
        TestRunner.suite("Dna");

        // ── AntisenseDNA: complement in place, no reversal, result lowercased ──
        TestRunner.eq("AntisenseDNA(acgt)", "tgca", Dna.AntisenseDNA("acgt"));
        TestRunner.eq("AntisenseDNA uppercase lowercases", "tgca", Dna.AntisenseDNA("ACGT"));

        // ── ComplementDNA: reverse complement (lowercased) ──
        TestRunner.eq("ComplementDNA(acgt)", "acgt", Dna.ComplementDNA("acgt"));
        TestRunner.eq("ComplementDNA(aaaa)", "tttt", Dna.ComplementDNA("aaaa"));
        TestRunner.eq("ComplementDNA(atgc)", "gcat", Dna.ComplementDNA("atgc"));

        // ── ComplementDNA2: reverse complement, original case preserved ──
        TestRunner.eq("ComplementDNA2 keeps case", "rynRYNacgtACGT",
                Dna.ComplementDNA2("ACGTacgtNRYnry"));
        TestRunner.eq("ComplementDNA2 is an involution", "acgtACGT",
                Dna.ComplementDNA2(Dna.ComplementDNA2("acgtACGT")));

        // odd length: the middle base must be complemented in place
        TestRunner.eq("ComplementDNA2 odd length", "gca", Dna.ComplementDNA2("tgc"));
        TestRunner.eq("ComplementDNA2 single base", "t", Dna.ComplementDNA2("a"));
        TestRunner.eq("ComplementDNA2 empty", "", Dna.ComplementDNA2(""));

        // ── Ambiguity codes complement to their partners ──
        TestRunner.eq("AntisenseDNA ambiguity r<->y", "yr", Dna.AntisenseDNA("ry"));
        TestRunner.eq("AntisenseDNA ambiguity k<->m", "mk", Dna.AntisenseDNA("km"));
        TestRunner.eq("AntisenseDNA self-complementary s/w/n", "swn", Dna.AntisenseDNA("swn"));

        // ── ReverseSeq ──
        TestRunner.eq("ReverseSeq(String)", "tgca", Dna.ReverseSeq("acgt"));
        TestRunner.eq("ReverseSeq(char[])", "tgca", Dna.ReverseSeq("acgt".toCharArray()));

        // ── DNA(): normalise to lowercase, drop non-nucleotide characters ──
        TestRunner.eq("DNA lowercases", "acgt", Dna.DNA("ACGT"));
        TestRunner.eq("DNA maps u->t", "t", Dna.DNA("U"));
        TestRunner.eq("DNA drops digits and spaces", "acgt", Dna.DNA("a1c 2g\t3t"));
        TestRunner.eq("DNA null is empty", "", Dna.DNA((String) null));
        TestRunner.eq("DNA empty is empty", "", Dna.DNA(""));
        TestRunner.eq("DNA of non-nucleotides is empty", "", Dna.DNA("123 456"));
        TestRunner.eq("DNA(byte[])", "acgt", Dna.DNA("ACGT".getBytes(US_ASCII)));

        // ── The static-table refactor must not leak state between calls ──
        repeatedCallsAreStable();
    }

    /**
     * The lookup tables are now static and shared. If any method mutated one in
     * place, a second call would differ from the first — including when calls to
     * different methods are interleaved.
     */
    private static void repeatedCallsAreStable() {
        String seq = "ACGTacgtNRYnrykmswbvdh";

        String anti1 = Dna.AntisenseDNA(seq);
        String comp1 = Dna.ComplementDNA(seq);
        String comp2a = Dna.ComplementDNA2(seq);
        String dna1 = Dna.DNA(seq);

        // interleave, then repeat every call
        Dna.DNA(seq);
        Dna.ComplementDNA2(seq);
        Dna.AntisenseDNA(seq);

        TestRunner.eq("AntisenseDNA stable across calls", anti1, Dna.AntisenseDNA(seq));
        TestRunner.eq("ComplementDNA stable across calls", comp1, Dna.ComplementDNA(seq));
        TestRunner.eq("ComplementDNA2 stable across calls", comp2a, Dna.ComplementDNA2(seq));
        TestRunner.eq("DNA stable across calls", dna1, Dna.DNA(seq));

        // DNA(byte[]) rewrites its argument in place; a fresh array must give the
        // same answer every time.
        String a = Dna.DNA(seq.getBytes(US_ASCII));
        String b = Dna.DNA(seq.getBytes(US_ASCII));
        TestRunner.eq("DNA(byte[]) stable across calls", a, b);
    }
}
