/**
 * Coordinate arithmetic of {@link SeqStore}, the virtual concatenation that lets
 * a combined run exceed the ~2.1 GB single-String limit.
 *
 * <p>This is the highest-risk arithmetic in the codebase: every global-to-local
 * conversion in the combined path goes through {@code segmentOf} and
 * {@code substring}, and an off-by-one at a segment boundary would silently
 * shift block coordinates rather than throw.
 */
final class SeqStoreTest {

    private SeqStoreTest() {
    }

    static void run() {
        TestRunner.suite("SeqStore");

        // segments:  "aaa" | "cccc" | "gg"      global: 0..8
        SeqStore s = new SeqStore(new String[]{"aaa", "cccc", "gg"});

        TestRunner.eq("length is the sum", 9, s.length());
        TestRunner.eq("segmentCount", 3, s.segmentCount());
        TestRunner.eq("segment(1)", "cccc", s.segment(1));

        // ── segmentStart ──
        TestRunner.eq("segmentStart(0)", 0, s.segmentStart(0));
        TestRunner.eq("segmentStart(1)", 3, s.segmentStart(1));
        TestRunner.eq("segmentStart(2)", 7, s.segmentStart(2));

        // ── segmentOf: every position, including both sides of each boundary ──
        TestRunner.eq("segmentOf(0)", 0, s.segmentOf(0));
        TestRunner.eq("segmentOf(2) last of seg0", 0, s.segmentOf(2));
        TestRunner.eq("segmentOf(3) first of seg1", 1, s.segmentOf(3));
        TestRunner.eq("segmentOf(6) last of seg1", 1, s.segmentOf(6));
        TestRunner.eq("segmentOf(7) first of seg2", 2, s.segmentOf(7));
        TestRunner.eq("segmentOf(8) last overall", 2, s.segmentOf(8));

        // ── charAt over the whole store ──
        StringBuilder all = new StringBuilder();
        for (long i = 0; i < s.length(); i++) {
            all.append(s.charAt(i));
        }
        TestRunner.eq("charAt walks the concatenation", "aaaccccgg", all.toString());

        // ── substring inside one segment (the fast path) ──
        TestRunner.eq("substring within seg0", "aa", s.substring(0, 2));
        TestRunner.eq("substring whole seg1", "cccc", s.substring(3, 7));
        TestRunner.eq("substring within seg1", "cc", s.substring(4, 6));
        TestRunner.eq("substring whole seg2", "gg", s.substring(7, 9));

        // ── substring across boundaries (the stitching path) ──
        TestRunner.eq("substring across one boundary", "aac", s.substring(1, 4));
        TestRunner.eq("substring across two boundaries", "accccg", s.substring(2, 8));
        TestRunner.eq("substring is the whole store", "aaaccccgg", s.substring(0, 9));

        // ── degenerate ranges ──
        TestRunner.eq("empty substring", "", s.substring(4, 4));
        TestRunner.eq("single char substring", "c", s.substring(3, 4));

        // ── a single segment behaves like a plain String ──
        SeqStore one = new SeqStore(new String[]{"acgtacgt"});
        TestRunner.eq("single segment length", 8, one.length());
        TestRunner.eq("single segment segmentOf", 0, one.segmentOf(7));
        TestRunner.eq("single segment substring", "gtac", one.substring(2, 6));

        // ── empty segments must not break the index ──
        SeqStore withEmpty = new SeqStore(new String[]{"ac", "", "gt"});
        TestRunner.eq("empty segment: length", 4, withEmpty.length());
        TestRunner.eq("empty segment: substring spans it", "acgt", withEmpty.substring(0, 4));
        TestRunner.eq("empty segment: charAt after it", 'g', withEmpty.charAt(2));
    }
}
