/**
 * Block extraction from a soft-masked sequence.
 *
 * <p>Includes the regression test for the {@code 'Z'}/{@code 'z'} off-by-one:
 * {@code ReadMask} classified letters with {@code c > 64 && c < 90} and
 * {@code c > 96 && c < 122}, so the two boundary letters were not counted as
 * letters at all. Every base after such a letter then carried a sequence index
 * one too low, silently shifting all following block coordinates. The check
 * below compares a sequence containing {@code z}/{@code Z} against the identical
 * sequence with an ordinary letter in that position: the coordinates must match.
 */
final class MaskResultTest {

    private MaskResultTest() {
    }

    private static byte[] noSsr(int len) {
        return new byte[len];               // all zero: nothing excluded as STR
    }

    static void run() {
        TestRunner.suite("MaskResult");

        // ── one lowercase run in the middle ──
        // index:            0123 4567 89
        MaskResult m1 = new MaskResult();
        int[] b1 = m1.ReadMask("AAAAccccGG", 0, 1, noSsr(10));
        TestRunner.eq("single run: [start, length]", new int[]{4, 4}, b1);
        TestRunner.eq("repeatslen counts lowercase", 4, m1.getRepeatsLen());

        // ── run reaching the end of the sequence (the tail path) ──
        MaskResult m2 = new MaskResult();
        int[] b2 = m2.ReadMask("AAAAcccc", 0, 1, noSsr(8));
        TestRunner.eq("tail run", new int[]{4, 4}, b2);

        // ── run starting at position 0 ──
        MaskResult m3 = new MaskResult();
        int[] b3 = m3.ReadMask("ccccAAAA", 0, 1, noSsr(8));
        TestRunner.eq("leading run", new int[]{0, 4}, b3);

        // ── no lowercase at all ──
        MaskResult m4 = new MaskResult();
        int[] b4 = m4.ReadMask("AAAAAAAA", 0, 1, noSsr(8));
        TestRunner.eq("no masked region", new int[]{}, b4);

        // ── minLenSeq filters short runs ──
        MaskResult m5 = new MaskResult();
        int[] b5 = m5.ReadMask("AAccAA", 0, 5, noSsr(6));
        TestRunner.eq("run shorter than minLenSeq is dropped", new int[]{}, b5);

        MaskResult m6 = new MaskResult();
        int[] b6 = m6.ReadMask("AAccccccAA", 0, 6, noSsr(10));
        TestRunner.eq("run exactly at minLenSeq is kept", new int[]{2, 6}, b6);

        // ── N bases are counted as gaps ──
        MaskResult m7 = new MaskResult();
        m7.ReadMask("NNAAccccAA", 0, 1, noSsr(10));
        TestRunner.eq("gap bases counted", 2, m7.getGaps());

        // ── an STR-masked position cannot open a run ──
        byte[] ssr = noSsr(10);
        ssr[4] = 1;                          // first lowercase base is STR
        MaskResult m8 = new MaskResult();
        int[] b8 = m8.ReadMask("AAAAccccGG", 0, 1, ssr);
        TestRunner.ok("STR position does not open a run at that base",
                b8.length == 0 || b8[0] != 4);

        // ── regression: 'Z' and 'z' must count as letters ──
        zBoundaryRegression();
    }

    /**
     * 'z' (122) and 'Z' (90) sit exactly on the boundaries the old range test
     * excluded. Substituting an ordinary letter of the same case must not change
     * the reported coordinates; under the old code it shifted them by one.
     */
    private static void zBoundaryRegression() {
        // lowercase boundary: 'z' in the masked run, vs 'c' in the same place
        int[] withZ = new MaskResult().ReadMask("AAAAcczcGG", 0, 1, noSsr(10));
        int[] withC = new MaskResult().ReadMask("AAAAccccGG", 0, 1, noSsr(10));
        TestRunner.eq("lowercase 'z' counts as a letter", withC, withZ);

        // uppercase boundary: 'Z' in the flanking upper-case region, vs 'A'
        int[] withUpZ = new MaskResult().ReadMask("AAZAccccGG", 0, 1, noSsr(10));
        int[] withUpA = new MaskResult().ReadMask("AAAAccccGG", 0, 1, noSsr(10));
        TestRunner.eq("uppercase 'Z' counts as a letter", withUpA, withUpZ);

        // 'Z' immediately before the run is the position where a dropped letter
        // would shift the run start.
        int[] zBefore = new MaskResult().ReadMask("AAAZccccGG", 0, 1, noSsr(10));
        TestRunner.eq("'Z' just before a run does not shift it",
                new int[]{4, 4}, zBefore);

        // 'z' as the last base of a trailing run
        int[] zTail = new MaskResult().ReadMask("AAAAcccz", 0, 1, noSsr(8));
        TestRunner.eq("'z' closing a tail run", new int[]{4, 4}, zTail);

        // ReadUpperMask has always used inclusive ranges; the two must now agree
        // on which characters are letters.
        int[] upper = new MaskResult().ReadUpperMask("ccccAAAZ", 0, 1);
        TestRunner.eq("ReadUpperMask treats 'Z' as an upper-case letter",
                new int[]{4, 4}, upper);
    }
}
