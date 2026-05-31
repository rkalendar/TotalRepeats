import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Detects and masks low-complexity / simple-sequence-repeat (SSR) regions
 * in a nucleotide sequence using a sliding di/trinucleotide diversity score.
 *
 * Sentinel values used in the working byte array:
 *   0–3  = ACGT (after dx2 remapping)
 *   4    = ambiguous / N base
 *   5    = masked (SSR) position
 */
public class LowComplexitySequence2 {

    // -----------------------------------------------------------------------
    // Constants
    // -----------------------------------------------------------------------

    /** Window length for the diversity score calculation. */
    private static final int WINDOW = 17;

    /** Minimum run length (bp) before a candidate region is masked. */
    private static final int SSR_MIN_RUN = 30;

    /** Minimum block length (bp) kept in the final SSR block list. */
    private static final int MIN_BLOCK_LEN = 50;

    /** Sentinel value written into the working array for masked bases. */
    private static final byte MASKED = 5;

    /** Sentinel value used for ambiguous (N) bases in the dx2 encoding. */
    private static final byte AMBIGUOUS = 4;

    // -----------------------------------------------------------------------
    // State
    // -----------------------------------------------------------------------

    private byte[] mapb   = new byte[0];
    private int[]  ibloks = new int[0];
    private int    totalssr = 0;

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /**
     * Scans {@code seq} for SSR / low-complexity blocks and populates the
     * internal state (map, block list, total length).
     *
     * @param seq          nucleotide sequence string
     * @param telomer      diversity threshold / gap-bridging length;
     *                     set to 0 to skip SSR detection
     * @param ssrDetection {@code false} forces telomer=0 (no masking)
     */
    public void FindAllSSRs(String seq, int telomer, boolean ssrDetection) {
        if (!ssrDetection) {
            telomer = 0;
        }

        // Work on a defensive copy so the caller's data is not mutated.
        final byte[] raw = seq.getBytes(StandardCharsets.US_ASCII);
        final byte[] b   = simpleRepeatsMasking(raw.clone(), telomer);

        mapb = new byte[b.length];
        final int n = b.length;
        final ArrayList<Integer> blocks = new ArrayList<>(128);

        int i = 0;
        while (i < n) {
            if (b[i] != MASKED) {
                i++;
                continue;
            }

            final int start = i;
            int e   = i + 1;
            int gap = 0;

            // Bridge short gaps between masked runs (up to telomer non-masked bases).
            while (e < n) {
                if (b[e] == MASKED) {
                    gap = 0;
                    e++;
                } else if (gap < telomer) {
                    gap++;
                    e++;
                } else {
                    break;
                }
            }

            final int len = e - start;
            if (len >= MIN_BLOCK_LEN) {
                blocks.add(start);
                blocks.add(len);
            }

            i = e;
        }

        // Store block pairs and build the map.
        ibloks = new int[blocks.size()];
        for (int k = 0; k < blocks.size(); k++) {
            ibloks[k] = blocks.get(k);
        }

        totalssr = 0;
        for (int k = 0; k + 1 < ibloks.length; k += 2) {
            final int s   = ibloks[k];
            final int len = ibloks[k + 1];
            totalssr += len;
            Arrays.fill(mapb, s, s + len, (byte) AMBIGUOUS);
        }
    }

    /** Returns the SSR blocks as a single-element list containing the flat int array. */
    public List<int[]> blocks() {
        return Collections.singletonList(ibloks.clone());
    }

    /** Returns a copy of the byte map (4 = SSR region, 0 = normal). */
    public byte[] MapBytes() {
        return mapb.clone();
    }

    /** Returns the flat array of (start, length) pairs for each SSR block. */
    public int[] IntBlocks() {
        return ibloks.clone();
    }

    /** Returns the total number of bases covered by detected SSR blocks. */
    public int GetTotalRepeats() {
        return totalssr;
    }

    // -----------------------------------------------------------------------
    // Internal masking logic
    // -----------------------------------------------------------------------

    /**
     * Masks low-complexity positions in {@code b} by writing {@link #MASKED}
     * into runs whose di/trinucleotide diversity score falls below
     * {@code threshold}.
     *
     * <p>The method modifies {@code b} in-place and returns it for convenience.
     * Callers should pass a copy if the original data must be preserved.
     *
     * @param b         working byte array (will be modified)
     * @param threshold minimum accepted diversity count (Kmax)
     * @return the same array {@code b}, with low-complexity bases set to 5
     */
    private byte[] simpleRepeatsMasking(byte[] b, int threshold) {
        final int l = b.length;
        if (l < WINDOW + 3) {
            return b;
        }

        // Using int[] avoids silent byte-overflow for frequency counts.
        final int[][] v2 = new int[5][5];
        final int[][][] v3 = new int[5][5][5];
        // k2[i] = number of distinct di+trinucleotides in [i, i+WINDOW).
        final int[] k2 = new int[l - WINDOW];

        // Remap raw ASCII bytes to 0-4 encoding via lookup table.
        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
        }

        // ---- Seed the first window ----------------------------------------
        for (int i = 0; i < WINDOW; i++) {
            final int a = b[i], c = b[i + 1];
            v2[a][c]++;
            if (v2[a][c] == 1) k2[0]++;
        }
        for (int i = 0; i < WINDOW - 1; i++) {
            final int a = b[i], c = b[i + 1], d = b[i + 2];
            v3[a][c][d]++;
            if (v3[a][c][d] == 1) k2[0]++;
        }

        // ---- Slide the window ---------------------------------------------
        for (int i = 1; i < l - WINDOW; i++) {
            k2[i] = k2[i - 1];

            // Remove the leftmost di- and trinucleotide.
            final int a = b[i - 1], c = b[i], d = b[i + 1];
            v2[a][c]--;
            if (v2[a][c] == 0) k2[i]--;
            v3[a][c][d]--;
            if (v3[a][c][d] == 0) k2[i]--;

            // Add the rightmost di- and trinucleotide.
            final int p = b[i + WINDOW - 2], q = b[i + WINDOW - 1], r = b[i + WINDOW];
            v2[q][r]++;
            if (v2[q][r] == 1) k2[i]++;
            v3[p][q][r]++;
            if (v3[p][q][r] == 1) k2[i]++;
        }

        // ---- Identify and mask low-diversity runs -------------------------
        for (int i = 0; i < l - WINDOW; i++) {
            if (k2[i] < threshold && b[i] < AMBIGUOUS) {
                final int x = i;
                int u = 0;
                int e = i + 1;

                while (e < l - WINDOW) {
                    if (k2[e] > threshold || b[e] == AMBIGUOUS) {
                        if (++u > WINDOW) break;
                    }
                    e++;
                }

                i = e - 1;   // advance outer loop past this candidate

                if (e - x > SSR_MIN_RUN) {
                    final int fillStart = x > 0 ? x - 1 : x;
                    // Guard against overflow at the end of the array.
                    final int fillEnd   = Math.min(e + WINDOW, l);
                    Arrays.fill(b, fillStart, fillEnd, MASKED);
                }
            }
        }

        return b;
    }
}
