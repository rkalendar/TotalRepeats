import java.util.ArrayList;
import java.util.Arrays;

public class LowComplexitySequence {

    private byte[] SimpleRepeatsMasking(byte b[], int telomer) {
        final int lmer = 17;
        final int Kmax = telomer;   // Kmax=11 for SSR  telomers = 15;
        final int l = b.length;
        if (l < lmer + 3) {
            return b;
        }
        final int lim = l - lmer;

        // Hoist the static mapping table into a local so the JIT keeps the
        // array reference in a register and avoids re-reading the static field.
        final byte[] dx2 = tables.dx2;
        for (int i = 0; i < l; i++) {
            b[i] = dx2[b[i]];
        }

        // Flattened di-/tri-nucleotide count tables (were byte[5][5] and
        // byte[5][5][5]); a single 1D index avoids the jagged-array pointer
        // chasing and the extra bounds checks of nested array access.
        // index2 = n1*5 + n2 ; index3 = (n1*5 + n2)*5 + n3
        final int[] v2 = new int[25];
        final int[] v3 = new int[125];
        final byte[] k2 = new byte[lim];

        // First window: count distinct 2-mers and 3-mers. 'score' is kept in a
        // register and only written to k2 once per position.
        int score = 0;
        for (int i = 0; i < lmer; i++) {
            int e2 = b[i] * 5 + b[i + 1];
            if (++v2[e2] == 1) {
                score++;
            }
        }
        for (int i = 0; i < lmer - 1; i++) {
            int e3 = (b[i] * 5 + b[i + 1]) * 5 + b[i + 2];
            if (++v3[e3] == 1) {
                score++;
            }
        }
        k2[0] = (byte) score;

        // Slide the window: drop the leftmost 2-/3-mer, add the rightmost.
        for (int i = 1; i < lim; i++) {
            int n1 = b[i - 1], n2 = b[i], n3 = b[i + 1];
            int r2 = n1 * 5 + n2;
            if (--v2[r2] == 0) {
                score--;
            }
            int r3 = r2 * 5 + n3;
            if (--v3[r3] == 0) {
                score--;
            }

            int m1 = b[i + lmer - 2], m2 = b[i + lmer - 1], m3 = b[i + lmer];
            int a2 = m2 * 5 + m3;
            if (++v2[a2] == 1) {
                score++;
            }
            int a3 = (m1 * 5 + m2) * 5 + m3;
            if (++v3[a3] == 1) {
                score++;
            }

            k2[i] = (byte) score;
        }

        // Mark low-complexity stretches as 5.
        for (int i = 0; i < lim; i++) {
            if (k2[i] < Kmax && b[i] < 4) {
                int x = i;
                int u = 0;
                int e = i + 1;
                while (e < lim) {
                    if (k2[e] > Kmax || b[e] == 4) {
                        if (++u > lmer) {
                            break;
                        }
                    }
                    e++;
                }
                i = e - 1;
                if (e - x > ssrlen) {
                    if (x > 0) {
                        x--;
                    }
                    final int end = e + lmer;
                    for (int h = x; h < end; h++) {
                        b[h] = 5;
                    }
                }
            }
        }
        return b;
    }

    public void FindAllSSRs(String seq, int telomer, boolean SSRdetection) {
        if (!SSRdetection) {
            telomer = 0;
        }
        final byte[] raw = seq.getBytes();
        final byte[] b = SimpleRepeatsMasking(raw, telomer);
        mapb = new byte[b.length];

        final int n = b.length;
        final ArrayList<Integer> blocks = new ArrayList<>(128);

        int i = 0;
        while (i < n) {
            if (b[i] != 5) { // not SSR start
                i++;
                continue;
            }

            int start = i;
            int e = i + 1;

            // cumulative gap bridging: allow up to 'telomer' consecutive non-5s
            int gap = 0;
            while (e < n) {
                if (b[e] == 5) {
                    gap = 0;      // reset gap inside SSR island
                    e++;
                } else if (gap < telomer) {
                    gap++;        // consume gap byte
                    e++;
                } else {
                    break;        // exceeded allowed gap
                }
            }

            int len = e - start;
            if (len >= minlenseq) {
                blocks.add(start);
                blocks.add(len);
            }

            i = e; // jump past the merged block
        }

        // Materialize ibloks and fill map
        ibloks = new int[blocks.size()];
        for (int k = 0; k < blocks.size(); k++) {
            ibloks[k] = blocks.get(k);
        }

        totalssr = 0;
        for (int k = 0; k < ibloks.length; k += 2) {
            int s = ibloks[k];
            int len = ibloks[k + 1];
            totalssr += len;
            Arrays.fill(mapb, s, s + len, (byte) 4); // mark merged block
        }
    }

    public ArrayList Blocks() {
        ArrayList<int[]> bb = new ArrayList<>();
        bb.add(ibloks);
        return bb;
    }

    public byte[] MapBytes() {
        return mapb;
    }

    public int[] IntBlocks() {
        return ibloks;
    }

    public int GetTotalRepeats() {
        return totalssr;
    }
    private byte[] mapb;
    private int[] ibloks;
    private int totalssr = 0;
    private final int ssrlen = 30;
    private final int minlenseq = 50;     // sequence length 
}
