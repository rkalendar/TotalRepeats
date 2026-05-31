import java.util.ArrayList;
import java.util.Arrays;

public class LowComplexitySequence {

    private byte[] SimpleRepeatsMasking(byte b[], int telomer) {
        int n1;
        int n2;
        int n3;
        int lmer = 17;
        int Kmax = telomer;   // Kmax=11 for SSR  telomers = 15;  
        int l = b.length;
        if (l < lmer + 3) {
            return b;
        }

        byte[][] v2 = new byte[5][5];
        byte[][][] v3 = new byte[5][5][5];
        byte[] k2 = new byte[l - lmer];

        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
        }

        for (int i = 0; i < lmer; i++) { //max=18
            n1 = b[i];        //16
            n2 = b[i + 1];    //17
            v2[n1][n2]++;
            if (v2[n1][n2] == 1) {
                k2[0]++;
            }
        }

        for (int i = 0; i < lmer - 1; i++) {
            n1 = b[i];       //15
            n2 = b[i + 1];   //16
            n3 = b[i + 2];   //17
            v3[n1][n2][n3]++;
            if (v3[n1][n2][n3] == 1) {
                k2[0]++;
            }
        }

        for (int i = 1; i < l - lmer; i++) {
            k2[i] = k2[i - 1];
            n1 = b[i - 1];
            n2 = b[i];
            n3 = b[i + 1];
            v2[n1][n2]--;
            if (v2[n1][n2] == 0) {
                k2[i]--;
            }

            v3[n1][n2][n3]--;
            if (v3[n1][n2][n3] == 0) {
                k2[i]--;
            }

            n1 = b[i + lmer - 2]; //16
            n2 = b[i + lmer - 1]; //17
            n3 = b[i + lmer];     //18
            v2[n2][n3]++;
            if (v2[n2][n3] == 1) {
                k2[i]++;
            }

            v3[n1][n2][n3]++;
            if (v3[n1][n2][n3] == 1) {
                k2[i]++;
            }
        }

        for (int i = 0; i < l - lmer; i++) {
            if (k2[i] < Kmax && b[i] < 4) {
                int x = i;
                int u = 0;
                int e = i + 1;
                while (e < l - lmer) {
                    if (k2[e] > Kmax || b[e] == 4) {
                        u++;
                        if (u > lmer) {
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
                    for (int h = x; h < e + lmer; h++) {
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
