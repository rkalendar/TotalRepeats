import java.nio.charset.StandardCharsets;

public final class Dna {

    /**
     * DNA-normalisation map (uppercase/lowercase → lowercase, extended codes
     * E/F/J/L, inosine I→g), built ONCE and shared. Previously every DNA()/
     * AntisenseDNA() call rebuilt a 128-byte table; hoisting it to a static
     * removes that per-call allocation. AntisenseDNA reuses {@link Tables#cdna}
     * (byte-identical to the table it used to build); DNA keeps its own map
     * because it differs from {@link Tables#cdn} (e.g. I→g here vs I→c there).
     */
    private static final byte[] DNA_MAP = buildDnaMap();

    private static byte[] buildDnaMap() {
        byte[] cdn = new byte[128];
        cdn[65] = 97;   // A
        cdn[66] = 98;   // B
        cdn[67] = 99;   // C
        cdn[68] = 100;  // D
        cdn[69] = 97;   // dA=E
        cdn[70] = 99;   // dC=F
        cdn[71] = 103;  // G
        cdn[72] = 104;  // H
        cdn[73] = 103;  // I
        cdn[74] = 103;  // dG=J
        cdn[75] = 107;  // K
        cdn[76] = 116;  // dT=L
        cdn[77] = 109;  // M
        cdn[78] = 110;  // N
        cdn[82] = 114;  // R
        cdn[83] = 115;  // S
        cdn[84] = 116;  // T
        cdn[85] = 116;  // U
        cdn[86] = 118;  // V
        cdn[87] = 119;  // W
        cdn[89] = 121;  // Y

        cdn[97] = 97;    // a
        cdn[98] = 98;    // b
        cdn[99] = 99;    // c
        cdn[100] = 100;  // d
        cdn[101] = 97;   // dA=E
        cdn[102] = 99;   // dC=F
        cdn[103] = 103;  // g
        cdn[104] = 104;  // h
        cdn[105] = 103;  // i
        cdn[106] = 103;  // dG=J
        cdn[107] = 107;  // k
        cdn[108] = 116;  // dT=L
        cdn[109] = 109;  // m
        cdn[110] = 110;  // n
        cdn[114] = 114;  // r
        cdn[115] = 115;  // s
        cdn[116] = 116;  // t
        cdn[117] = 116;  // u
        cdn[118] = 118;  // v
        cdn[119] = 119;  // w
        cdn[121] = 121;  // y
        return cdn;
    }

    public static String AntisenseDNA(String source) {
        final byte[] cdna = Tables.cdna;   // shared static complement table
        byte[] b = source.getBytes(StandardCharsets.US_ASCII);
        for (int i = 0; i < b.length; i++) {
            if (cdna[b[i]] > 0) {
                if (b[i] < 97) {
                    b[i] = (byte) (b[i] + 32);
                }
                b[i] = cdna[b[i]];
            }
        }
        return new String(b, StandardCharsets.US_ASCII);
    }

    public static String ReverseSeq(String source) {
        return new StringBuilder(source).reverse().toString();
    }

    public static String ReverseSeq(char[] source) {
        StringBuilder s = new StringBuilder();
        return s.append(source).reverse().toString();
    }

    public static String ComplementDNA(String source) {
        return new StringBuilder(AntisenseDNA(source)).reverse().toString();
    }

    public static String ComplementDNA2(String source) {
        byte[] b = source.getBytes(StandardCharsets.US_ASCII);
        int l = source.length();
        int n = l / 2;
        for (int i = 0; i < n; i++) {
            if (Tables.cdna[b[i]] > 0) {
                byte t = Tables.cdna[b[l - i - 1]];
                b[l - i - 1] = Tables.cdna[b[i]];
                b[i] = t;
            }
        }
        if ((l % 2) == 1) {
            if (Tables.cdna[b[n]] > 0) {
                b[n] = Tables.cdna[b[n]];
            }
        }
        return new String(b, StandardCharsets.US_ASCII);
    }

    public static String DNA(byte[] b) {
        int l = b.length;
        final byte[] cdn = DNA_MAP;
        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (cdn[b[i]] > 0) {
                    b[++n] = cdn[b[i]];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1, StandardCharsets.US_ASCII)) : "";
    }

    public static String DNA(String source) {
        if (source == null || source.isEmpty()) {
            return "";
        }
        final byte[] cdn = DNA_MAP;
        int l = source.length();
        byte[] b = source.getBytes(StandardCharsets.US_ASCII);
        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (cdn[b[i]] > 0) {
                    b[++n] = cdn[b[i]];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1, StandardCharsets.US_ASCII)) : "";
    }

}
