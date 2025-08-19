import java.util.Arrays;

public final class ReadingSequencesFiles {

    public ReadingSequencesFiles(byte[] s) {
        LoadTable();
        ReadingSequences(s);
    }

    public ReadingSequencesFiles(String s) {
        LoadTable();
        ReadingSequences(s.getBytes());
    }

    private void LoadTable() {
        dnl = new byte[128];         //https://www.rapidtables.com/code/text/ascii-table.html
//M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T    
//a   b   c   d   g   h   i   k   m   n   r   s   t   u   v   w   y
//97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 121            
        dnl[65] = 97;   // A
        dnl[66] = 98;   // B
        dnl[67] = 99;   // C
        dnl[68] = 100;  // D
        dnl[71] = 103;  // G
        dnl[72] = 104;  // H
        dnl[73] = 99;   // I
        dnl[75] = 107;  // K
        dnl[77] = 109;  // M
        dnl[78] = 110;  // N
        dnl[82] = 114;  // R
        dnl[83] = 115;  // S
        dnl[84] = 116;  // T
        dnl[85] = 116;  // U
        dnl[86] = 118;  // V
        dnl[87] = 119;  // W
        dnl[89] = 121;  // Y        
        dnl[97] = 97;   // a
        dnl[98] = 98;    // b
        dnl[99] = 99;    // c
        dnl[100] = 100;  // d
        dnl[103] = 103;  // g
        dnl[104] = 104;  // h
        dnl[105] = 99;   // i
        dnl[107] = 107;  // k
        dnl[109] = 109;  // m
        dnl[110] = 110;  // n
        dnl[114] = 114;  // r
        dnl[115] = 115;  // s
        dnl[116] = 116;  // t
        dnl[117] = 116;  // u
        dnl[118] = 118;  // v
        dnl[119] = 119;  // w
        dnl[121] = 121;  // y 
    }

    public ReadingSequencesFiles() {
        LoadTable();
    }

    public String[] getSequences() {
        if (ns == 0) {
            return null;
        }
        return sequence;
    }

    public String[] getNames() {
        if (ns == 0) {
            return null;
        }
        return name_seq;
    }

    public int getNseq() {
        return ns;
    }

    public void ReadingMaskSequencesFiles(byte[] source) {
        if (source == null) {
            return;
        }
        int l = source.length;
        int s = 0; // total length
        ns = 0;    // amount fasta sequences

        for (int i = 0; i < l; i++) {
            if (source[i] < 9) {
                ns = 0;
                break;
            }
            if (dnl[source[i]] > 0) {
                s++;
            }
            if (source[i] == 62) {
                ns++;
            }
        }
        if (s == 0 || ns == 0) {
            return;
        }
        name_seq = new String[ns];
        sequence = new String[ns];
        int n = -1;
        int t = 0;

        for (int i = 0; i < l; i++) {
            if (source[i] == 62) {
                if (t > 0) {

                    byte[] d = Arrays.copyOfRange(source, t, i - 1); //public static short[] copyOfRange(short[] original, int from, int to)
                    int x = 0;
                    for (int j = 0; j < d.length; j++) {
                        if (dnl[d[j]] > 0) {
                            d[x++] = d[j];
                        }
                    }
                    sequence[n] = new String(Arrays.copyOfRange(d, 0, x));
                    lSeqs += x;
                }
                n++;
                for (int j = i + 1; j < l; j++) {
                    if (source[j] == 10 || source[j] == 13) {
                        name_seq[n] = new String(source, i + 1, j - i).trim();
                        i = j;
                        t = j + 1;
                        break;
                    }
                }
            }
        }
        byte[] d = Arrays.copyOfRange(source, t, l);
        int x = 0;
        for (int j = 0; j < d.length; j++) {
            if (dnl[d[j]] > 0) {
                d[x++] = d[j];
            }
        }

        sequence[n] = new String(Arrays.copyOfRange(d, 0, x));
        lSeqs += x;
    }

    private void ReadingSequences(byte[] source) {
        if (source == null) {
            return;
        }

        int l = source.length;
        int s = 0; // total length
        ns = 0;    // amount fasta sequences

        for (int i = 0; i < l; i++) {
            if (source[i] < 9) {
                ns = 0;
                break;
            }
            if (dnl[source[i]] > 0) {
                s++;
            }
            if (source[i] == 62) {
                ns++;
            }
        }
        if (s == 0 || ns == 0) {
            return;
        }
        name_seq = new String[ns];
        sequence = new String[ns];
        int n = -1;
        int t = 0;

        for (int i = 0; i < l; i++) {
            if (source[i] == 62) {
                if (t > 0) {
                    byte[] d = Arrays.copyOfRange(source, t, i - 1); //public static short[] copyOfRange(short[] original, int from, int to)
                    int x = 0;
                    for (int j = 0; j < d.length; j++) {
                        if (dnl[d[j]] > 0) {
                            d[x++] = dnl[d[j]];
                        }
                    }
                    sequence[n] = new String(Arrays.copyOfRange(d, 0, x));
                    lSeqs += x;
                }
                n++;
                for (int j = i + 1; j < l; j++) {
                    if (source[j] == 10 || source[j] == 13) {
                        name_seq[n] = new String(source, i + 1, j - i).trim();
                        i = j;
                        t = j + 1;
                        break;
                    }
                }
            }
        }
        byte[] d = Arrays.copyOfRange(source, t, l);
        int x = 0;
        for (int j = 0; j < d.length; j++) {
            if (dnl[d[j]] > 0) {
                d[x++] = dnl[d[j]];
            }
        }
        lSeqs = lSeqs + x;
        sequence[n] = new String(Arrays.copyOfRange(d, 0, x));
    }

    public int getLength() {
        return lSeqs;
    }

    private String[] name_seq;
    private String[] sequence;
    private int lSeqs = 0;
    private int ns = 0;
    private byte[] dnl;
}
