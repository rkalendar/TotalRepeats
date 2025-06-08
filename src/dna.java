public final class dna {

    public static String AntisenseDNA(String source) {
        byte[] cdna = new byte[128];
        cdna[65] = 84;  //A
        cdna[66] = 86;  //B
        cdna[67] = 71;  //C
        cdna[68] = 72;  //D
        cdna[71] = 67;  //G
        cdna[72] = 68;  //H
        cdna[73] = 71;  //I
        cdna[75] = 77;  //K
        cdna[77] = 75;  //M
        cdna[78] = 78;  //N
        cdna[82] = 89;  //R
        cdna[83] = 83;  //S
        cdna[84] = 65;  //T
        cdna[85] = 65;  //U
        cdna[86] = 66;  //V
        cdna[87] = 87;  //W
        cdna[89] = 82;  //Y

        cdna[97] = 116;//  t <- a
        cdna[98] = 118;//  v <- b
        cdna[99] = 103;//  g <- c
        cdna[100] = 104;// h <- d
        cdna[103] = 99; // c <- g
        cdna[104] = 100;// d <- h
        cdna[105] = 99; // i <- g
        cdna[107] = 109;// m <- k
        cdna[109] = 107;// k <- m
        cdna[110] = 110;// n <- n
        cdna[114] = 121;// y <- r
        cdna[115] = 115;// s <- s
        cdna[116] = 97; // a <- t
        cdna[117] = 97; // a <- u
        cdna[118] = 98; // b <- v
        cdna[119] = 119;// w <- w
        cdna[121] = 114;// r <- y

        byte[] b = source.getBytes();
        for (int i = 0; i < source.length(); i++) {
            if (cdna[b[i]] > 0) {
                if (b[i] < 97) {
                    b[i] = (byte) (b[i] + 32);
                }
                b[i] = cdna[b[i]];
            }
        }
        return new String(b);
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
        byte[] b = source.getBytes();
        int l = source.length();
        int n = l / 2;
        for (int i = 0; i < n; i++) {
            if (tables.cdna[b[i]] > 0) {
                byte t = tables.cdna[b[l - i - 1]];
                b[l - i - 1] = tables.cdna[b[i]];
                b[i] = t;
            }
        }
        if ((l % 2) == 1) {
            if (tables.cdna[b[n]] > 0) {
                b[n] = tables.cdna[b[n]];
            }
        }
        return new String(b);
    }

    public static String DNA(byte[] b) {
        int l = b.length;
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

        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (cdn[b[i]] > 0) {
                    b[++n] = cdn[b[i]];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1)) : "";
    }

    public static String DNA(String source) {
        if (source == null || source.isEmpty()) {
            return "";
        }
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

        int l = source.length();
        byte[] b = source.getBytes();
        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (cdn[b[i]] > 0) {
                    b[++n] = cdn[b[i]];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1)) : "";
    }
}
