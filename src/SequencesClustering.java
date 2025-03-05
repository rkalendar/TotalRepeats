
import java.util.Arrays;
import java.util.HashMap;

public final class SequencesClustering {

    public SequencesClustering(String seq, int nkmer, int[] x, int similarity) {
        int nseq = x.length / 2;
        if (nseq < 1) {
            return;
        }
        if (similarity < 76) {
            similarity = 76;
        }
        if (similarity > 95) {
            similarity = 95;
        }

        d = new int[nseq][2];
        for (int j = 0; j < nseq; j++) {
            int p = j * 2;
            d[j][0] = x[p];
            d[j][1] = x[p + 1];
        }

        Arrays.sort(d, (int[] a, int[] b) -> {
            return Integer.compare(b[1], a[1]);
        });

        if (nkmer < 2) {
            ncl = QClustering(seq, nseq, similarity, nkmer);
        } else {
            cx = new int[nseq];      // cluster ID  
            int x1 = 2;
            for (int j = 0; j < 4; j++) {
                int g = Clustering2(seq, nseq, similarity, j * nkmer, nkmer * (1 + j), x1);
                if (g == x1) {
                    break;
                }
                x1 = g;
            }
            ncl = x1;
        }
    }

    private int QClustering(String seq, int nseq, int sim, int t) {
        String[] kmers = {"aatt", "acgt", "agct", "ttaa", "tgca", "tcga", "ccgg", "catg", "ctag", "ggcc", "gatc", "gtac"};

        int th = kmers.length;
        if (t == 1) {
            th = kmers.length;
        }
        if (t == 0) {
            th = 1;
        }

        HashMap<String, Integer> pt = new HashMap<>();
        for (int i = 0; i < th; i++) {
            pt.put(kmers[i], i);
        }

        int nkmers = pt.size();
        int kmer = 4;

        int[][] m2 = new int[nseq][nkmers];
        for (int j = 0; j < nseq; j++) {
            for (int i = 0; i < d[j][1] - kmer + 1; i++) {
                String s = seq.substring(d[j][0] + i, d[j][0] + i + kmer);
                if (pt.containsKey(s)) {
                    m2[j][pt.get(s)]++;
                }
            }
        }

        // compare sequences         
        int n = 2;               // amount of clusters
        cx = new int[nseq];      // cluster ID           
        for (int i = 0; i < nseq; i++) {

            if (cx[i] == 0) {
                n++;
                cx[i] = n;
                int sng = 0;

                for (int j = i + 1; j < nseq; j++) {
                    if (cx[j] == 0) {

                        int[] m = new int[nkmers + 1];
                        for (int k = 0; k < nkmers; k++) {
                            if (m2[i][k] > 1 && m2[j][k] > 1) {
                                if (m[0] == 0) {
                                    m[0] = 1;
                                    m[m[0]] = k;
                                } else {
                                    m[0]++;
                                    m[m[0]] = k;
                                }
                            }
                        }

                        int v = 0; // practice matches 
                        int z = 0; // theoretically maximum possible matches  
                        for (int k = 1; k < m[0]; k++) {
                            for (int y = k + 1; y < 1 + m[0]; y++) {

                                int di, dj;
                                z++;
                                if (m2[i][m[k]] < m2[i][m[y]]) {
                                    di = (100 * m2[i][m[k]]) / m2[i][m[y]];
                                    dj = (100 * m2[j][m[k]]) / m2[j][m[y]];
                                } else {
                                    di = (100 * m2[i][m[y]]) / m2[i][m[k]];
                                    dj = (100 * m2[j][m[y]]) / m2[j][m[k]];
                                }
                                if (di == dj) {
                                    v++;
                                }
                                if (di > dj) {
                                    if (di <= (dj * dif)) {
                                        v++;
                                    }
                                }
                                if (di < dj) {
                                    if (dj <= (di * dif)) {
                                        v++;
                                    }
                                }
                            }
                        }
                        if (v > 0 && ((100 * v) / z) > sim) {
                            cx[j] = n;
                            sng++;
                        }
                    }
                }
                if (sng == 0) {
                    n--;
                    cx[i] = 2;
                }
            }
        }
        return n;
    }

    private int Clustering2(String seq, int nseq, int sim, int strt, int size, int fn) {
        String[] kmers = {"tttc", "gaaa", "attt", "aaat", "tttg", "caaa", "tgaa", "ttct", "aaca", "ttca", "agaa", "tgtt", "tctt", "ttga", "aaga", "tcaa", "ttgt", "aaac", "gttt", "acaa", "aaag", "tgat", "cttt", "atca", "tcat", "atga", "catt", "acat", "attg", "caat", "atgt", "aatg", "attc", "gaat", "gatt", "aatc", "tatt", "aata", "atct", "ttcc", "agtt", "aact", "gaag", "cttc", "agat", "aagt", "ggaa", "actt", "ccaa", "tcca", "ttgg", "tgga", "ttgc", "tctc", "gcaa", "gatg", "catc", "tcct", "caac", "gttg", "agga", "caag", "tgta", "aagg", "taca", "agca", "tatc", "cctt", "gata", "tgct", "acca", "atgg", "tggt", "atac", "gtat", "cttg", "ccat", "cata", "gtga", "tatg", "tctg", "gctt", "ttat", "gcat", "gttc", "ggat", "agag", "ataa", "atcc", "aagc", "ggtt", "gaac", "caga", "aacc", "atgc", "ttac", "gtaa", "ctga", "tcag", "ctgt", "atag", "ctat", "agta", "aggt", "tact", "ctcc", "atta", "taat", "tgtc", "tcac", "tgag", "acct", "acag", "ctca", "ggag", "gtca", "agtg", "taga", "tcta", "cctc", "taaa", "tgac", "ttta", "gaca", "gagg", "actc", "gagt", "actg", "tggc", "cagt", "taac", "gtta", "gtgg", "ggga", "tagt", "ccac", "gctg", "acta", "cacc", "cagc", "tccc", "ctgc", "gcag", "ggta", "gtct", "tacc", "ccag", "cact", "gtag", "ctgg", "ggtg", "agac", "tgcc", "tggg", "agtc", "cagg", "gact", "ggca", "cctg", "aggg", "taag", "ctaa", "gcca", "tagc", "ccct", "ctta", "gctc", "gcta", "ttag", "gagc", "ggct", "ctac", "ttcg", "gcac", "tagg", "cgaa", "ccta", "gtgc", "gggt", "tcgt", "aggc", "cgtt", "acga", "gcct", "aacg", "atcg", "cgat", "gacc", "ggtc", "ggac", "accc", "ccca", "gtcc", "ctcg", "agcc", "ccga", "cgag", "acgg", "tcgc", "cgct", "cgga", "tccg", "cgtc", "cacg", "gacg", "tcgg", "agcg", "gcga", "cgcc", "gccc", "accg", "cgta", "gccg", "tacg", "cgtg", "ccgc", "cgac", "ggcg", "gggc", "cggt", "cggc", "gtcg", "ccgt", "cgca", "tgcg", "gcgg", "cggg", "acgc", "gcgt", "cccg"};
        /* 
        "aaat", "aaac", "aaag", "aata", "aaca", "aacc", "aaga", "aagg", "ataa", "atta", "attt", "acaa", "acca", "accc", "agaa", "agag", "agga", "aggg", "taaa", "taat", "tatt", "ttat", "ttta", "tttc", "tttg", "ttct", "ttcc", "ttgt", "ttgg", "tctt", "tctc", "tcct", "tccc", "tgtt", "tggt", "tggg", "caaa", "caac", "cacc", "cttt", "cttc", "ctcc", "ccaa", "ccac", "cctt", "cctc", "ccca", "ccct", "cccg", "ccgc", "cgcc", "cggc", "cggg", "gaaa", "gaag", "gagg", "gttt", "gttg", "gtgg", "gccc", "gccg", "gcgg", "ggaa", "ggag", "ggtt", "ggtg", "ggcg", "ggga", "gggt", "gggc",
        "ttca", "tgaa", "atct", "ctca", "caga", "tctg", "cctg", "aatg", "acag", "cagg", "ccag", "acat", "tgag", "ctgg", "tcat", "tcaa", "ctgt", "ttga", "tcca", "ctga", "agca", "atga", "catt", "tcag", "atgt", "actt", "cact", "aact", "attc", "tgga", "cagc", "aagt", "taca", "atca", "tgct", "tgta", "tgat", "gaat", "ccat", "gcag", "agtt", "gcct", "agat", "ctgc", "tcac", "cagt", "aggc", "gctg", "agtg", "caag", "atgg", "actg", "taga", "aatc", "agcc", "caat", "acct", "cttg", "gcca", "cata", "tcta", "gatt", "attg", "tggc", "gtga", "tgcc", "ggct", "ggca", "tatg", "agac", "ctaa", "aggt", "gcaa", "gaca", "actc", "aagc", "catc", "ttag", "tact", "gtct",
        "gctt", "agta", "ttgc", "ctta", "tgtc", "atag", "gatg", "ctat", "gcat", "ttac", "gagt", "taag", "gtaa", "atcc", "atgc", "acta", "ggat", "atac", "gtat", "tgac", "tatc", "tagt", "gagc", "gata", "gtca", "gaac", "gctc", "agtc", "gact", "ccta", "gttc", "taac", "gcac", "tagg", "ctac", "gtag", "gtgc", "gtta", "tagc", "gcta", "tacc", "ggac", "gacc", "ggta", "gtcc", "ggtc", "cacg", "cgtg", "ctcg", "cgag", "acgg", "cgtt", "aacg", "ccgt", "ttcg", "cgct", "agcg", "tcgt", "cgga", "tcgg", "ccga", "tccg", "cgca", "acgc", "acga", "cgaa", "gcgt", "tgcg", "gacg", "cgtc", "atcg", "cgat", "tcgc", "accg", "cggt", "gcga", "cgta", "tacg", "cgac", "gtcg"
         */
//   };

        if (size > kmers.length) {
            size = kmers.length;
        }
        int q = 0;
        HashMap<String, Integer> pt = new HashMap<>();
        for (int i = strt; i < size; i++) {
            pt.put(kmers[i], q++);
        }

        int nkmers = pt.size();
        int kmer = 4;

        int[][] m2 = new int[nseq][nkmers];
        int[][] m3 = new int[nseq][nkmers];

        for (int j = 0; j < nseq; j++) {
            if (cx[j] == 0) {
                for (int i = 0; i < d[j][1] - kmer + 1; i++) {
                    String s = seq.substring(d[j][0] + i, d[j][0] + i + kmer);
                    if (pt.containsKey(s)) {
                        m2[j][pt.get(s)]++;
                    }
                    s = dna.ComplementDNA2(s);
                    if (pt.containsKey(s)) {
                        m3[j][pt.get(s)]++;
                    }
                }
            }
        }

        // compare sequences         
        int n = fn;               // amount of clusters

        for (int i = 0; i < nseq; i++) {

            if (cx[i] == 0) {
                n++;
                cx[i] = n;
                int sng = 0;

                //Reverse direction
                for (int j = i + 1; j < nseq; j++) {
                    if (cx[j] == 0) {
                        int[] m = new int[nkmers + 1];
                        for (int k = 0; k < nkmers; k++) {
                            if (m2[i][k] > 1 && m3[j][k] > 1) {
                                if (m[0] == 0) {
                                    m[0] = 1;
                                    m[m[0]] = k;
                                } else {
                                    m[0]++;
                                    m[m[0]] = k;
                                }
                            }
                        }

                        int v = 0; // practice matches 
                        int z = 0; // theoretically maximum possible matches  
                        for (int k = 1; k < m[0]; k++) {
                            for (int y = k + 1; y < 1 + m[0]; y++) {

                                int di, dj;
                                z++;
                                if (m2[i][m[k]] < m2[i][m[y]]) {
                                    di = (100 * m2[i][m[k]]) / m2[i][m[y]];
                                    dj = (100 * m3[j][m[k]]) / m3[j][m[y]];
                                } else {
                                    di = (100 * m2[i][m[y]]) / m2[i][m[k]];
                                    dj = (100 * m3[j][m[y]]) / m3[j][m[k]];
                                }
                                if (di == dj) {
                                    v++;
                                }
                                if (di > dj) {
                                    if (di <= (dj * dif)) {
                                        v++;
                                    }
                                }
                                if (di < dj) {
                                    if (dj <= (di * dif)) {
                                        v++;
                                    }
                                }
                            }
                        }
                        if (v > 0 && ((100 * v) / z) > sim) {
                            cx[j] = -n;
                            sng++;
                        }
                    }
                }

                //Forward direction       
                for (int j = i + 1; j < nseq; j++) {
                    if (cx[j] == 0) {
                        int[] m = new int[nkmers + 1];
                        for (int k = 0; k < nkmers; k++) {
                            if (m2[i][k] > 1 && m2[j][k] > 1) {
                                if (m[0] == 0) {
                                    m[0] = 1;
                                    m[m[0]] = k;
                                } else {
                                    m[0]++;
                                    m[m[0]] = k;
                                }
                            }
                        }

                        int v = 0; // practice matches 
                        int z = 0; // theoretically maximum possible matches  
                        for (int k = 1; k < m[0]; k++) {
                            for (int y = k + 1; y < 1 + m[0]; y++) {

                                int di, dj;
                                z++;
                                if (m2[i][m[k]] < m2[i][m[y]]) {
                                    di = (100 * m2[i][m[k]]) / m2[i][m[y]];
                                    dj = (100 * m2[j][m[k]]) / m2[j][m[y]];
                                } else {
                                    di = (100 * m2[i][m[y]]) / m2[i][m[k]];
                                    dj = (100 * m2[j][m[y]]) / m2[j][m[k]];
                                }
                                if (di == dj) {
                                    v++;
                                }
                                if (di > dj) {
                                    if (di <= (dj * dif)) {
                                        v++;
                                    }
                                }
                                if (di < dj) {
                                    if (dj <= (di * dif)) {
                                        v++;
                                    }
                                }
                            }
                        }
                        if (v > 0 && ((100 * v) / z) > sim) {
                            cx[j] = n;
                            sng++;
                        }
                    }
                }

                if (sng == 0) {
                    n--;
                    cx[i] = 0;
                }
            }
        }
        return n;
    }

    public int getNcl() {
        return ncl;
    }

    public int[] Result() {
        return cx;
    }

    public int[][] ResultArray() {
        return d;
    }

    private final double dif = 1.2; //Deviation
    private int ncl;
    private int[] cx;
    private int[][] d;
}
