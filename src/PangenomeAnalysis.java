import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Pangenome analysis of the repeat families produced by a combined
 * (multi-file) TotalRepeats run.
 *
 * <p>The combined clustering groups homologous repeat instances detected across
 * the whole input set into families (CRP clusters). Every instance is mapped
 * back to the sequence record it came from (via the cumulative record-boundary
 * array {@code seqslen}) and from there to its <b>genome</b>. A genome is, by
 * default, one input file — all records of a multi-record FASTA (chromosomes,
 * contigs, plasmids) count as one genome — or, on request, every record on its
 * own. Each family then has a copy-number profile over the genomes, from which
 * the analysis derives:
 *
 * <ul>
 *   <li>the classic <b>core / accessory / unique</b> split (present in all /
 *       some / exactly one genome) and the four-tier classification used by
 *       bacterial pangenome tools (core ≥ 99 %, soft-core 95–99 %, shell
 *       15–95 %, cloud &lt; 15 % of genomes);</li>
 *   <li>the <b>family frequency spectrum</b> (the pangenome U-curve);</li>
 *   <li><b>rarefaction curves</b> — the expected pangenome and core sizes as
 *       genomes are added — computed <i>exactly</i> over all orderings of the
 *       genomes (no random permutations, so the result is reproducible), and a
 *       <b>Heaps' law</b> fit of the number of new families per added genome,
 *       which tells an open pangenome from a closed one;</li>
 *   <li><b>pairwise genome distances</b> — Jaccard on presence/absence (all
 *       families and accessory-only) and Bray–Curtis on copy numbers — and a
 *       <b>UPGMA tree</b> of the genomes in Newick format;</li>
 *   <li><b>presence patterns</b>: which combinations of genomes share which
 *       families (the input of an UpSet plot);</li>
 *   <li><b>copy-number variation</b>: families whose copy number varies most
 *       between genomes, and per-genome family expansions.</li>
 * </ul>
 *
 * <p>The class is self-contained: it needs only the cluster blocks, the record
 * boundaries, the names and (optionally) the record-to-genome map, so every
 * combined run mode can use it. All numbers are formatted with
 * {@link Locale#ROOT}, so the decimal separator is a point whatever the system
 * locale is — the files are meant to be parsed.
 */
public final class PangenomeAnalysis {

    // ── Tier thresholds (fractions of the genome count), Roary/Panaroo convention ──
    static final double CORE_FRACTION = 0.99;
    static final double SOFT_CORE_FRACTION = 0.95;
    static final double SHELL_FRACTION = 0.15;

    /** Rows listed per section in the text report (the TSV files are complete). */
    private static final int TOP_PATTERNS = 50;
    private static final int TOP_CNV = 20;

    /** One repeat family (CRP cluster) with its per-genome copy-number profile. */
    public static final class Family {

        public final int clusterId;        // matches the ClusterID column of the GFF report
        public final String refLabel;      // reference annotation, or "-"
        public final int[] count;          // instances per genome
        public final long[] bp;            // covered bp per genome
        public int frequency;              // number of genomes in which it is present
        public int totalMembers;
        public long totalBp;
        public String category;            // core | accessory | unique
        public String tier;                // core | soft-core | shell | cloud
        public String pattern;             // presence string over the genomes, e.g. "1101"
        public int minCopies, maxCopies;   // over all genomes (absent = 0)
        public double meanCopies, cv;      // mean and coefficient of variation over all genomes

        Family(int clusterId, String refLabel, int nGenomes) {
            this.clusterId = clusterId;
            this.refLabel = refLabel;
            this.count = new int[nGenomes];
            this.bp = new long[nGenomes];
        }
    }

    /** A set of genomes sharing exactly the same families: one row of an UpSet plot. */
    public static final class PatternGroup {

        public final String pattern;
        public final int frequency;
        public int families;
        public long bp;

        PatternGroup(String pattern, int frequency) {
            this.pattern = pattern;
            this.frequency = frequency;
        }
    }

    /** A node of the UPGMA genome tree; a leaf has {@code leaf >= 0} and no children. */
    public static final class Node {

        public final int leaf;             // genome index, or -1 for an internal node
        public final Node left, right;
        public final double height;        // distance from the leaves (UPGMA: half the merge distance)
        public final int size;             // leaves below this node

        Node(int leaf) {
            this.leaf = leaf;
            this.left = null;
            this.right = null;
            this.height = 0;
            this.size = 1;
        }

        Node(Node left, Node right, double height) {
            this.leaf = -1;
            this.left = left;
            this.right = right;
            this.height = height;
            this.size = left.size + right.size;
        }
    }

    // ── Input ────────────────────────────────────────────────────────────────
    private final int nGenomes;
    private final String[] names;
    private final long[] seqslen;          // cumulative end coordinates per record (global, long)
    private final int[] recordGenome;      // record index -> genome index
    private final String[] recordNames;
    private final boolean byFile;          // one genome = one input file (else one record)
    private final long[] genomeLength;     // bp per genome (sum of its records)

    // ── Families and aggregates ──────────────────────────────────────────────
    private final List<Family> families = new ArrayList<>();
    private final int[] ucrpCount;         // genome-specific unclustered repeats
    private final long[] ucrpBp;
    private final int[] strCount;          // short tandem repeats
    private final long[] strBp;

    private int coreCount, accessoryCount, uniqueCount;
    private long coreBp, accessoryBp, uniqueBp;
    private final int[] tierCount = new int[4];
    private final long[] tierBp = new long[4];
    private static final String[] TIER_NAMES = {"core", "soft-core", "shell", "cloud"};

    // ── Derived statistics ──────────────────────────────────────────────────
    private int[] specCount;               // families present in exactly k genomes (index k)
    private long[] specBp;
    private double[] expPan, expCore, expNew;   // rarefaction, index n = 1..N (0 unused)
    private double heapsAlpha = Double.NaN, heapsKappa = Double.NaN, heapsGamma = Double.NaN;
    private int[][] shared, union;
    private double[][] jaccard, jaccardAccessory, brayCurtis;
    private Node tree;
    private int[] leafOrder;
    private List<PatternGroup> patterns;

    // ── Construction ────────────────────────────────────────────────────────
    /** Every sequence record is its own genome (the original behaviour). */
    public PangenomeAnalysis(ArrayList<long[]> bb, long[] seqslen, String[] names,
                             int[] refclust, String[] refsname) {
        this(bb, seqslen, names, refclust, refsname, null, null);
    }

    /**
     * @param recordGenome genome index of every record, or {@code null} to treat
     *                     every record as its own genome
     * @param genomeNames  one name per genome when {@code recordGenome} is given
     */
    public PangenomeAnalysis(ArrayList<long[]> bb, long[] seqslen, String[] names,
                             int[] refclust, String[] refsname,
                             int[] recordGenome, String[] genomeNames) {
        int nRecords = seqslen.length;
        this.seqslen = seqslen;
        this.recordNames = uniqueNames(names);
        if (recordGenome == null || genomeNames == null) {
            this.recordGenome = new int[nRecords];
            for (int i = 0; i < nRecords; i++) {
                this.recordGenome[i] = i;
            }
            this.names = this.recordNames;
            this.byFile = false;
        } else {
            if (recordGenome.length < nRecords) {
                throw new IllegalArgumentException("recordGenome covers " + recordGenome.length
                        + " records but there are " + nRecords);
            }
            this.recordGenome = recordGenome;
            this.names = uniqueNames(genomeNames);
            this.byFile = true;
        }
        this.nGenomes = this.names.length;

        this.genomeLength = new long[nGenomes];
        long prev = 0;
        for (int i = 0; i < nRecords; i++) {
            genomeLength[this.recordGenome[i]] += seqslen[i] - prev;
            prev = seqslen[i];
        }

        this.ucrpCount = new int[nGenomes];
        this.ucrpBp = new long[nGenomes];
        this.strCount = new int[nGenomes];
        this.strBp = new long[nGenomes];

        // bb index 0 = STR, index 1 = UCRP, index >= 2 = CRP families.
        for (int i = 0; i < bb.size(); i++) {
            long[] z = bb.get(i);
            if (z == null) {
                continue;
            }
            int clusterId = i + 1;       // same numbering as the GFF ClusterID column

            if (i == 0) {                // STR aggregate
                for (int j = 0; j + 1 < z.length; j += 2) {
                    int g = genomeOf(z[j]);
                    strCount[g]++;
                    strBp[g] += Math.abs(z[j + 1]);
                }
                continue;
            }
            if (i == 1) {                // UCRP: unclustered repeats
                for (int j = 0; j + 1 < z.length; j += 2) {
                    int g = genomeOf(z[j]);
                    ucrpCount[g]++;
                    ucrpBp[g] += Math.abs(z[j + 1]);
                }
                continue;
            }

            // CRP homology family
            String refLabel = "-";
            if (refclust != null && clusterId < refclust.length && refclust[clusterId] > 0
                    && refsname != null && refclust[clusterId] - 1 < refsname.length) {
                refLabel = refclust[clusterId] + ":" + refsname[refclust[clusterId] - 1];
            }
            Family f = new Family(clusterId, refLabel, nGenomes);
            for (int j = 0; j + 1 < z.length; j += 2) {
                int g = genomeOf(z[j]);
                f.count[g]++;
                f.bp[g] += Math.abs(z[j + 1]);
            }
            finishFamily(f);
            if (f.frequency == 0) {
                continue;               // defensive: empty family
            }
            classify(f);
            families.add(f);
        }

        computeSpectrum();
        computeRarefaction();
        computeDistances();
        computeTree();
        computePatterns();
    }

    /** Appends an index to any name that occurs more than once, so columns stay distinguishable. */
    private static String[] uniqueNames(String[] in) {
        String[] out = new String[in.length];
        Map<String, Integer> seen = new HashMap<>();
        for (String s : in) {
            seen.merge(s == null ? "" : s, 1, Integer::sum);
        }
        for (int i = 0; i < in.length; i++) {
            String key = (in[i] == null) ? "" : in[i];
            String s = key.isEmpty() ? "seq" : key;
            out[i] = (seen.get(key) > 1) ? s + "#" + (i + 1) : s;
        }
        return out;
    }

    /** Maps a global start coordinate to its genome index (binary search over the record ends). */
    private int genomeOf(long start) {
        int lo = 0, hi = seqslen.length - 1;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (seqslen[mid] > start) {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        return recordGenome[lo];        // lo = first record whose end exceeds start (or the last one)
    }

    private void finishFamily(Family f) {
        StringBuilder pat = new StringBuilder(nGenomes);
        f.minCopies = Integer.MAX_VALUE;
        double sum = 0, sumSq = 0;
        for (int g = 0; g < nGenomes; g++) {
            int c = f.count[g];
            if (c > 0) {
                f.frequency++;
            }
            pat.append(c > 0 ? '1' : '0');
            f.totalMembers += c;
            f.totalBp += f.bp[g];
            f.minCopies = Math.min(f.minCopies, c);
            f.maxCopies = Math.max(f.maxCopies, c);
            sum += c;
            sumSq += (double) c * c;
        }
        f.pattern = pat.toString();
        f.meanCopies = sum / nGenomes;
        double var = sumSq / nGenomes - f.meanCopies * f.meanCopies;
        f.cv = (f.meanCopies > 0 && var > 0) ? Math.sqrt(var) / f.meanCopies : 0.0;
    }

    private void classify(Family f) {
        if (f.frequency == nGenomes) {
            f.category = "core";
            coreCount++;
            coreBp += f.totalBp;
        } else if (f.frequency == 1) {
            f.category = "unique";
            uniqueCount++;
            uniqueBp += f.totalBp;
        } else {
            f.category = "accessory";
            accessoryCount++;
            accessoryBp += f.totalBp;
        }
        int t = tierOf(f.frequency, nGenomes);
        f.tier = TIER_NAMES[t];
        tierCount[t]++;
        tierBp[t] += f.totalBp;
    }

    /** 0 = core, 1 = soft-core, 2 = shell, 3 = cloud. */
    static int tierOf(int frequency, int nGenomes) {
        double frac = (double) frequency / nGenomes;
        if (frac >= CORE_FRACTION) {
            return 0;
        }
        if (frac >= SOFT_CORE_FRACTION) {
            return 1;
        }
        if (frac >= SHELL_FRACTION) {
            return 2;
        }
        return 3;
    }

    // ── Derived statistics ──────────────────────────────────────────────────
    private void computeSpectrum() {
        specCount = new int[nGenomes + 1];
        specBp = new long[nGenomes + 1];
        for (Family f : families) {
            specCount[f.frequency]++;
            specBp[f.frequency] += f.totalBp;
        }
    }

    /**
     * Expected pangenome / core size when {@code n} of the {@code N} genomes are
     * sampled uniformly at random, computed exactly from the frequency spectrum:
     * a family present in {@code f} genomes is missed by a sample of {@code n}
     * with probability C(N−f, n)/C(N, n) and present in every sampled genome
     * with probability C(f, n)/C(N, n). Summing over the spectrum gives the
     * average over all N! orderings without enumerating or sampling them.
     */
    private void computeRarefaction() {
        int N = nGenomes;
        expPan = new double[N + 1];
        expCore = new double[N + 1];
        expNew = new double[N + 1];
        for (int n = 1; n <= N; n++) {
            double pan = 0, core = 0;
            for (int f = 1; f <= N; f++) {
                int cnt = specCount[f];
                if (cnt == 0) {
                    continue;
                }
                double pAbsent = 1.0;       // C(N-f, n) / C(N, n)
                for (int i = 0; i < n && pAbsent > 0; i++) {
                    pAbsent *= (double) (N - f - i) / (N - i);
                }
                double pAll = 1.0;          // C(f, n) / C(N, n)
                for (int i = 0; i < n && pAll > 0; i++) {
                    pAll *= (double) (f - i) / (N - i);
                }
                pan += cnt * (1.0 - pAbsent);
                core += cnt * pAll;
            }
            expPan[n] = pan;
            expCore[n] = core;
            expNew[n] = pan - expPan[n - 1];
        }

        // Heaps' law (Tettelin et al. 2008): new(n) = kappa * n^-alpha, fitted on n = 2..N.
        double[] x = new double[N + 1];
        for (int n = 0; n <= N; n++) {
            x[n] = n;
        }
        double[] fitNew = powerFit(x, expNew, 2, N);
        if (fitNew != null) {
            heapsAlpha = -fitNew[0];
            heapsKappa = fitNew[1];
        }
        double[] fitPan = powerFit(x, expPan, 1, N);   // pan(n) = kappa' * n^gamma
        if (fitPan != null) {
            heapsGamma = fitPan[0];
        }
    }

    /** Least-squares fit of y = k·x^e in log–log space over x[from..to]; {e, k}, or null if underdetermined. */
    static double[] powerFit(double[] x, double[] y, int from, int to) {
        int n = 0;
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = from; i <= to; i++) {
            if (y[i] <= 0 || x[i] <= 0) {
                continue;
            }
            double lx = Math.log(x[i]), ly = Math.log(y[i]);
            n++;
            sx += lx;
            sy += ly;
            sxx += lx * lx;
            sxy += lx * ly;
        }
        if (n < 2) {
            return null;
        }
        double den = n * sxx - sx * sx;
        if (den <= 1e-12) {
            return null;
        }
        double slope = (n * sxy - sx * sy) / den;
        double icpt = (sy - slope * sx) / n;
        return new double[]{slope, Math.exp(icpt)};
    }

    private void computeDistances() {
        int N = nGenomes;
        int F = families.size();
        int words = (F + 63) >>> 6;
        long[][] bits = new long[N][words];
        long[] accessoryMask = new long[words];     // families that are not core
        for (int i = 0; i < F; i++) {
            Family f = families.get(i);
            for (int g = 0; g < N; g++) {
                if (f.count[g] > 0) {
                    bits[g][i >>> 6] |= 1L << i;
                }
            }
            if (!"core".equals(f.category)) {
                accessoryMask[i >>> 6] |= 1L << i;
            }
        }

        shared = new int[N][N];
        union = new int[N][N];
        jaccard = new double[N][N];
        jaccardAccessory = new double[N][N];
        brayCurtis = new double[N][N];
        for (int a = 0; a < N; a++) {
            jaccard[a][a] = 1.0;
            jaccardAccessory[a][a] = 1.0;
            shared[a][a] = popcount(bits[a]);
            union[a][a] = shared[a][a];
            for (int b = a + 1; b < N; b++) {
                int s = 0, u = 0, sa = 0, ua = 0;
                for (int w = 0; w < words; w++) {
                    long and = bits[a][w] & bits[b][w];
                    long or = bits[a][w] | bits[b][w];
                    s += Long.bitCount(and);
                    u += Long.bitCount(or);
                    sa += Long.bitCount(and & accessoryMask[w]);
                    ua += Long.bitCount(or & accessoryMask[w]);
                }
                shared[a][b] = shared[b][a] = s;
                union[a][b] = union[b][a] = u;
                // Two genomes with no families at all are identical: J = 1, not 0.
                jaccard[a][b] = jaccard[b][a] = (u > 0) ? (double) s / u : 1.0;
                jaccardAccessory[a][b] = jaccardAccessory[b][a] = (ua > 0) ? (double) sa / ua : 1.0;

                long diff = 0, sum = 0;
                for (Family f : families) {
                    diff += Math.abs(f.count[a] - f.count[b]);
                    sum += f.count[a] + f.count[b];
                }
                brayCurtis[a][b] = brayCurtis[b][a] = (sum > 0) ? (double) diff / sum : 0.0;
            }
        }
    }

    private static int popcount(long[] w) {
        int c = 0;
        for (long x : w) {
            c += Long.bitCount(x);
        }
        return c;
    }

    private void computeTree() {
        double[][] d = new double[nGenomes][nGenomes];
        for (int a = 0; a < nGenomes; a++) {
            for (int b = 0; b < nGenomes; b++) {
                d[a][b] = 1.0 - jaccard[a][b];
            }
        }
        tree = upgma(d);
        leafOrder = new int[nGenomes];
        collectLeaves(tree, leafOrder, new int[1]);
    }

    /**
     * UPGMA (average linkage) on a symmetric distance matrix. Ties are broken by
     * the lowest pair of indices, so the tree is deterministic. Node heights are
     * half the merge distance; branch lengths in the Newick output are height
     * differences, so the tree is ultrametric.
     */
    static Node upgma(double[][] dist) {
        int n = dist.length;
        List<Node> nodes = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            nodes.add(new Node(i));
        }
        double[][] d = new double[n][];
        for (int i = 0; i < n; i++) {
            d[i] = dist[i].clone();
        }
        boolean[] dead = new boolean[n];
        int alive = n;
        while (alive > 1) {
            int bi = -1, bj = -1;
            double best = Double.POSITIVE_INFINITY;
            for (int i = 0; i < n; i++) {
                if (dead[i]) {
                    continue;
                }
                for (int j = i + 1; j < n; j++) {
                    if (!dead[j] && d[i][j] < best) {
                        best = d[i][j];
                        bi = i;
                        bj = j;
                    }
                }
            }
            Node a = nodes.get(bi), b = nodes.get(bj);
            Node merged = new Node(a, b, best / 2.0);
            for (int k = 0; k < n; k++) {
                if (dead[k] || k == bi || k == bj) {
                    continue;
                }
                double nd = (a.size * d[bi][k] + b.size * d[bj][k]) / (a.size + b.size);
                d[bi][k] = d[k][bi] = nd;
            }
            nodes.set(bi, merged);
            dead[bj] = true;
            alive--;
        }
        for (int i = 0; i < n; i++) {
            if (!dead[i]) {
                return nodes.get(i);
            }
        }
        return null;
    }

    private static void collectLeaves(Node t, int[] out, int[] pos) {
        if (t == null) {
            return;
        }
        if (t.leaf >= 0) {
            out[pos[0]++] = t.leaf;
            return;
        }
        collectLeaves(t.left, out, pos);
        collectLeaves(t.right, out, pos);
    }

    /** Newick serialisation with branch lengths; names are made Newick-safe. */
    static String newick(Node t, String[] names) {
        StringBuilder sb = new StringBuilder();
        appendNewick(t, names, sb);
        return sb.append(';').toString();
    }

    private static void appendNewick(Node t, String[] names, StringBuilder sb) {
        if (t.leaf >= 0) {
            sb.append(newickSafe(names[t.leaf]));
            return;
        }
        sb.append('(');
        appendNewick(t.left, names, sb);
        sb.append(':').append(f6(t.height - t.left.height)).append(',');
        appendNewick(t.right, names, sb);
        sb.append(':').append(f6(t.height - t.right.height)).append(')');
    }

    private static String newickSafe(String s) {
        StringBuilder b = new StringBuilder(s.length());
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            b.append((Character.isWhitespace(c) || c == '(' || c == ')' || c == ':' || c == ';'
                    || c == ',' || c == '\'' || c == '[' || c == ']') ? '_' : c);
        }
        return b.toString();
    }

    private void computePatterns() {
        Map<String, PatternGroup> map = new HashMap<>();
        for (Family f : families) {
            PatternGroup p = map.computeIfAbsent(f.pattern, k -> new PatternGroup(k, f.frequency));
            p.families++;
            p.bp += f.totalBp;
        }
        patterns = new ArrayList<>(map.values());
        patterns.sort(Comparator.comparingInt((PatternGroup p) -> -p.families)
                .thenComparingInt(p -> -p.frequency)
                .thenComparing(p -> p.pattern, Comparator.reverseOrder()));
    }

    // ── Accessors (used by the figure and the tests) ─────────────────────────
    public int genomeCount() {
        return nGenomes;
    }

    public String[] genomeNames() {
        return names.clone();
    }

    public long[] genomeLengths() {
        return genomeLength.clone();
    }

    public List<Family> families() {
        return Collections.unmodifiableList(families);
    }

    public int coreCount() {
        return coreCount;
    }

    public int accessoryCount() {
        return accessoryCount;
    }

    public int uniqueCount() {
        return uniqueCount;
    }

    /** Families present in exactly k genomes, index k = 1..N (index 0 unused). */
    public int[] spectrum() {
        return specCount.clone();
    }

    /** Expected pangenome size after sampling n genomes, index n = 1..N. */
    public double[] expectedPan() {
        return expPan.clone();
    }

    /** Expected core size after sampling n genomes, index n = 1..N. */
    public double[] expectedCore() {
        return expCore.clone();
    }

    /** Expected number of new families contributed by the n-th genome, index n = 1..N. */
    public double[] expectedNew() {
        return expNew.clone();
    }

    /** Heaps' law exponent (NaN when fewer than two usable points, i.e. N &lt; 3). */
    public double heapsAlpha() {
        return heapsAlpha;
    }

    public double heapsGamma() {
        return heapsGamma;
    }

    public double[][] jaccard() {
        return copy(jaccard);
    }

    public double[][] brayCurtis() {
        return copy(brayCurtis);
    }

    public Node tree() {
        return tree;
    }

    /** Genome indices in the left-to-right order of the UPGMA tree. */
    public int[] leafOrder() {
        return leafOrder.clone();
    }

    public String newick() {
        return newick(tree, names);
    }

    public List<PatternGroup> patterns() {
        return Collections.unmodifiableList(patterns);
    }

    /** Per genome: {families present, core, accessory, unique} counts. */
    public int[][] genomeCategoryCounts() {
        int[][] r = new int[nGenomes][4];
        for (Family f : families) {
            int c = switch (f.category) {
                case "core" -> 1;
                case "accessory" -> 2;
                default -> 3;
            };
            for (int g = 0; g < nGenomes; g++) {
                if (f.count[g] > 0) {
                    r[g][0]++;
                    r[g][c]++;
                }
            }
        }
        return r;
    }

    private static double[][] copy(double[][] m) {
        double[][] r = new double[m.length][];
        for (int i = 0; i < m.length; i++) {
            r[i] = m[i].clone();
        }
        return r;
    }

    // ── Formatting ──────────────────────────────────────────────────────────
    private static String pct(long part, long whole) {
        return (whole > 0) ? f2(100.0 * part / whole) + "%" : "0.00%";
    }

    static String f2(double v) {
        return String.format(Locale.ROOT, "%.2f", v);
    }

    static String f4(double v) {
        return String.format(Locale.ROOT, "%.4f", v);
    }

    static String f6(double v) {
        return String.format(Locale.ROOT, "%.6f", v);
    }

    private static String heapsVerdict(double alpha) {
        if (Double.isNaN(alpha)) {
            return "n/a";
        }
        return (alpha > 1.0) ? "closed" : "open";
    }

    /** Lists the genomes of a presence pattern, abbreviated after {@code max} names. */
    private String genomesOf(String pattern, int max) {
        StringBuilder sb = new StringBuilder();
        int listed = 0, total = 0;
        for (int g = 0; g < pattern.length(); g++) {
            if (pattern.charAt(g) == '1') {
                total++;
                if (listed < max) {
                    if (listed > 0) {
                        sb.append(',');
                    }
                    sb.append(names[g]);
                    listed++;
                }
            }
        }
        if (total > listed) {
            sb.append(",...(+").append(total - listed).append(')');
        }
        return sb.toString();
    }

    // ── Output ──────────────────────────────────────────────────────────────
    /**
     * Writes the whole report set next to {@code base}:
     * {@code _pangenome.txt} (summary), {@code _pangenome.tsv} (copy-number
     * matrix), {@code _pangenome_bp.tsv} (bp matrix),
     * {@code _pangenome_distance.tsv} (Jaccard distance matrix),
     * {@code _pangenome.nwk} (UPGMA tree) and {@code _pangenome.svg} (figure).
     */
    public void write(String base) throws IOException {
        writeSummary(base + "_pangenome.txt", base);
        writeMatrix(base + "_pangenome.tsv");
        writeBpMatrix(base + "_pangenome_bp.tsv");
        writeDistanceMatrix(base + "_pangenome_distance.tsv");
        writeNewick(base + "_pangenome.nwk");
        new PangenomeFigure(this).write(base + "_pangenome.svg");
    }

    private void writeSummary(String file, String base) throws IOException {
        int total = families.size();
        int N = nGenomes;
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome report: " + file);

            w.write("#TotalRepeats " + TotalRepeats.VERSION + " pangenome analysis\n");
            w.write("#Genomes (N=" + N + "), one genome = one " + (byFile
                    ? "input file (all records of a file count as one genome; -pangroup=record treats every record separately)"
                    : "sequence record (-pangroup=file groups the records of a file into one genome)") + ":\n");
            for (int g = 0; g < N; g++) {
                StringBuilder recs = new StringBuilder();
                int n = 0;
                for (int r = 0; r < recordGenome.length; r++) {
                    if (recordGenome[r] == g) {
                        if (n < 8) {
                            if (n > 0) {
                                recs.append(", ");
                            }
                            recs.append(recordNames[r]);
                        } else if (n == 8) {
                            recs.append(", ...");
                        }
                        n++;
                    }
                }
                w.write("#  " + (g + 1) + "\t" + names[g] + "\tlength=" + genomeLength[g]
                        + "\trecords=" + n + (byFile ? ": " + recs : "") + "\n");
            }
            w.write("#\n");
            w.write("#Definitions:\n");
            w.write("#  Family (CRP) = cluster of homologous repeat instances across the analysed genomes; ClusterID matches the .gff report.\n");
            w.write("#  Present in a genome = the family has >= 1 instance located in that genome.\n");
            String accRange = (N > 2) ? "present in 2.." + (N - 1) : "not possible with " + N + " genome" + (N == 1 ? "" : "s");
            w.write("#  core = present in all " + N + " genomes; accessory = " + accRange
                    + "; unique = present in exactly 1.\n");
            w.write("#  Tier = the four-level classification used by bacterial pangenome tools: core >= 99% of genomes, soft-core 95-99%, shell 15-95%, cloud < 15%.\n");
            w.write("#  Rarefaction = expected pangenome / core size when n genomes are sampled, averaged exactly over all orderings of the genomes.\n");
            w.write("#  Heaps' law (Tettelin et al. 2008): new families per added genome n(N) = kappa * N^-alpha; alpha <= 1 marks an open pangenome, alpha > 1 a closed one.\n");
            w.write("#  Jaccard = shared / union of families (1 = identical sets); jaccard_accessory excludes the core families; bray_curtis = sum|ci-cj| / sum(ci+cj) over copy numbers (0 = identical).\n");
            w.write("#  A repeat is detected only where it has >= 2 copies in the search space (per record in -collate/-combinemask, the whole set in -joint),\n");
            w.write("#  so single-copy sequence shared between genomes does not form a family and is not counted here.\n");
            w.write("#  Decimal separator is '.' regardless of the system locale.\n#\n");

            // ── Summary ──
            w.write("== Summary ==\n");
            w.write("Genomes: " + N + "\n");
            w.write("Repeat families (CRP clusters): " + total + "\n");
            w.write("  core      (in " + N + " genomes): " + coreCount
                    + "\t(" + pct(coreCount, total) + ")\tbp=" + coreBp + "\n");
            w.write("  accessory (" + ((N > 2) ? "in 2.." + (N - 1) : "n/a") + "):       " + accessoryCount
                    + "\t(" + pct(accessoryCount, total) + ")\tbp=" + accessoryBp + "\n");
            w.write("  unique    (in 1 genome):       " + uniqueCount
                    + "\t(" + pct(uniqueCount, total) + ")\tbp=" + uniqueBp + "\n");
            w.write("Tiers:");
            String[] tierRange = {">=99%", "95-99%", "15-95%", "<15%"};
            for (int t = 0; t < 4; t++) {
                w.write("  " + TIER_NAMES[t] + " (" + tierRange[t] + ")=" + tierCount[t]
                        + " (" + pct(tierCount[t], total) + ", bp=" + tierBp[t] + ")");
            }
            w.write("\n");
            double avgGenome = expPan[1];
            w.write("Pangenome size (families): " + total + "\tAverage genome (families): " + f2(avgGenome)
                    + "\tCore share of an average genome: " + (avgGenome > 0 ? f2(100.0 * coreCount / avgGenome) + "%" : "n/a") + "\n");
            if (Double.isNaN(heapsAlpha)) {
                w.write("Openness (Heaps' law): n/a (needs >= 3 genomes and new families after the second genome)\n");
            } else {
                w.write("Openness (Heaps' law): alpha=" + f4(heapsAlpha) + "\tkappa=" + f2(heapsKappa)
                        + "\t-> " + heapsVerdict(heapsAlpha) + " pangenome"
                        + (N < 5 ? " (estimate from only " + N + " genomes)" : "")
                        + "\tpan(n) ~ n^gamma, gamma=" + (Double.isNaN(heapsGamma) ? "n/a" : f4(heapsGamma)) + "\n");
            }
            w.write("\n");

            // ── Spectrum ──
            w.write("== Family frequency spectrum ==\n");
            w.write("Present_in_genomes\tFamilies\tPct\tBp\tTier\n");
            for (int k = 1; k <= N; k++) {
                w.write(k + "\t" + specCount[k] + "\t" + pct(specCount[k], total) + "\t" + specBp[k]
                        + "\t" + TIER_NAMES[tierOf(k, N)] + "\n");
            }
            w.write("\n");

            // ── Rarefaction ──
            w.write("== Rarefaction (expected over all genome orderings) ==\n");
            w.write("Genomes_sampled\tPan_families\tCore_families\tNew_families\n");
            for (int n = 1; n <= N; n++) {
                w.write(n + "\t" + f2(expPan[n]) + "\t" + f2(expCore[n]) + "\t" + f2(expNew[n]) + "\n");
            }
            w.write("\n");

            // ── Per-genome summary ──
            w.write("== Per-genome summary ==\n");
            w.write("Genome\tLength_bp\tFamilies_present\tCore\tAccessory\tUnique_to_this"
                    + "\tCore_bp\tAccessory_bp\tUnique_bp\tCRP_bp\tCRP_pct\tUCRP_n\tUCRP_bp\tUCRP_pct\tSTR_n\tSTR_bp\tSTR_pct\n");
            for (int g = 0; g < N; g++) {
                int present = 0, core = 0, acc = 0, uniq = 0;
                long cBp = 0, aBp = 0, uBp = 0;
                for (Family f : families) {
                    if (f.count[g] > 0) {
                        present++;
                        switch (f.category) {
                            case "core" -> {
                                core++;
                                cBp += f.bp[g];
                            }
                            case "accessory" -> {
                                acc++;
                                aBp += f.bp[g];
                            }
                            default -> {
                                uniq++;
                                uBp += f.bp[g];
                            }
                        }
                    }
                }
                long crpBp = cBp + aBp + uBp;
                long L = genomeLength[g];
                w.write(names[g] + "\t" + L + "\t" + present + "\t" + core + "\t" + acc + "\t" + uniq
                        + "\t" + cBp + "\t" + aBp + "\t" + uBp
                        + "\t" + crpBp + "\t" + pct(crpBp, L)
                        + "\t" + ucrpCount[g] + "\t" + ucrpBp[g] + "\t" + pct(ucrpBp[g], L)
                        + "\t" + strCount[g] + "\t" + strBp[g] + "\t" + pct(strBp[g], L) + "\n");
            }
            w.write("\n");

            // ── Pairwise ──
            w.write("== Pairwise genome comparison ==\n");
            w.write("genome_i\tgenome_j\tshared\tunion\tjaccard\tjaccard_accessory\tbray_curtis\n");
            for (int a = 0; a < N; a++) {
                for (int b = a + 1; b < N; b++) {
                    w.write(names[a] + "\t" + names[b] + "\t" + shared[a][b] + "\t" + union[a][b]
                            + "\t" + f4(jaccard[a][b]) + "\t" + f4(jaccardAccessory[a][b])
                            + "\t" + f4(brayCurtis[a][b]) + "\n");
                }
            }
            w.write("\n");

            // ── Tree ──
            w.write("== Genome clustering (UPGMA on Jaccard distance = 1 - jaccard) ==\n");
            w.write(newick() + "\n");
            w.write("(also written to " + new File(base + "_pangenome.nwk").getName()
                    + "; the full distance matrix is in " + new File(base + "_pangenome_distance.tsv").getName() + ")\n\n");

            // ── Patterns ──
            w.write("== Presence patterns (top " + Math.min(TOP_PATTERNS, patterns.size()) + " of "
                    + patterns.size() + "; 1 = present, genome order as listed in the header) ==\n");
            w.write("Pattern\tGenomes\tFamilies\tBp\tGenome_names\n");
            int shown = 0;
            for (PatternGroup p : patterns) {
                if (shown++ >= TOP_PATTERNS) {
                    break;
                }
                w.write(p.pattern + "\t" + p.frequency + "\t" + p.families + "\t" + p.bp
                        + "\t" + genomesOf(p.pattern, 6) + "\n");
            }
            w.write("\n");

            // ── Copy-number variation ──
            w.write("== Copy-number variation ==\n");
            w.write("[core families with the most variable copy number (top " + TOP_CNV + " by CV)]\n");
            w.write("ClusterID\tRefLabel\tMinCopies\tMaxCopies\tMeanCopies\tCV\tFold\n");
            List<Family> cores = new ArrayList<>();
            for (Family f : families) {
                if ("core".equals(f.category)) {
                    cores.add(f);
                }
            }
            cores.sort(Comparator.comparingDouble((Family f) -> -f.cv).thenComparingInt(f -> f.clusterId));
            shown = 0;
            for (Family f : cores) {
                if (shown++ >= TOP_CNV) {
                    break;
                }
                w.write(f.clusterId + "\t" + f.refLabel + "\t" + f.minCopies + "\t" + f.maxCopies
                        + "\t" + f2(f.meanCopies) + "\t" + f4(f.cv)
                        + "\t" + f2((double) f.maxCopies / Math.max(1, f.minCopies)) + "\n");
            }
            w.write("\n[per-genome expansions: >= 2x the mean copy number of the other genomes and >= 3 extra copies; families present in >= 2 genomes]\n");
            w.write("Genome\tClusterID\tRefLabel\tCopies_here\tMean_others\tFold\tExtra_copies\n");
            for (int g = 0; g < N; g++) {
                List<double[]> exp = new ArrayList<>();   // {familyIndex, copies, meanOthers}
                for (int i = 0; i < families.size(); i++) {
                    Family f = families.get(i);
                    if (f.frequency < 2 || N < 2) {
                        continue;
                    }
                    double others = (double) (f.totalMembers - f.count[g]) / (N - 1);
                    if (f.count[g] >= 2.0 * others && f.count[g] - others >= 3.0) {
                        exp.add(new double[]{i, f.count[g], others});
                    }
                }
                exp.sort((p, q) -> Double.compare((q[1] - q[2]), (p[1] - p[2])));
                shown = 0;
                for (double[] e : exp) {
                    if (shown++ >= TOP_CNV) {
                        break;
                    }
                    Family f = families.get((int) e[0]);
                    String fold = (e[2] > 0) ? f2(e[1] / e[2]) : "inf";
                    w.write(names[g] + "\t" + f.clusterId + "\t" + f.refLabel + "\t" + (int) e[1]
                            + "\t" + f2(e[2]) + "\t" + fold + "\t" + f2(e[1] - e[2]) + "\n");
                }
            }
            w.write("\n");

            // ── Core families ──
            w.write("== Core families (common to all genomes; by bp) ==\n");
            w.write("ClusterID\tRefLabel\tMembers\tBp\tMinCopies\tMaxCopies\tCV\n");
            cores.sort(Comparator.comparingLong((Family f) -> -f.totalBp).thenComparingInt(f -> f.clusterId));
            for (Family f : cores) {
                w.write(f.clusterId + "\t" + f.refLabel + "\t" + f.totalMembers + "\t" + f.totalBp
                        + "\t" + f.minCopies + "\t" + f.maxCopies + "\t" + f4(f.cv) + "\n");
            }
            w.write("\n");

            // ── Unique families ──
            w.write("== Genome-specific (unique) families (by bp) ==\n");
            for (int g = 0; g < N; g++) {
                final int gg = g;
                w.write("[" + names[g] + "]\n");
                w.write("ClusterID\tRefLabel\tMembers\tBp\n");
                List<Family> u = new ArrayList<>();
                for (Family f : families) {
                    if ("unique".equals(f.category) && f.count[g] > 0) {
                        u.add(f);
                    }
                }
                u.sort(Comparator.comparingLong((Family f) -> -f.bp[gg]).thenComparingInt(f -> f.clusterId));
                for (Family f : u) {
                    w.write(f.clusterId + "\t" + f.refLabel + "\t" + f.count[g] + "\t" + f.bp[g] + "\n");
                }
                w.write("\n");
            }
        }
    }

    /**
     * Machine-readable copy-number matrix: one row per CRP family (plus the STR
     * and UCRP aggregates) and one column per genome. The first columns keep
     * the layout of earlier releases; the new ones are appended at the end.
     */
    private void writeMatrix(String file) throws IOException {
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome matrix: " + file);
            StringBuilder header = new StringBuilder("ClusterID\tType\tCategory\tFrequency\tRefLabel");
            for (int g = 0; g < nGenomes; g++) {
                header.append("\t").append(names[g]);
            }
            header.append("\tTotalMembers\tTotalBp\tTier\tPattern\tMinCopies\tMaxCopies\tMeanCopies\tCV\n");
            w.write(header.toString());

            writeAggregateRow(w, 1, "STR", strCount, strBp);
            writeAggregateRow(w, 2, "UCRP", ucrpCount, ucrpBp);
            for (Family f : families) {
                StringBuilder row = new StringBuilder();
                row.append(f.clusterId).append("\tCRP\t").append(f.category)
                        .append("\t").append(f.frequency).append("\t").append(f.refLabel);
                for (int g = 0; g < nGenomes; g++) {
                    row.append("\t").append(f.count[g]);
                }
                row.append("\t").append(f.totalMembers).append("\t").append(f.totalBp)
                        .append("\t").append(f.tier).append("\t").append(f.pattern)
                        .append("\t").append(f.minCopies).append("\t").append(f.maxCopies)
                        .append("\t").append(f2(f.meanCopies)).append("\t").append(f4(f.cv)).append("\n");
                w.write(row.toString());
            }
        }
    }

    private void writeAggregateRow(BufferedWriter w, int id, String type, int[] cnt, long[] bp) throws IOException {
        StringBuilder row = new StringBuilder();
        StringBuilder pat = new StringBuilder();
        int freq = 0, members = 0, min = Integer.MAX_VALUE, max = 0;
        long total = 0;
        for (int g = 0; g < nGenomes; g++) {
            if (cnt[g] > 0) {
                freq++;
            }
            pat.append(cnt[g] > 0 ? '1' : '0');
            members += cnt[g];
            total += bp[g];
            min = Math.min(min, cnt[g]);
            max = Math.max(max, cnt[g]);
        }
        row.append(id).append("\t").append(type).append("\taggregate\t").append(freq).append("\t-");
        for (int g = 0; g < nGenomes; g++) {
            row.append("\t").append(cnt[g]);
        }
        row.append("\t").append(members).append("\t").append(total).append("\t-\t").append(pat)
                .append("\t").append(min).append("\t").append(max)
                .append("\t").append(f2((double) members / nGenomes)).append("\t-\n");
        w.write(row.toString());
    }

    /** Covered bp per family and genome, same row layout as the copy-number matrix. */
    private void writeBpMatrix(String file) throws IOException {
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome bp matrix: " + file);
            StringBuilder header = new StringBuilder("ClusterID\tType\tCategory\tFrequency\tRefLabel");
            for (int g = 0; g < nGenomes; g++) {
                header.append("\t").append(names[g]);
            }
            header.append("\tTotalBp\n");
            w.write(header.toString());
            writeAggregateBpRow(w, 1, "STR", strCount, strBp);
            writeAggregateBpRow(w, 2, "UCRP", ucrpCount, ucrpBp);
            for (Family f : families) {
                StringBuilder row = new StringBuilder();
                row.append(f.clusterId).append("\tCRP\t").append(f.category)
                        .append("\t").append(f.frequency).append("\t").append(f.refLabel);
                for (int g = 0; g < nGenomes; g++) {
                    row.append("\t").append(f.bp[g]);
                }
                row.append("\t").append(f.totalBp).append("\n");
                w.write(row.toString());
            }
        }
    }

    private void writeAggregateBpRow(BufferedWriter w, int id, String type, int[] cnt, long[] bp) throws IOException {
        StringBuilder row = new StringBuilder();
        int freq = 0;
        long total = 0;
        for (int g = 0; g < nGenomes; g++) {
            if (cnt[g] > 0) {
                freq++;
            }
            total += bp[g];
        }
        row.append(id).append("\t").append(type).append("\taggregate\t").append(freq).append("\t-");
        for (int g = 0; g < nGenomes; g++) {
            row.append("\t").append(bp[g]);
        }
        row.append("\t").append(total).append("\n");
        w.write(row.toString());
    }

    /** Square Jaccard distance matrix (1 - jaccard over all families), genome names on both axes. */
    private void writeDistanceMatrix(String file) throws IOException {
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome distance matrix: " + file);
            StringBuilder sb = new StringBuilder("JaccardDistance");
            for (int g = 0; g < nGenomes; g++) {
                sb.append("\t").append(names[g]);
            }
            sb.append("\n");
            for (int a = 0; a < nGenomes; a++) {
                sb.append(names[a]);
                for (int b = 0; b < nGenomes; b++) {
                    sb.append("\t").append(f4(1.0 - jaccard[a][b]));
                }
                sb.append("\n");
            }
            w.write(sb.toString());
        }
    }

    private void writeNewick(String file) throws IOException {
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome tree: " + file);
            w.write(newick());
            w.write("\n");
        }
    }
}
