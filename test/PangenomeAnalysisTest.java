import java.util.ArrayList;
import java.util.List;

/**
 * Pangenome statistics on a hand-built cluster table, where every expected
 * value can be worked out by hand: classification, spectrum, exact
 * rarefaction (checked against brute-force enumeration of genome subsets),
 * Heaps' law on a synthetic open pangenome, Jaccard / Bray–Curtis distances,
 * the UPGMA tree, presence patterns, copy-number statistics and the grouping of
 * records into genomes.
 */
final class PangenomeAnalysisTest {

    private PangenomeAnalysisTest() {
    }

    static void run() {
        TestRunner.suite("PangenomeAnalysis");
        classification();
        rarefactionExact();
        heapsOpenPangenome();
        heapsClosedPangenome();
        distancesAndTree();
        upgmaKnownTree();
        patternsAndCopyNumber();
        groupingByFile();
        duplicateNames();
        tiers();
    }

    // ── Fixture ─────────────────────────────────────────────────────────────
    // Three records of 1000 bp each: a = [0,1000), b = [1000,2000), c = [2000,3000).
    private static final long[] SEQSLEN = {1000, 2000, 3000};
    private static final String[] NAMES = {"a", "b", "c"};

    /** Blocks as (start, length) pairs; a negative length marks the reverse strand. */
    private static long[] blocks(long... v) {
        return v;
    }

    /**
     * bb[0] = STR, bb[1] = UCRP, bb[2..] = CRP families:
     *   id3: a×2, b×1, c×1  (core)
     *   id4: a×2            (unique to a)
     *   id5: a×1, b×3       (accessory)
     *   id6: null           (placeholder, skipped)
     *   id7: b×1, c×1       (accessory)
     */
    private static ArrayList<long[]> fixture() {
        ArrayList<long[]> bb = new ArrayList<>();
        bb.add(blocks(10, 20, 1010, 30, 2010, 40));                     // STR: one per record
        bb.add(blocks(100, 50, 2100, 60));                              // UCRP: a, c
        bb.add(blocks(200, 100, 300, -100, 1200, 100, 2200, 100));      // family 3: a,a,b,c
        bb.add(blocks(500, 80, 600, 80));                               // family 4: a,a
        bb.add(blocks(700, 50, 1500, 50, 1600, -50, 1700, 50));         // family 5: a,b,b,b
        bb.add(null);                                                   // family 6: empty
        bb.add(blocks(1800, 70, 2800, 70));                             // family 7: b,c
        return bb;
    }

    private static PangenomeAnalysis analysis() {
        return new PangenomeAnalysis(fixture(), SEQSLEN, NAMES, null, null);
    }

    private static PangenomeAnalysis.Family family(PangenomeAnalysis pa, int clusterId) {
        for (PangenomeAnalysis.Family f : pa.families()) {
            if (f.clusterId == clusterId) {
                return f;
            }
        }
        return null;
    }

    private static boolean near(double a, double b) {
        return Math.abs(a - b) < 1e-9;
    }

    // ── Checks ──────────────────────────────────────────────────────────────
    private static void classification() {
        PangenomeAnalysis pa = analysis();
        TestRunner.eq("three records are three genomes", 3, pa.genomeCount());
        TestRunner.eq("empty cluster rows are skipped", 4, pa.families().size());
        TestRunner.eq("core count", 1, pa.coreCount());
        TestRunner.eq("accessory count", 2, pa.accessoryCount());
        TestRunner.eq("unique count", 1, pa.uniqueCount());

        PangenomeAnalysis.Family f3 = family(pa, 3);
        TestRunner.eq("family 3 is core", "core", f3.category);
        TestRunner.eq("family 3 copies in a", 2, f3.count[0]);
        TestRunner.eq("family 3 bp counts reverse-strand copies by |length|", 400, f3.totalBp);
        TestRunner.eq("family 3 pattern", "111", f3.pattern);
        TestRunner.eq("family 4 is unique", "unique", family(pa, 4).category);
        TestRunner.eq("family 4 pattern", "100", family(pa, 4).pattern);
        TestRunner.eq("family 5 is accessory", "accessory", family(pa, 5).category);
        TestRunner.eq("family 7 pattern", "011", family(pa, 7).pattern);

        TestRunner.eq("spectrum", new int[]{0, 1, 2, 1}, pa.spectrum());
        int[][] cc = pa.genomeCategoryCounts();
        TestRunner.eq("genome a: families present", 3, cc[0][0]);
        TestRunner.eq("genome a: unique", 1, cc[0][3]);
        TestRunner.eq("genome c: families present", 2, cc[2][0]);
        TestRunner.eq("genome lengths follow the record boundaries", new long[]{1000, 1000, 1000}, pa.genomeLengths());
    }

    /** The closed-form expectation must equal the average over every subset of genomes. */
    private static void rarefactionExact() {
        PangenomeAnalysis pa = analysis();
        List<PangenomeAnalysis.Family> fams = pa.families();
        int N = pa.genomeCount();
        double[] pan = pa.expectedPan();
        double[] core = pa.expectedCore();
        for (int n = 1; n <= N; n++) {
            double sumPan = 0, sumCore = 0;
            int subsets = 0;
            for (int mask = 1; mask < (1 << N); mask++) {
                if (Integer.bitCount(mask) != n) {
                    continue;
                }
                subsets++;
                for (PangenomeAnalysis.Family f : fams) {
                    boolean any = false, all = true;
                    for (int g = 0; g < N; g++) {
                        if ((mask & (1 << g)) == 0) {
                            continue;
                        }
                        if (f.count[g] > 0) {
                            any = true;
                        } else {
                            all = false;
                        }
                    }
                    if (any) {
                        sumPan++;
                    }
                    if (all) {
                        sumCore++;
                    }
                }
            }
            TestRunner.ok("expected pan size at n=" + n + " matches brute force",
                    near(pan[n], sumPan / subsets));
            TestRunner.ok("expected core size at n=" + n + " matches brute force",
                    near(core[n], sumCore / subsets));
        }
        TestRunner.ok("pan(N) is the whole pangenome", near(pan[N], 4));
        TestRunner.ok("core(N) is the core", near(core[N], 1));
        TestRunner.ok("pan(1) is the average genome", near(pan[1], (3 + 3 + 2) / 3.0));
        double[] nw = pa.expectedNew();
        TestRunner.ok("new(1) equals pan(1)", near(nw[1], pan[1]));
        TestRunner.ok("new(n) telescopes to pan(N)", near(nw[1] + nw[2] + nw[3], pan[N]));
    }

    /** Every genome brings only its own families: new(n) is constant, so alpha = 0 (open). */
    private static void heapsOpenPangenome() {
        int N = 5;
        long[] ends = new long[N];
        String[] names = new String[N];
        ArrayList<long[]> bb = new ArrayList<>();
        bb.add(new long[0]);
        bb.add(new long[0]);
        for (int g = 0; g < N; g++) {
            ends[g] = 1000L * (g + 1);
            names[g] = "g" + g;
            long base = 1000L * g;
            for (int k = 0; k < 3; k++) {                       // three private families per genome
                bb.add(blocks(base + 10 + 40 * k, 20, base + 30 + 40 * k, 20));
            }
        }
        PangenomeAnalysis pa = new PangenomeAnalysis(bb, ends, names, null, null);
        TestRunner.eq("all families unique", 15, pa.uniqueCount());
        double[] nw = pa.expectedNew();
        TestRunner.ok("new families per genome is constant", near(nw[2], 3) && near(nw[5], 3));
        TestRunner.ok("Heaps' alpha is 0 for a fully open pangenome", near(pa.heapsAlpha(), 0));
        TestRunner.ok("gamma is 1: pan grows linearly", near(pa.heapsGamma(), 1));
        TestRunner.ok("Jaccard between disjoint genomes is 0", near(pa.jaccard()[0][1], 0));
    }

    /** Every family is core: nothing new after the first genome, so alpha is undefined. */
    private static void heapsClosedPangenome() {
        int N = 4;
        long[] ends = new long[N];
        String[] names = new String[N];
        for (int g = 0; g < N; g++) {
            ends[g] = 1000L * (g + 1);
            names[g] = "g" + g;
        }
        ArrayList<long[]> bb = new ArrayList<>();
        bb.add(new long[0]);
        bb.add(new long[0]);
        long[] fam = new long[2 * N];
        for (int g = 0; g < N; g++) {
            fam[2 * g] = 1000L * g + 50;
            fam[2 * g + 1] = 100;
        }
        bb.add(fam);
        PangenomeAnalysis pa = new PangenomeAnalysis(bb, ends, names, null, null);
        TestRunner.eq("one core family", 1, pa.coreCount());
        TestRunner.ok("expected core stays 1 at every n", near(pa.expectedCore()[1], 1) && near(pa.expectedCore()[N], 1));
        TestRunner.ok("no new families after the first genome", near(pa.expectedNew()[2], 0));
        TestRunner.ok("alpha is undefined without new families", Double.isNaN(pa.heapsAlpha()));
        TestRunner.ok("identical genomes are at Jaccard 1", near(pa.jaccard()[0][3], 1));
        TestRunner.ok("identical genomes have Bray-Curtis 0", near(pa.brayCurtis()[0][3], 0));
        TestRunner.ok("the tree joins identical genomes at height 0",
                pa.newick().contains(":0.000000"));
    }

    private static void distancesAndTree() {
        PangenomeAnalysis pa = analysis();
        double[][] j = pa.jaccard();
        // a = {3,4,5}, b = {3,5,7}, c = {3,7}
        TestRunner.ok("Jaccard(a,b) = 2/4", near(j[0][1], 0.5));
        TestRunner.ok("Jaccard(a,c) = 1/4", near(j[0][2], 0.25));
        TestRunner.ok("Jaccard(b,c) = 2/3", near(j[1][2], 2.0 / 3));
        TestRunner.ok("Jaccard is symmetric", near(j[0][1], j[1][0]));
        TestRunner.ok("Jaccard of a genome with itself is 1", near(j[2][2], 1));
        double[][] bc = pa.brayCurtis();
        // copies a = (2,2,1,0), b = (1,0,3,1): |diff| = 1+2+2+1 = 6, sum = 10
        TestRunner.ok("Bray-Curtis(a,b) = 6/10", near(bc[0][1], 0.6));

        String nwk = pa.newick();
        TestRunner.ok("Newick ends with a semicolon", nwk.endsWith(";"));
        TestRunner.ok("Newick names every genome", nwk.contains("a:") && nwk.contains("b:") && nwk.contains("c:"));
        TestRunner.ok("b and c (closest pair) are sisters", nwk.contains("(b:") || nwk.contains(",b:"));
        TestRunner.ok("the outer group is a", nwk.startsWith("(a:") || nwk.endsWith("a:0.375000);") || nwk.contains(",a:"));
        int[] order = pa.leafOrder();
        TestRunner.eq("leaf order lists every genome once", 3, order.length);
        TestRunner.ok("b and c are adjacent in leaf order",
                Math.abs(indexOf(order, 1) - indexOf(order, 2)) == 1);
    }

    private static int indexOf(int[] a, int v) {
        for (int i = 0; i < a.length; i++) {
            if (a[i] == v) {
                return i;
            }
        }
        return -1;
    }

    /** UPGMA on a hand-solved matrix: d(a,b)=0.2, d(a,c)=d(b,c)=0.8 -> ((a:0.1,b:0.1):0.3,c:0.4). */
    private static void upgmaKnownTree() {
        double[][] d = {
            {0.0, 0.2, 0.8},
            {0.2, 0.0, 0.8},
            {0.8, 0.8, 0.0}
        };
        PangenomeAnalysis.Node t = PangenomeAnalysis.upgma(d);
        TestRunner.ok("root height is half the final merge distance", near(t.height, 0.4));
        TestRunner.eq("UPGMA Newick with branch lengths",
                "((a:0.100000,b:0.100000):0.300000,c:0.400000);",
                PangenomeAnalysis.newick(t, new String[]{"a", "b", "c"}));
        TestRunner.eq("Newick-unsafe characters in names are replaced",
                "((x_1:0.100000,y_2:0.100000):0.300000,z_3:0.400000);",
                PangenomeAnalysis.newick(t, new String[]{"x 1", "y(2", "z:3"}));

        // Unequal cluster sizes must be weighted in the distance update.
        double[][] e = {
            {0.0, 0.1, 0.6, 0.9},
            {0.1, 0.0, 0.4, 0.9},
            {0.6, 0.4, 0.0, 0.9},
            {0.9, 0.9, 0.9, 0.0}
        };
        PangenomeAnalysis.Node u = PangenomeAnalysis.upgma(e);
        // (a,b) at 0.1 -> d(ab,c) = (0.6+0.4)/2 = 0.5 -> (ab,c) at 0.5 -> d(abc,d) = 0.9
        TestRunner.ok("weighted update: root at 0.45", near(u.height, 0.45));
        TestRunner.ok("weighted update: (ab,c) joins at 0.25", near(u.left.height, 0.25) || near(u.right.height, 0.25));
    }

    private static void patternsAndCopyNumber() {
        PangenomeAnalysis pa = analysis();
        List<PangenomeAnalysis.PatternGroup> pats = pa.patterns();
        TestRunner.eq("four distinct patterns", 4, pats.size());
        for (PangenomeAnalysis.PatternGroup p : pats) {
            TestRunner.eq("pattern " + p.pattern + " holds one family", 1, p.families);
        }
        TestRunner.ok("patterns are ordered by family count then frequency",
                pats.get(0).frequency >= pats.get(1).frequency);

        PangenomeAnalysis.Family f5 = family(pa, 5);          // copies (1, 3, 0)
        TestRunner.eq("min copies over all genomes includes absence", 0, f5.minCopies);
        TestRunner.eq("max copies", 3, f5.maxCopies);
        TestRunner.ok("mean copies", near(f5.meanCopies, 4.0 / 3));
        double var = (1 + 9 + 0) / 3.0 - (4.0 / 3) * (4.0 / 3);
        TestRunner.ok("coefficient of variation", near(f5.cv, Math.sqrt(var) / (4.0 / 3)));
        PangenomeAnalysis.Family f3 = family(pa, 3);          // copies (2, 1, 1)
        TestRunner.eq("core family min copies", 1, f3.minCopies);
        TestRunner.eq("core family max copies", 2, f3.maxCopies);
        TestRunner.eq("tier of a 3/3 family", "core", f3.tier);
    }

    /** Records a and b belong to one file: families 4 (a) and 7 (b,c) change meaning. */
    private static void groupingByFile() {
        PangenomeAnalysis pa = new PangenomeAnalysis(fixture(), SEQSLEN, NAMES, null, null,
                new int[]{0, 0, 1}, new String[]{"one.fa", "two.fa"});
        TestRunner.eq("two files are two genomes", 2, pa.genomeCount());
        TestRunner.eq("genome names come from the files", "one.fa", pa.genomeNames()[0]);
        TestRunner.eq("genome length sums its records", new long[]{2000, 1000}, pa.genomeLengths());
        TestRunner.eq("family 3 is still core", "core", family(pa, 3).category);
        TestRunner.eq("family 4 is unique to the first file", "10", family(pa, 4).pattern);
        TestRunner.eq("family 5 (a,b) becomes unique to the first file", "unique", family(pa, 5).category);
        TestRunner.eq("family 7 (b,c) spans both files", "core", family(pa, 7).category);
        TestRunner.eq("family 7 copies in the first file", 1, family(pa, 7).count[0]);
        TestRunner.eq("spectrum over files", new int[]{0, 2, 2}, pa.spectrum());
    }

    private static void duplicateNames() {
        PangenomeAnalysis pa = new PangenomeAnalysis(fixture(), SEQSLEN, new String[]{"chr1", "chr1", "chr2"}, null, null);
        String[] n = pa.genomeNames();
        TestRunner.eq("duplicate record names are indexed", "chr1#1", n[0]);
        TestRunner.eq("second duplicate", "chr1#2", n[1]);
        TestRunner.eq("a unique name is left alone", "chr2", n[2]);
    }

    private static void tiers() {
        TestRunner.eq("100/100 is core", 0, PangenomeAnalysis.tierOf(100, 100));
        TestRunner.eq("99/100 is core", 0, PangenomeAnalysis.tierOf(99, 100));
        TestRunner.eq("96/100 is soft-core", 1, PangenomeAnalysis.tierOf(96, 100));
        TestRunner.eq("50/100 is shell", 2, PangenomeAnalysis.tierOf(50, 100));
        TestRunner.eq("15/100 is shell", 2, PangenomeAnalysis.tierOf(15, 100));
        TestRunner.eq("14/100 is cloud", 3, PangenomeAnalysis.tierOf(14, 100));
        TestRunner.eq("1/3 is shell (small N)", 2, PangenomeAnalysis.tierOf(1, 3));
        TestRunner.eq("3/3 is core", 0, PangenomeAnalysis.tierOf(3, 3));

        double[] x = {0, 1, 2, 4, 8};
        double[] y = {0, 2, 1, 0.5, 0.25};         // y = 2 * x^-1
        double[] fit = PangenomeAnalysis.powerFit(x, y, 1, 4);
        TestRunner.ok("power-law fit recovers the exponent", near(fit[0], -1));
        TestRunner.ok("power-law fit recovers the coefficient", near(fit[1], 2));
        TestRunner.ok("a single point cannot be fitted", PangenomeAnalysis.powerFit(x, y, 1, 1) == null);
    }
}
