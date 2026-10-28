import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * Clusters DNA sequences based on k-mer frequency vector similarity. Sequences
 * are represented as 4-mer frequency vectors and clustered by comparing
 * pairwise ratio profiles. Both forward and reverse-complement orientations are
 * considered. Optionally matches clusters against a set of reference sequences.
 *
 * Coordinates are global {@code long} offsets into a {@link SeqStore} (a
 * segmented view over many chromosomes), so the combined data may exceed the
 * ~2.1 Gb String limit. Each clustered item is a repeat block that lies fully
 * within one segment.
 *
 * <h2>How the profile measure works, and what it used to get wrong</h2>
 * A block's profile is its raw 4-mer count vector over a 238-entry palette, held as
 * sparse CSR rows. Two blocks are compared by the COSINE of those vectors, and a
 * candidate joins a seed when the cosine clears a chance floor that is calibrated,
 * per length scale, from this run's own data ({@link #calibrateProfileFloors}). When
 * the seed is much longer than the candidate, the comparison is made against the best
 * scale-matched WINDOW of the seed rather than against the seed as a whole
 * ({@link #buildWindowDense}).
 *
 * <p>All three of those are corrections to a measure that was, when checked against
 * ground truth independent of it (exact 25-mer sharing, adjudicated by Smith-Waterman
 * against a shuffled null), wrong in both directions at once — it missed most real
 * homology AND ~45% of what it grouped was noise. On MF782455 (1.39 Mb, 765 blocks):
 *
 * <pre>
 *                          blocks   families   real members   verified (SW)
 *   -vector, before          205       79            67            52%
 *   -vector, after           436       74           354            98%
 *   -contain (untouched)     465       92           365           100%
 * </pre>
 *
 * The same holds away from that genome. On E. coli MG1655 (NC_000913, 4.6 Mb) -vector goes
 * from 93 blocks at 60% real to 294 at 100% real (100 of 100 sampled pairs confirmed by
 * alignment) — MORE than -contain finds there (210). On a 7.5 Mb fungal chromosome, 384
 * blocks at 50% real to 600 at 88%. Across 42 Mb of E. coli with -collate, 2293 blocks to
 * 5594, in 85 s against the old 62 s.
 *
 * <p>In those runs only half to three-quarters of the members share an exact 25-mer with
 * their seed, yet alignment confirms them — composition sees copies that have diverged past
 * any exact k-mer match, which containment by construction cannot. That is the case for
 * keeping a composition measure alongside containment rather than replacing it.
 *
 * The three causes, each measured rather than reasoned:
 * <ul>
 * <li><b>No chance control.</b> The containment measure always had one
 * ({@link #containmentFloor}); the profile compared every pair at every length against
 * one flat threshold. The null score of two UNRELATED 250 bp windows has p99 = 65, so
 * ~1% of unrelated short pairs cleared a bar of 70 per comparison — and the greedy loop
 * offers every candidate to ~90 seeds in turn. Meanwhile at 2 kb the null never once
 * exceeded -1 across 20000 draws, so the same bar was discarding real homology by a huge
 * margin. No flat threshold can be right for both.</li>
 * <li><b>The wrong statistic.</b> The old ratio-of-ratios similarity read only the
 * dimensions two vectors shared, discarding which 4-mers were ABSENT — the dominant
 * signal for a short block. It also needed a count &gt; 3 to see a slot at all, which
 * left 40% of blocks with too few non-zeros to be compared at all. See
 * {@link #cosineMerge}.</li>
 * <li><b>Whole-block averaging.</b> Real repeat homology is partial; averaging a shared
 * region into a whole block dilutes it away before any comparison happens. See
 * {@link #buildWindowDense}.</li>
 * </ul>
 *
 * <h2>Representation notes</h2>
 * <ul>
 * <li><b>Sparse CSR vectors.</b> Instead of {@code int[numSeqs][256]} dense rows we keep
 * three flat arrays (idx, val, ptr) per orientation. Cosine merges two sparse vectors in
 * O(nnz) instead of scanning all 256 slots, and the flat layout removes the per-row
 * allocations / pointer chasing that pressured GC and memory bandwidth (the real limiter
 * for multi-core scaling).</li>
 * <li><b>Fused forward+reverse scan.</b> One pass per seed instead of two: each candidate
 * is tested forward first (priority preserved), and only reverse if forward failed.</li>
 * <li><b>Active set.</b> The greedy outer loop stays sequential and deterministic, but
 * candidates are drawn from a compacted list of still-unclustered indices, so each scan
 * touches only live candidates.</li>
 * <li><b>Reference matching.</b> The forward vectors built during clustering are reused
 * (no per-member recomputation from the store), and members are grouped by cluster in a
 * single pass instead of re-scanning all sequences per cluster.</li>
 * </ul>
 */
public final class SequencesClustering {

    // ── Tuning constants ────────────────────────────────────────────────────────
    private static final int KMER_SIZE = 4;
    private static final int MIN_REF_LENGTH = 80;     // min length of a reference sequence

    // ── Containment backend tuning (the -contain mode) ──────────────────────────
    // A separate, asymmetric measure for sequences of very different length: a
    // short block is absorbed into a longer seed when most of its canonical
    // k-mers occur in that seed. Ported from GeneDistance's -contain. The 4-mer
    // ratio-profile path above is left untouched.
    private static final int CONTAIN_K_MIN = 11;         // odd: a k-mer cannot be its own reverse complement
    private static final int CONTAIN_K_MAX = 24;         // 2*k <= 48 bits, always safe in a long
    private static final double MIN_CONTAINMENT = 0.50;  // share of a block's k-mers that must occur in the seed
    private static final int MIN_SHARED_KMERS = 5;       // fewer shared k-mers than this is treated as chance
    private static final long[] EMPTY_L = new long[0];

    /**
     * How many standard deviations clear of the chance level a containment ratio
     * must stand to count as evidence. One-sided 5 sigma is ~3e-7 per test, which
     * holds the expected false positives well under one across the ~1e5-1e6
     * candidate tests a large run makes. See {@link #containmentFloor}.
     */
    private static final double CHANCE_SIGMA = 5.0;

    // ── Profile chance floor and local windowing (the -vector mode) ─────────────
    // Three coupled corrections to the profile measure — a per-scale chance floor
    // (calibrateProfileFloors), scale-matched seed windows (buildWindowDense) and cosine
    // over raw counts (cosineMerge). They exist because the profile was measured, against
    // ground truth independent of it, to be wrong in BOTH directions at once: it missed
    // most real homology AND ~45% of what it grouped was noise. Each method's javadoc
    // carries the measurement that motivates it.

    /**
     * Number of windows drawn per length class when calibrating the chance floor. They
     * are compared all-against-all, so this yields n*(n-1)/2 null scores — 384 windows
     * give ~73k — while paying only 384 vector builds. That ratio is the whole reason the
     * calibration is affordable: building vectors is the cost, and scoring pairs is nearly
     * free by comparison.
     *
     * <p>This is not just a precision knob, it is THE precision knob: the floor is the
     * maximum of the null sample, so the maximum of n*(n-1)/2 draws estimates the
     * ~1-2/n^2 quantile, and raising n raises the bar. It should be read against how many
     * comparisons a run makes — every candidate is offered to every seed, so a few hundred
     * blocks already means ~1e4-1e5 comparisons, and the floor has to be a tail quantile
     * of that order or chance matches get through. Swept end-to-end on two genomes
     * (verified = member shares exact 25-mers with its seed):
     *
     * <pre>
     *          ch3_MG8 (7.5 Mb)              MF782455 (1.4 Mb)
     *   n      blocks  verified  real        blocks  verified  real     time (ch3_MG8)
     *   192      596     73.4%    369          525     92.7%    416         12 s
     *   384      602     83.9%    421          454     96.1%    370         17 s
     *   640      604      —        —           445      —        —          30 s
     * </pre>
     *
     * 384 is the chosen point: on ch3_MG8 it is strictly better than 192 on BOTH axes —
     * more precise and MORE real homology found, which looks contradictory until you
     * remember the greedy loop is not monotone in the threshold. A lax floor lets an early
     * promiscuous seed swallow the hub blocks that would otherwise have nucleated their
     * own families, so tightening it can raise recall. On MF782455 it trades 46 real
     * members for 3.4 points of precision. 640 costs nearly twice the time for no further
     * gain.
     */
    private static final int NULL_WINDOWS = 384;

    /**
     * Cap on draws per usable null window. Short classes reject most draws (at 128 bp
     * only ~3% of windows yield a usable vector), so the draw loop needs headroom; the
     * cap stops it spinning when a scale is simply unusable. Draws are cheap — a
     * substring, a shuffle and one vector build.
     */
    private static final int NULL_DRAW_LIMIT = 80;

    /**
     * Fewest usable null windows that still make a floor meaningful. Below this the
     * class is declared unusable rather than given a floor estimated from a handful of
     * pairs — an underestimated floor is exactly how noise gets in.
     */
    private static final int MIN_NULL_WINDOWS = 32;

    /**
     * Fewest non-homologous null PAIRS that still make a floor meaningful. Distinct from
     * {@link #MIN_NULL_WINDOWS}: enough windows can be drawn and still leave almost no
     * usable pairs, because a repetitive genome has most window pairs sharing a k-mer and
     * therefore screened out of the null. The floor is a maximum over the surviving pairs,
     * so a handful of them estimates nothing.
     */
    private static final int MIN_NULL_PAIRS = 256;

    /**
     * k for the exact-match screen that keeps homologous pairs OUT of the null. Two
     * windows sharing even one exact 25-mer might be related (4^25 makes that essentially
     * impossible by chance), so the pair is dropped rather than counted as chance.
     *
     * <p>This screen is why the null is drawn from REAL windows instead of shuffled ones.
     * Shuffling is the obvious way to get a homology-free null and it is valid for some
     * statistics — it was checked against the real thing for the old ratio statistic on
     * two genomes and tracked it closely. But it is NOT valid in general: D2* centres
     * counts on their compositional expectation, which is exactly what a shuffle flattens,
     * so its shuffled null collapses and its threshold comes out far too lax (6% false
     * positives, versus 0% on a real screened null). Screening real windows is honest for
     * any statistic, and costs one small k-mer set per calibration window.
     */
    private static final int NULL_SCREEN_K = 25;

    /**
     * Fewest distinct 4-mers a vector needs before a cosine against it means anything.
     * This replaces the old MIN_PAIRS_NNZ bound (9), which existed because the old
     * ratio similarity needed more than MIN_VECTOR_PAIRS dimension pairs to return
     * anything. Combined with the old rule that a slot needed a count &gt; 3 to survive
     * at all, it silently excluded 40% of blocks on the benchmark genome before any
     * comparison ran — the single largest cause of the mode's blindness. Cosine reads
     * raw counts, so the only real requirement is that the vector not be degenerate.
     */
    private static final int MIN_COSINE_NNZ = 8;

    /**
     * Fixed PRNG seed for the calibration draw. The floors, and therefore the cluster
     * set, must not depend on run-to-run randomness: with this fixed, the whole mode
     * stays deterministic and thread-count-independent, exactly as before.
     */
    private static final long NULL_SEED = 0x9E3779B97F4A7C15L;

    /** Shortest length class; shorter blocks all fall in class 0 and share its floor. */
    private static final int PROFILE_CLASS_BASE = 128;

    /** Number of sqrt(2)-spaced length classes: 128 … ~46 kb. */
    private static final int PROFILE_CLASS_COUNT = 18;

    private static final double LN2 = Math.log(2);

    /**
     * A seed is scanned window-by-window only once it is at least this many times
     * longer than the candidate. Below the ratio the whole-seed vector already IS a
     * scale-matched comparison, so windowing would only add cost and multiple testing.
     */
    private static final double WINDOW_MIN_RATIO = 1.5;

    /**
     * How much of the gap between the measured chance level and 1.0 to add to the floor,
     * per doubling of the number of windows searched — a multiple-comparisons charge for
     * keeping the best of W window comparisons.
     *
     * <p>Measured to be ZERO. It looks like it should be positive: searching W windows is
     * W chances at a chance match. But {@link #calibrateProfileFloors} already sets the
     * floor at the MAXIMUM over ~18k screened null pairs (~p99.995), which leaves so
     * little room above it that the extra charge buys almost no precision and costs a
     * great deal of recall. On MF782455, sweeping this constant end-to-end:
     *
     * <pre>
     *   penalty   blocks   families   verified   real members
     *     0.0       525       76        92.7%        416
     *     0.25      437      127        93.2%        289
     *     0.5       409      116        87.4%        256
     *     1.0       399      113        87.1%        249
     * </pre>
     *
     * 0.25 buys half a point of precision for 127 real members — a bad trade, and
     * Smith-Waterman put both operating points at 98% real regardless. So it is set to 0.
     *
     * <p>The mechanism is kept rather than deleted because the measurement has a known
     * boundary: the benchmark's longest block is ~6 kb, so W never exceeded ~46 here. A
     * run whose blocks are orders of magnitude longer would search far more windows per
     * seed, and if such a run ever shows chance absorptions this is the knob to raise —
     * with the sweep above as the guide to what it costs.
     */
    private static final double WINDOW_PENALTY_PER_DOUBLING = 0.0;

    /**
     * Minimum candidate span before a per-seed similarity scan is run in
     * parallel.
     */
    private static final int PARALLEL_SCAN_THRESHOLD = 256;

    /**
     * Canonical 4-mer palette: 238 entries — all 256 4-mers over {a,c,g,t} except
     * 18 low-complexity ones, which are deliberately absent so that microsatellites
     * cannot drive the profile (SSRs are detected separately). The excluded set is
     * every 4-mer of period 1 or 2 — the four homopolymers (aaaa, cccc, gggg, tttt)
     * and the twelve dinucleotide repeats (acac/gtgt, caca/tgtg, atat, cgcg, gcgc,
     * tata, agag/ctct, gaga/tctc) — plus the two AT-only palindromes aatt and ttaa.
     *
     * <p>A code with no palette entry maps to -1 in FWD_INDEX/RC_INDEX and is simply
     * not counted, so the exclusions are a filter, not a remapping.
     *
     * <p>The set is closed under reverse complement, which matters: it is what keeps
     * the forward and reverse profiles of a window mirror images of each other. It
     * was NOT closed before — ctct and gaga were dropped while their reverse
     * complements agag and tctc were kept, leaving the (AG)n/(CT)n family the only
     * dinucleotide family half-suppressed and the two strands measuring subtly
     * different things. That was corrected by also dropping agag and tctc.
     *
     * <p>aatt and ttaa are the two exclusions the period rule does not explain (both
     * have period 4). They are excluded together so the AT-only palindromes are
     * treated symmetrically; previously only aatt was dropped.
     */
   private static final String[] KMERS = {"acgt", "agct", "tgca", "tcga", "catg", "ctag", "gatc", "gtac", "ccgg", "ggcc", "tttc", "gaaa", "attt", "aaat", "tttg", "caaa", "tgaa", "ttct", "aaca", "ttca", "agaa", "tgtt", "tctt", "ttga", "aaga", "tcaa", "ttgt", "aaac", "gttt", "acaa", "aaag", "tgat", "cttt", "atca", "tcat", "atga", "catt", "acat", "attg", "caat", "atgt", "aatg", "attc", "gaat", "gatt", "aatc", "tatt", "aata", "atct", "ttcc", "agtt", "aact", "gaag", "cttc", "agat", "aagt", "ggaa", "actt", "ccaa", "tcca", "ttgg", "tgga", "ttgc", "gcaa", "gatg", "catc", "tcct", "caac", "gttg", "agga", "caag", "tgta", "aagg", "taca", "agca", "tatc", "cctt", "gata", "tgct", "acca", "atgg", "tggt", "atac", "gtat", "cttg", "ccat", "cata", "gtga", "tatg", "tctg", "gctt", "ttat", "gcat", "gttc", "ggat", "ataa", "atcc", "aagc", "ggtt", "gaac", "caga", "aacc", "atgc", "ttac", "gtaa", "ctga", "tcag", "ctgt", "atag", "ctat", "agta", "aggt", "tact", "ctcc", "atta", "taat", "tgtc", "tcac", "tgag", "acct", "acag", "ctca", "ggag", "gtca", "agtg", "taga", "tcta", "cctc", "taaa", "tgac", "ttta", "gaca", "gagg", "actc", "gagt", "actg", "tggc", "cagt", "taac", "gtta", "gtgg", "ggga", "tagt", "ccac", "gctg", "acta", "cacc", "cagc", "tccc", "ctgc", "gcag", "ggta", "gtct", "tacc", "ccag", "cact", "gtag", "ctgg", "ggtg", "agac", "tgcc", "tggg", "agtc", "cagg", "gact", "ggca", "cctg", "aggg", "taag", "ctaa", "gcca", "tagc", "ccct", "ctta", "gctc", "gcta", "ttag", "gagc", "ggct", "ctac", "ttcg", "gcac", "tagg", "cgaa", "ccta", "gtgc", "gggt", "tcgt", "aggc", "cgtt", "acga", "gcct", "aacg", "atcg", "cgat", "gacc", "ggtc", "ggac", "accc", "ccca", "gtcc", "ctcg", "agcc", "ccga", "cgag", "acgg", "tcgc", "cgct", "cgga", "tccg", "cgtc", "cacg", "gacg", "tcgg", "agcg", "gcga", "cgcc", "gccc", "accg", "cgta", "gccg", "tacg", "cgtg", "ccgc", "cgac", "ggcg", "gggc", "cggt", "cggc", "gtcg", "ccgt", "cgca", "tgcg", "gcgg", "cggg", "acgc", "gcgt", "cccg"};

    /**
     * Pre-built immutable lookup: kmer string → vector index.
     */
    private static final Map<String, Integer> KMER_INDEX;
    private static final int NUM_KMERS;

    /**
     * Encoded-kmer lookup tables (replace the per-position String/HashMap
     * path). A 4-mer over {a,c,g,t} is packed into 8 bits (2 bits per base),
     * giving a direct index in [0,256). FWD_INDEX maps that code to the forward
     * vector slot (or -1), RC_INDEX maps it to the reverse-complement vector
     * slot (or -1), and BASE_CODE maps an ASCII base to its 2-bit code (or -1
     * for any non-acgt char). The reverse-complement table is built once via
     * {@code Dna.ComplementDNA2} so its exact semantics are preserved without
     * per-position string allocation.
     */
    private static final int[] FWD_INDEX = new int[256];
    private static final int[] RC_INDEX = new int[256];
    private static final int[] BASE_CODE = new int[128];

    static {
        Map<String, Integer> map = new HashMap<>(KMERS.length * 2);
        for (int i = 0; i < KMERS.length; i++) {
            map.put(KMERS[i], i);
        }
        KMER_INDEX = Collections.unmodifiableMap(map);
        NUM_KMERS = KMER_INDEX.size();

        Arrays.fill(BASE_CODE, -1);
        BASE_CODE['a'] = 0;
        BASE_CODE['c'] = 1;
        BASE_CODE['g'] = 2;
        BASE_CODE['t'] = 3;

        Arrays.fill(FWD_INDEX, -1);
        Arrays.fill(RC_INDEX, -1);
        final char[] bases = {'a', 'c', 'g', 't'};
        for (int b0 = 0; b0 < 4; b0++) {
            for (int b1 = 0; b1 < 4; b1++) {
                for (int b2 = 0; b2 < 4; b2++) {
                    for (int b3 = 0; b3 < 4; b3++) {
                        int code = (b0 << 6) | (b1 << 4) | (b2 << 2) | b3;
                        String s = new String(new char[]{bases[b0], bases[b1], bases[b2], bases[b3]});
                        Integer f = KMER_INDEX.get(s);
                        if (f != null) {
                            FWD_INDEX[code] = f;
                        }
                        Integer r = KMER_INDEX.get(Dna.ComplementDNA2(s));
                        if (r != null) {
                            RC_INDEX[code] = r;
                        }
                    }
                }
            }
        }
    }

    /**
     * Per-thread scratch: the per-row count accumulators used only while building
     * vectors. Cosine needs no comparison scratch — it merges two CSR rows in place.
     */
    private static final ThreadLocal<int[]> DENSE_FWD = ThreadLocal.withInitial(() -> new int[NUM_KMERS]);
    private static final ThreadLocal<int[]> DENSE_REV = ThreadLocal.withInitial(() -> new int[NUM_KMERS]);

    private static final int[] EMPTY = new int[0];

    // ── Instance state ──────────────────────────────────────────────────────────
    private long[][] sequences;          // [i] = {globalOffset, length}, sorted by length desc
    private final int[] clusterIds;      // cluster assignment per sequence
    private int[] refClusterIds;         // cluster → reference-sequence mapping
    private int clusterCount;
    private final SeqStore store;        // segmented view over the combined data
    private final String[] referenceSeqs;
    private final int numReferenceSeqs;

    // Sequence vectors in CSR form (built once, retained for reference matching).
    // Row j occupies [fwdPtr[j], fwdPtr[j+1]); entries are (index, count), index-ascending.
    // Values are RAW 4-mer counts: see cosineMerge for why the old thresholded,
    // reciprocal-frequency encoding was dropped.
    private int[] fwdIdx, fwdVal, fwdPtr;
    private int[] revIdx, revVal, revPtr;
    private double[] fwdNorm, revNorm;

    // Reference vectors in CSR form.
    private int[] refFwdIdx, refFwdVal, refFwdPtr;
    private int[] refRevIdx, refRevVal, refRevPtr;
    private double[] refFwdNorm, refRevNorm;

    /**
     * Per length class, the smallest profile score that is not a chance match at that
     * scale, calibrated from this run's own data. Integer.MAX_VALUE marks a class that
     * was never calibrated (no block of that length exists), so nothing can pass it.
     * See {@link #calibrateProfileFloors}.
     */
    private double[] profileNullFloor;

    // Containment backend state (only built when useContainment): per block, the
    // sorted distinct forward k-mers. Storing the forward set alone, rather than one
    // canonical min(fwd, rev) set, is what lets a match recover its orientation — the
    // opposite-strand set is derived per seed by reverseComplementSet, not kept here
    // (see scanSeedContainment). Cost ~8 bytes per distinct k-mer per block.
    private long[][] cFwd;
    private int containK;   // the exact k the containment backend uses this run

    // ── Clustering measure ──────────────────────────────────────────────────────
    /** Cluster by the symmetric 4-mer ratio profile (whole-block composition). This is the default. */
    public static final int MODE_PROFILE = 0;
    /** Cluster by asymmetric k-mer containment. */
    public static final int MODE_CONTAIN = 1;

    // ── Constructor ─────────────────────────────────────────────────────────────
    /**
     * @param store segmented view over the concatenated sequence data
     * @param refSeqs reference sequences for cluster annotation (nullable)
     * @param offsets flat array of [offset0, length0, offset1, length1, …]
     * (global long coords)
     * @param useParallelClustering if true, uses the parallel clustering
     * strategy
     */
    public SequencesClustering(SeqStore store, String[] refSeqs, long[] offsets, boolean useParallelClustering) {
        this(store, refSeqs, offsets, useParallelClustering, MODE_PROFILE);
    }

    /**
     * @param useContainment kept for callers written against the old two-measure
     * API; maps to {@link #MODE_CONTAIN} / {@link #MODE_PROFILE}.
     */
    public SequencesClustering(SeqStore store, String[] refSeqs, long[] offsets,
            boolean useParallelClustering, boolean useContainment) {
        this(store, refSeqs, offsets, useParallelClustering,
                useContainment ? MODE_CONTAIN : MODE_PROFILE);
    }

    /**
     * @param mode which measure clusters the blocks: {@link #MODE_PROFILE} or
     * {@link #MODE_CONTAIN}. The output contract — signed cluster ids, longest-first
     * offsets, exclusive count, refClusterIds shape — is identical whichever is chosen.
     */
    public SequencesClustering(SeqStore store, String[] refSeqs, long[] offsets,
            boolean useParallelClustering, int mode) {
        this.store = store;

        if (refSeqs != null) {
            this.referenceSeqs = refSeqs;
            this.numReferenceSeqs = refSeqs.length;
        } else {
            this.referenceSeqs = new String[0];
            this.numReferenceSeqs = 0;
        }

        int numSeqs = offsets.length / 2;
        if (numSeqs < 1) {
            this.clusterIds = new int[0];
            this.refClusterIds = new int[0];
            return;
        }

        this.sequences = parseAndSortByLengthDesc(offsets, numSeqs);
        this.clusterIds = new int[numSeqs];

        this.clusterCount = (mode == MODE_CONTAIN)
                ? clusterContainment(numSeqs)
                : useParallelClustering ? clusterParallel(numSeqs) : clusterSequential(numSeqs);

        this.refClusterIds = new int[clusterCount];
        if (numReferenceSeqs > 0) {
            // Containment mode annotates off the containment k-mer index (built for every
            // block during clustering); the profile mode uses the forward CSR vectors.
            this.refClusterIds = (mode == MODE_PROFILE)
                    ? matchClustersToReferences(numSeqs, clusterCount)
                    : matchClustersToReferencesContainment(numSeqs, clusterCount);
        }
    }

    // ── Public API ──────────────────────────────────────────────────────────────
    public int getClusterCount() {
        return clusterCount;
    }

    public int[] getClusterIds() {
        return clusterIds;
    }

    public int[] getReferenceIds() {
        return refClusterIds;
    }

    public long[][] getSequenceOffsets() {
        return sequences;
    }

    // ── Parsing ─────────────────────────────────────────────────────────────────
    private static long[][] parseAndSortByLengthDesc(long[] offsets, int numSeqs) {
        long[][] seqs = new long[numSeqs][2];
        for (int j = 0; j < numSeqs; j++) {
            seqs[j][0] = offsets[j * 2];
            seqs[j][1] = offsets[j * 2 + 1];
        }
        Arrays.sort(seqs, (a, b) -> Long.compare(b[1], a[1]));
        return seqs;
    }

    // ── K-mer vector construction (CSR) ─────────────────────────────────────────
    /**
     * Builds forward and reverse-complement 4-mer count vectors for every sequence and
     * stores them in the instance CSR arrays, together with each row's Euclidean norm so
     * that {@link #cosineMerge} is a single merge pass.
     *
     * <p>Values are RAW COUNTS. Every slot the block actually contains is kept — there is
     * no minimum-occurrence rule, because cosine reads the whole vector and the pattern of
     * rare and absent 4-mers is the main thing that distinguishes two short blocks.
     */
    private void buildKmerVectors(int numSeqs) {
        final int[][] tFidx = new int[numSeqs][];
        final int[][] tFval = new int[numSeqs][];
        final int[][] tRidx = new int[numSeqs][];
        final int[][] tRval = new int[numSeqs][];
        final long[][] seqData = sequences;
        final SeqStore st = store;

        IntStream.range(0, numSeqs).parallel().forEach(j -> {
            long offset = seqData[j][0];
            int length = (int) seqData[j][1];   // a repeat block always fits in an int

            int[] fd = DENSE_FWD.get();
            int[] rd = DENSE_REV.get();
            Arrays.fill(fd, 0);
            Arrays.fill(rd, 0);

            accumulateKmers(st, offset, length, fd, rd);

            int[][] f = compactRow(fd);
            int[][] r = compactRow(rd);
            tFidx[j] = f[0];
            tFval[j] = f[1];
            tRidx[j] = r[0];
            tRval[j] = r[1];
        });

        int[][] fwd = flatten(tFidx, tFval, numSeqs);
        int[][] rev = flatten(tRidx, tRval, numSeqs);
        fwdIdx = fwd[0];
        fwdVal = fwd[1];
        fwdPtr = fwd[2];
        revIdx = rev[0];
        revVal = rev[1];
        revPtr = rev[2];
        fwdNorm = rowNorms(fwdVal, fwdPtr, numSeqs);
        revNorm = rowNorms(revVal, revPtr, numSeqs);
    }

    /** Euclidean norm of every CSR row, precomputed once so cosine stays a single pass. */
    private static double[] rowNorms(int[] val, int[] ptr, int rows) {
        double[] n = new double[rows];
        for (int j = 0; j < rows; j++) {
            long acc = 0;
            for (int t = ptr[j]; t < ptr[j + 1]; t++) {
                acc += (long) val[t] * val[t];
            }
            n[j] = Math.sqrt((double) acc);
        }
        return n;
    }

    /**
     * Compacts a normalised dense row into ascending (index, value) arrays,
     * value &gt; 0 only.
     */
    private static int[][] compactRow(int[] dense) {
        int m = 0;
        for (int k = 0; k < NUM_KMERS; k++) {
            if (dense[k] > 0) {
                m++;
            }
        }
        if (m == 0) {
            return new int[][]{EMPTY, EMPTY};
        }
        int[] idx = new int[m];
        int[] val = new int[m];
        int p = 0;
        for (int k = 0; k < NUM_KMERS; k++) {
            int v = dense[k];
            if (v > 0) {
                idx[p] = k;
                val[p] = v;
                p++;
            }
        }
        return new int[][]{idx, val};
    }

    /**
     * Concatenates per-row sparse arrays into flat CSR {idx, val, ptr}.
     */
    private static int[][] flatten(int[][] rowIdx, int[][] rowVal, int numRows) {
        int[] ptr = new int[numRows + 1];
        for (int j = 0; j < numRows; j++) {
            ptr[j + 1] = ptr[j] + rowIdx[j].length;
        }
        int total = ptr[numRows];
        int[] idx = new int[total];
        int[] val = new int[total];
        for (int j = 0; j < numRows; j++) {
            int len = rowIdx[j].length;
            if (len > 0) {
                System.arraycopy(rowIdx[j], 0, idx, ptr[j], len);
                System.arraycopy(rowVal[j], 0, val, ptr[j], len);
            }
        }
        return new int[][]{idx, val, ptr};
    }

    /**
     * Accumulates forward and reverse-complement 4-mer counts over the global
     * window {@code [offset, offset+length)} of {@code store}. A repeat block
     * lies inside a single segment, so the common path resolves that segment
     * once and delegates to the String loop with int-local coordinates; the
     * rare spanning case falls back to per-base global access for correctness.
     */
    private static void accumulateKmers(SeqStore store, long offset, int length, int[] fwdRow, int[] revRow) {
        if (length < KMER_SIZE) {
            return;
        }
        int seg = store.segmentOf(offset);
        long segStart = store.segmentStart(seg);
        String s = store.segment(seg);
        long local = offset - segStart;

        if (local + length <= s.length()) {                       // within one segment (the norm)
            accumulateKmers(s, (int) local, length, fwdRow, revRow);
        } else {                                                  // rare: block crosses a boundary
            accumulateKmersSpanning(store, offset, length, fwdRow, revRow);
        }
    }

    /**
     * Core rolling-code accumulator over {@code seq[offset, offset+length)}
     * using direct table lookups — no per-position String, HashMap or boxing. A
     * non-acgt character resets the rolling window (the kmer is skipped). Pass
     * {@code revRow == null} for forward-only.
     */
    private static void accumulateKmers(String seq, int offset, int length, int[] fwdRow, int[] revRow) {
        int code = 0;
        int valid = 0;
        final int end = offset + length;
        for (int p = offset; p < end; p++) {
            char ch = seq.charAt(p);
            int b = (ch < 128) ? BASE_CODE[ch] : -1;
            if (b < 0) {
                valid = 0;                       // reset window on a non-acgt base
                continue;
            }
            code = ((code << 2) | b) & 0xFF;
            if (++valid >= KMER_SIZE) {
                int fi = FWD_INDEX[code];
                if (fi >= 0) {
                    fwdRow[fi]++;
                }
                if (revRow != null) {
                    int ri = RC_INDEX[code];
                    if (ri >= 0) {
                        revRow[ri]++;
                    }
                }
            }
        }
    }

    /**
     * Boundary-spanning fallback: identical semantics via global per-base
     * access.
     */
    private static void accumulateKmersSpanning(SeqStore store, long offset, int length, int[] fwdRow, int[] revRow) {
        int code = 0;
        int valid = 0;
        final long end = offset + length;
        for (long p = offset; p < end; p++) {
            char ch = store.charAt(p);
            int b = (ch < 128) ? BASE_CODE[ch] : -1;
            if (b < 0) {
                valid = 0;
                continue;
            }
            code = ((code << 2) | b) & 0xFF;
            if (++valid >= KMER_SIZE) {
                int fi = FWD_INDEX[code];
                if (fi >= 0) {
                    fwdRow[fi]++;
                }
                if (revRow != null) {
                    int ri = RC_INDEX[code];
                    if (ri >= 0) {
                        revRow[ri]++;
                    }
                }
            }
        }
    }

    /**
     * Builds CSR k-mer vectors for the reference sequences.
     */
    private void buildReferenceKmerVectors() {
        final int n = numReferenceSeqs;
        final int[][] tFidx = new int[n][];
        final int[][] tFval = new int[n][];
        final int[][] tRidx = new int[n][];
        final int[][] tRval = new int[n][];
        final String[] refs = referenceSeqs;

        IntStream.range(0, n).parallel().forEach(j -> {
            int length = refs[j].length();
            if (length <= MIN_REF_LENGTH) {
                tFidx[j] = EMPTY;
                tFval[j] = EMPTY;
                tRidx[j] = EMPTY;
                tRval[j] = EMPTY;
                return;
            }
            int[] fd = DENSE_FWD.get();
            int[] rd = DENSE_REV.get();
            Arrays.fill(fd, 0);
            Arrays.fill(rd, 0);

            accumulateKmers(refs[j], 0, length, fd, rd);

            int[][] f = compactRow(fd);
            int[][] r = compactRow(rd);
            tFidx[j] = f[0];
            tFval[j] = f[1];
            tRidx[j] = r[0];
            tRval[j] = r[1];
        });

        int[][] fwd = flatten(tFidx, tFval, n);
        int[][] rev = flatten(tRidx, tRval, n);
        refFwdIdx = fwd[0];
        refFwdVal = fwd[1];
        refFwdPtr = fwd[2];
        refRevIdx = rev[0];
        refRevVal = rev[1];
        refRevPtr = rev[2];
        refFwdNorm = rowNorms(refFwdVal, refFwdPtr, n);
        refRevNorm = rowNorms(refRevVal, refRevPtr, n);
    }

    /**
     * Cosine similarity of two raw 4-mer count vectors held as CSR rows, in [0, 1].
     *
     * <h3>Why this replaced the ratio-of-ratios similarity</h3>
     * The old {@code similarityCore} compared, for every pair of dimensions (k,h), how
     * alike the two vectors express the ratio between k and h — but ONLY over the
     * dimensions both vectors share. Two things follow, and both were measured to be
     * severe:
     *
     * <ul>
     * <li><b>It throws away the support pattern.</b> Which 4-mers are ABSENT is the
     * dominant signal for a short block: a 200 bp block has ~197 4-mers spread over a
     * 238-slot palette, so its vector is mostly zeros and ones, and two homologous 200 bp
     * blocks agree slot-by-slot including the zeros. Intersecting the supports and
     * comparing ratios inside the intersection discards exactly that. Cosine reads the
     * whole vector.</li>
     * <li><b>It needed a count &gt; 3 to see a slot at all</b> (the old
     * {@code normalizeVector}'s MIN_KMER_HITS), which is why 40% of blocks on the
     * benchmark genome produced fewer than the 9 non-zeros it required and were skipped
     * before any comparison ran. Raw counts have no such floor.</li>
     * </ul>
     *
     * Measured head-to-head on MF782455, each statistic judged at ITS OWN zero-chance
     * threshold (the maximum over a null of real, verified-non-homologous window pairs at
     * the same scale), over 214 verified-homologous and 427 verified-non-homologous pairs:
     *
     * <pre>
     *   statistic          recall     false positives
     *   ratio (old)          30%           1.0%
     *   cosine               77%           0.2%
     *   D2* (centred)        62%           0.0%
     * </pre>
     *
     * Cosine is both more sensitive and more specific, so this is not a recall/precision
     * trade — it is strictly better on both. It is also cheaper: O(nnz) on one merge pass
     * versus the old O(m^2) double loop over the shared support, which is why the mode got
     * faster rather than slower despite now scanning windows.
     *
     * <p>D2* (counts centred on their expectation under each window's own base composition)
     * scored better than the old statistic too, and is the more principled measure in the
     * alignment-free literature, but it lost to cosine here and carries a trap worth
     * recording: its null CANNOT be estimated by shuffling, because a shuffle flattens the
     * very compositional deviation D2* centres on. A shuffled null put its threshold far
     * too low and produced 6% false positives; only a real, homology-screened null gives it
     * an honest one. Cosine's null does not have that sensitivity.
     *
     * <p>Scale-invariance is a free bonus: cosine divides by both norms, so a window that
     * is a little longer than the candidate does not bias the score. That removes the
     * length artefact documented on the old {@code normalizeVector} — the {@code +1} pseudocount
     * in count space, which shrank short blocks toward uniformity — rather than working
     * around it.
     */
    private static double cosineMerge(int[] idxA, int[] valA, int fromA, int toA, double normA,
            int[] idxB, int[] valB, int fromB, int toB, double normB) {
        if (normA <= 0 || normB <= 0) {
            return 0.0;
        }
        long dot = 0;
        int pa = fromA, pb = fromB;
        while (pa < toA && pb < toB) {
            int ia = idxA[pa], ib = idxB[pb];
            if (ia == ib) {
                dot += (long) valA[pa] * valB[pb];
                pa++;
                pb++;
            } else if (ia < ib) {
                pa++;
            } else {
                pb++;
            }
        }
        return dot / (normA * normB);
    }

    // ── Profile chance floor: length classes ────────────────────────────────────
    /**
     * The sqrt(2)-spaced length class of a block. Classes are geometric because the
     * profile's noise scales with the LOG of length, not with length: the gap between
     * 200 bp and 400 bp matters enormously, the gap between 20 kb and 40 kb not at all.
     */
    private static int lengthClass(int len) {
        if (len <= PROFILE_CLASS_BASE) {
            return 0;
        }
        int c = (int) Math.floor(2.0 * (Math.log(len / (double) PROFILE_CLASS_BASE) / LN2));
        return Math.max(0, Math.min(PROFILE_CLASS_COUNT - 1, c));
    }

    /** Representative length of a class: the shortest block that falls in it. */
    private static int classLength(int c) {
        return (int) Math.round(PROFILE_CLASS_BASE * Math.pow(2, c / 2.0));
    }

    /**
     * The cosine a candidate in length class {@code cls} must EXCEED to be believed, when
     * the decision was taken as the best of {@code nWindows} window comparisons. Returns
     * a value &gt; 1 for a scale where chance could not be measured, which no cosine can
     * clear.
     *
     * <p>Takes the CLASS, not a length, deliberately. An earlier version took the window
     * length and recovered the class with {@code lengthClass(classLength(c))}, which is
     * not the identity: {@code classLength} rounds to an int, and for every odd class the
     * rounded value falls a hair below the class boundary, so the round trip returned
     * {@code c-1} — 7 of the 18 classes silently read a floor calibrated for a different
     * scale. Passing the class removes the hazard rather than papering over it with an
     * epsilon.
     */
    private double profileFloor(int cls, int nWindows) {
        double base = profileNullFloor[cls];
        if (base >= 1.0) {
            return Double.MAX_VALUE;
        }
        // Searching W windows and keeping the best is W chances at a chance match. Charge
        // for it by moving the floor a fraction of the way from the measured chance level
        // to 1, geometrically in the number of windows.
        double penalty = (nWindows <= 1) ? 0.0
                : WINDOW_PENALTY_PER_DOUBLING * (Math.log(nWindows) / LN2) * (1.0 - base);
        return Math.min(0.999, base + penalty);
    }

    // ── Profile chance floor: calibration ───────────────────────────────────────
    /**
     * Calibrates, per length class, the highest profile score that mere chance produces
     * at that scale in THIS run's data — the profile's answer to
     * {@link #containmentFloor}.
     *
     * <h3>Why this exists</h3>
     * The containment measure has always had a chance control ({@link #containmentFloor},
     * 5 sigma clear of the null). The profile had NONE: it compared every pair, at every
     * length, against one flat threshold (MIN_SIMILARITY, 70). That is wrong in both
     * directions at once, and both were measured on MF782455 (1.39 Mb, 765 blocks) against ground
     * truth independent of either measure — pairs sharing >= 3 exact canonical 25-mers
     * (chance ~ 4^-25), cross-checked by Smith-Waterman against a shuffled null:
     *
     * <ul>
     * <li><b>Far too lax when blocks are short.</b> The null score of two NON-homologous
     * real 250 bp windows has p99 = 65 and p99.9 = 83 — so ~1% of unrelated 250 bp pairs
     * clear a flat bar of 70 per comparison, and the greedy loop offers every candidate
     * to ~90 seeds in turn. End to end that dominated the output: of the 205 blocks
     * -vector clustered, only 53% showed ANY homology to their family seed; Smith-Waterman
     * put the noise fraction at 45%, and those pairs scored at or BELOW their own shuffled
     * null (score/null 0.73–0.96), so they were not diverged homology — they were nothing.</li>
     * <li><b>Far too strict when blocks are long.</b> At 2 kb the null never once exceeded
     * -1 across 20000 draws. A flat bar of 70 there discards real homology by a huge
     * margin: blocks at 1313571 and 1320896 (5354 bp and 4415 bp) share 2028 exact 31-mers
     * and a 600 bp exact common substring — unambiguous homology — and score 62. They were
     * rejected. Both are among the 11 of the 12 longest blocks in the genome that -vector
     * left unclustered while -contain placed every one of them in a family.</li>
     * </ul>
     *
     * No flat threshold can fix this, and the constant was not simply mistuned: sweeping
     * the old MIN_SIMILARITY end-to-end trades noise for signal at nearly 1:1 (70 → 205
     * blocks at 53% verified; 80 → 102 at 75%; 90 → 29 at 100%). The bar has to depend on
     * length, because the null does.
     *
     * <h3>Why it is calibrated per run rather than tabulated</h3>
     * The null is a property of the GENOME, not of the algorithm. At 600 bp its p99 is 45
     * on MF782455 but 85 on E. coli CP006627 — a compositionally homogeneous replicon
     * makes unrelated windows resemble each other. A constant table fitted to one genome
     * would be badly wrong on the next, so the floor is measured from the data in hand,
     * exactly as containmentFloor derives its chance level from the actual set sizes.
     *
     * <h3>Method</h3>
     * For each length class that some block occupies, draw {@link #NULL_WINDOWS} windows
     * of that length from the store and score them all against each other, keeping only
     * pairs that share no exact {@link #NULL_SCREEN_K}-mer — i.e. pairs that are certainly
     * not homologous. The floor is the maximum cosine any such pair reaches.
     *
     * <p>The maximum of ~18k screened pairs is about the p99.995 point, which is
     * deliberately matched to the size of a run: a candidate is offered to every seed in
     * turn, so a run makes ~1e4-1e5 comparisons and a per-comparison false-positive rate
     * of ~5e-5 keeps expected chance absorptions to a handful rather than the hundreds a
     * flat bar admits. Drawing 192 windows and pairing them all is what makes this
     * affordable: the cost is 192 vector builds, not 18k.
     *
     * <p>Cost is a few hundred vector builds and ~18k cosines per occupied class — a small
     * fraction of the clustering itself. The k-mer sets are built for the calibration
     * sample only (a few hundred windows) and dropped immediately, so -vector keeps the
     * small, flat memory profile that is its reason to exist; it does not acquire
     * -contain's per-block k-mer index.
     */
    private void calibrateProfileFloors(int numSeqs) {
        final double[] floors = new double[PROFILE_CLASS_COUNT];
        Arrays.fill(floors, 2.0);   // > 1: unusable until measured

        final boolean[] need = new boolean[PROFILE_CLASS_COUNT];
        for (int j = 0; j < numSeqs; j++) {
            need[lengthClass((int) sequences[j][1])] = true;
        }

        final long total = store.length();
        IntStream.range(0, PROFILE_CLASS_COUNT).parallel().forEach(c -> {
            if (!need[c]) {
                return;
            }
            final int wlen = classLength(c);
            if (wlen < KMER_SIZE || wlen >= total) {
                return;
            }
            floors[c] = nullFloorForLength(wlen, total, NULL_SEED + c);
        });

        this.profileNullFloor = floors;
    }

    /**
     * The highest cosine {@link #NULL_WINDOWS} non-homologous windows of {@code wlen}
     * produce against each other — the chance level at that scale.
     *
     * <p>The null must be conditioned on the SAME event the real test is: a comparison
     * only happens when both vectors clear {@link #MIN_COSINE_NNZ}, so only windows that
     * clear it may enter the null. That is why this draws until it has
     * {@link #NULL_WINDOWS} USABLE windows rather than drawing NULL_WINDOWS and keeping
     * whatever survives.
     * The distinction is not academic — it was measured, and getting it wrong is a
     * silent, one-directional failure. At 128 bp only 5 windows in 192 yield a usable
     * vector, so the naive version took its maximum over ~10 pairs: that is the p90 of
     * the null, not its tail, and it set the floor for the 128 bp class at 42 when the
     * conditional null runs far higher. Short blocks are the bulk of a typical run
     * (452 of 765 here), so this one bug would have re-admitted most of the noise the
     * floor exists to exclude.
     *
     * <p>Returns Integer.MAX_VALUE when the draw cap is hit without enough usable
     * windows. That correctly makes the class unusable: if chance cannot be measured at
     * a scale, nothing at that scale can be believed.
     */
    private double nullFloorForLength(int wlen, long total, long seed) {
        final Random rnd = new Random(seed);
        final int[][] idx = new int[NULL_WINDOWS][];
        final int[][] val = new int[NULL_WINDOWS][];
        final double[] nrm = new double[NULL_WINDOWS];
        final long[][] key = new long[NULL_WINDOWS][];

        final int maxDraws = NULL_WINDOWS * NULL_DRAW_LIMIT;
        int built = 0;
        for (int draw = 0; draw < maxDraws && built < NULL_WINDOWS; draw++) {
            long off = (long) (rnd.nextDouble() * (total - wlen));
            String w = store.substring(off, off + wlen);

            int[] dense = new int[NUM_KMERS];
            accumulateKmers(w, 0, wlen, dense, null);
            int[][] row = compactRow(dense);
            if (row[0].length < MIN_COSINE_NNZ) {
                continue;                     // would never be compared: not part of the null
            }
            idx[built] = row[0];
            val[built] = row[1];
            long acc = 0;
            for (int v : row[1]) {
                acc += (long) v * v;
            }
            nrm[built] = Math.sqrt((double) acc);
            key[built] = forwardSet(w, 0, wlen, NULL_SCREEN_K);
            built++;
        }
        if (built < MIN_NULL_WINDOWS) {
            return 2.0;   // chance is not measurable here: nothing at this scale can be believed
        }

        double best = 0.0;
        long pairs = 0;
        for (int i = 0; i < built; i++) {
            for (int j = i + 1; j < built; j++) {
                if (intersectionCount(key[i], key[j]) > 0) {
                    continue;                 // the two windows may be homologous: not part of the null
                }
                pairs++;
                double s = cosineMerge(idx[i], val[i], 0, idx[i].length, nrm[i],
                        idx[j], val[j], 0, idx[j].length, nrm[j]);
                if (s > best) {
                    best = s;
                }
            }
        }
        // A highly repetitive genome can leave almost every drawn pair sharing a k-mer, and
        // so screened out. Without this the floor would fall back to its 0.0 initialiser and
        // admit EVERYTHING at that scale — the failure would be silent and maximally wrong.
        // Too few surviving pairs means chance was not measured, so the scale is unusable.
        if (pairs < MIN_NULL_PAIRS) {
            return 2.0;
        }
        return best;
    }

    // ── Profile chance floor: scale-matched seed windows ────────────────────────
    /**
     * The seed's own profile, recomputed over each window of {@code wlen} bases as a DENSE
     * palette-slot vector (one {@code int[NUM_KMERS]} per window, NUM_KMERS = 238). Windows step by half their
     * length so that a homologous stretch cannot fall between two of them.
     *
     * <h3>Why windowing</h3>
     * The whole-block profile answers "do these two blocks have the same overall
     * composition?", but real repeat homology is PARTIAL — a short element copied from
     * part of a long one, or two elements sharing only some modules. Averaging the shared
     * part into the whole block dilutes it away before any comparison happens, which is
     * why the profile went blind to length-disparate pairs. Comparing the candidate
     * against the best window of the seed AT THE CANDIDATE'S OWN SCALE restores it, and
     * has the second virtue of making both sides the same length — so the length artefact
     * baked into the old {@code normalizeVector}'s pseudocount cancels, and one chance floor per
     * scale is exactly the right null.
     *
     * <p>Measured on MF782455 against 25-mer-verified homology: of pairs whose candidate
     * is long enough to have an estimable profile, the whole-block comparison finds 15%
     * and the windowed one 35%. Broken out by candidate length the effect is stark — in
     * the 200-300 bp band, 5% found before, 81% after. At every length ratio above 1.5x
     * the old comparison found essentially nothing (0-1%).
     *
     * <h3>Why dense</h3>
     * The dense form lets a candidate's sparse vector be dotted against a window by direct
     * indexing ({@link #dotDenseSparse}): a branch-free O(candidate nnz) pass rather than the
     * sorted merge {@link #cosineMerge} runs, which for the near-dense windows of a long seed
     * was the dominant cost of the scan (about 4-6x faster here). The score is byte-identical
     * — same integer dot, same norm.
     *
     * <p>Cost is O(seedLen) per class regardless of window count, because the windows tile
     * the seed: 2*seedLen/wlen windows of wlen bases each, built once per (seed, class) and
     * shared by every candidate in that class.
     */
    private int[][] buildWindowDense(long seedOff, int seedLen, int wlen) {
        final int step = Math.max(1, wlen / 2);
        final int nWin = (seedLen - wlen) / step + 1;
        final int[][] windows = new int[nWin][NUM_KMERS];
        for (int w = 0; w < nWin; w++) {
            accumulateKmers(store, seedOff + (long) w * step, wlen, windows[w], null);
        }
        return windows;
    }

    /** Scatters a sparse CSR row into a dense {@code int[NUM_KMERS]} vector. */
    private static int[] denseFromSparse(int[] idx, int[] val, int from, int to) {
        final int[] d = new int[NUM_KMERS];
        for (int k = from; k < to; k++) {
            d[idx[k]] = val[k];
        }
        return d;
    }

    /** Euclidean norm of a dense {@code int[NUM_KMERS]} vector. */
    private static double denseNorm(int[] dense) {
        long acc = 0;
        for (int v : dense) {
            acc += (long) v * v;
        }
        return Math.sqrt((double) acc);
    }

    /** Dot of a dense window against a sparse candidate row: O(candidate nnz), branch-free. */
    private static long dotDenseSparse(int[] dense, int[] idx, int[] val, int from, int to) {
        long dot = 0;
        for (int k = from; k < to; k++) {
            dot += (long) dense[idx[k]] * val[k];
        }
        return dot;
    }

    // ── Clustering (shared greedy loop, sequential or parallel scan) ─────────────
    private int clusterSequential(int numSeqs) {
        return cluster(numSeqs, /* parallel */ false);
    }

    private int clusterParallel(int numSeqs) {
        return cluster(numSeqs, /* parallel */ true);
    }

    /**
     * Greedy seed-and-absorb clustering. The outer loop is sequential and
     * deterministic (a sequence becomes a new seed only if no earlier, longer
     * seed absorbed it — an order-dependent decision). The per-seed candidate
     * scan is the parallel part. Candidates are pulled from a compacted active
     * set of still- unclustered indices, kept ascending, so the seed is always
     * {@code active[0]}.
     *
     * <p>
     * Ids are COMPACT (a rejected singleton's id is reused) and the returned
     * value is the EXCLUSIVE upper bound (maxAssignedId + 1), matching the
     * original.
     */
    private int cluster(int numSeqs, boolean parallel) {
        final int[] active = new int[numSeqs];
        for (int j = 0; j < numSeqs; j++) {
            active[j] = j;
        }
        return cluster(numSeqs, parallel, active, numSeqs, 2);
    }

    /**
     * {@link #cluster} over a GIVEN candidate set, numbering on from a given id.
     *
     * <p>{@code active[0..activeLen)} must hold block indices whose clusterIds are
     * still 0, ascending (so the seed, {@code active[0]}, is the longest of them —
     * the blocks are length-desc sorted). Blocks outside that set are never read or
     * written, so the pass can run over an arbitrary subset of blocks. {@code startCluster}
     * is the highest id already handed out, so the first cluster this pass forms gets
     * {@code startCluster + 1}.
     *
     * <p>Vectors are built for ALL blocks regardless: the CSR rows are indexed by
     * block, and the cost is bounded by the palette, not the candidate count.
     */
    private int cluster(int numSeqs, boolean parallel,
            int[] active, int activeLen, int startCluster) {
        buildKmerVectors(numSeqs);
        calibrateProfileFloors(numSeqs);

        final int[] assignBuf = new int[numSeqs];     // reused; one slot per candidate
        int nextCluster = startCluster;

        while (activeLen > 0) {
            int seed = active[0];

            nextCluster++;
            clusterIds[seed] = nextCluster;

            int matched = scanSeed(seed, active, activeLen, nextCluster, parallel, assignBuf);

            // Remove singleton clusters; reuse the id so numbering stays gap-free.
            if (matched == 0) {
                nextCluster--;
                clusterIds[seed] = 0;
            }

            // Compact the active set: drop the seed (slot 0) and everything absorbed.
            int w = 0;
            for (int t = 1; t < activeLen; t++) {
                int idx = active[t];
                if (clusterIds[idx] == 0) {
                    active[w++] = idx;
                }
            }
            activeLen = w;
        }
        return nextCluster + 1;
    }

    /**
     * Fused forward+reverse candidate scan for one seed. Each candidate
     * {@code active[1..activeLen)} is judged forward first; if it clears the
     * threshold it is assigned {@code +clusterId}, otherwise it is judged
     * reverse and, on success, assigned {@code -clusterId}. This is exactly the
     * original "all forward, then all remaining reverse" outcome (forward keeps
     * priority), in a single pass.
     *
     * <p>
     * The parallel stage only writes to {@code assignBuf} (one private slot per
     * candidate) and reads shared state that does not change during the scan,
     * so it is race-free; the cluster-id writes are applied afterwards on the
     * caller.
     */
    private int scanSeed(int seedIdx, int[] active, int activeLen, int clusterId,
            boolean parallel, int[] assignBuf) {
        final int from = 1, to = activeLen;
        if (to <= from) {
            return 0;
        }

        final int nnzSeed = fwdPtr[seedIdx + 1] - fwdPtr[seedIdx];
        if (nnzSeed < MIN_COSINE_NNZ) {
            return 0;   // seed vector degenerate: nothing to compare against
        }

        final int seedLen = (int) sequences[seedIdx][1];
        final long seedOff = sequences[seedIdx][0];
        final int posId = clusterId, negId = -clusterId;
        final int[] act = active;
        final int[] fI = fwdIdx, fV = fwdVal, fP = fwdPtr;
        final int[] rI = revIdx, rV = revVal, rP = revPtr;
        final double[] fN = fwdNorm, rN = revNorm;
        final int[] out = assignBuf;
        Arrays.fill(out, 0, to - from, 0);

        // Bucket the candidates by length class. Every candidate in a class is judged at
        // the same scale, against the same seed windows and the same calibrated floor, so
        // the windows are built once per class rather than once per candidate.
        final int[] clsOf = new int[to - from];
        final int[] clsCount = new int[PROFILE_CLASS_COUNT];
        for (int p = from; p < to; p++) {
            int c = lengthClass((int) sequences[act[p]][1]);
            clsOf[p - from] = c;
            clsCount[c]++;
        }

        for (int c = 0; c < PROFILE_CLASS_COUNT; c++) {
            if (clsCount[c] == 0) {
                continue;
            }
            final int wlen = classLength(c);
            final boolean windowed = seedLen >= WINDOW_MIN_RATIO * wlen;

            // Windows of the seed at this class's scale as DENSE palette-slot vectors, or the
            // seed's own vector when the two are already the same size (then the whole-block
            // comparison IS scale-matched). Dense windows let each candidate be scored by a
            // branch-free dot (dotDenseSparse) instead of the sorted merge — identical score.
            final int[][] windows = windowed
                    ? buildWindowDense(seedOff, seedLen, wlen)
                    : new int[][]{denseFromSparse(fI, fV, fP[seedIdx], fP[seedIdx + 1])};
            final int nWin = windows.length;
            final double[] wNorm = new double[nWin];
            for (int q = 0; q < nWin; q++) {
                wNorm[q] = denseNorm(windows[q]);
            }
            final double floor = profileFloor(c, nWin);
            if (floor > 1.0) {
                continue;   // chance is unmeasurable at this scale: nothing here can be believed
            }

            final int[] positions = new int[clsCount[c]];
            int w = 0;
            for (int p = from; p < to; p++) {
                if (clsOf[p - from] == c) {
                    positions[w++] = p;
                }
            }

            IntStream stream = IntStream.range(0, positions.length);
            if (parallel && positions.length >= PARALLEL_SCAN_THRESHOLD) {
                stream = stream.parallel();
            }
            stream.forEach(t -> {
                final int p = positions[t];
                final int j = act[p];
                int v = 0;
                // Forward keeps priority over reverse, as before: every window is tried
                // forward before any is tried reverse. dot/(wNorm*cNorm) reproduces
                // cosineMerge exactly, including its "norm <= 0 -> 0" guard.
                int fFrom = fP[j], fTo = fP[j + 1];
                if ((fTo - fFrom) >= MIN_COSINE_NNZ && fN[j] > 0) {
                    for (int q = 0; q < nWin && v == 0; q++) {
                        if (wNorm[q] > 0
                                && dotDenseSparse(windows[q], fI, fV, fFrom, fTo) / (wNorm[q] * fN[j]) > floor) {
                            v = posId;
                        }
                    }
                }
                if (v == 0) {
                    int rFrom = rP[j], rTo = rP[j + 1];
                    if ((rTo - rFrom) >= MIN_COSINE_NNZ && rN[j] > 0) {
                        for (int q = 0; q < nWin && v == 0; q++) {
                            if (wNorm[q] > 0
                                    && dotDenseSparse(windows[q], rI, rV, rFrom, rTo) / (wNorm[q] * rN[j]) > floor) {
                                v = negId;
                            }
                        }
                    }
                }
                out[p - from] = v;
            });
        }

        // Apply assignments single-threaded.
        int matched = 0;
        for (int p = from; p < to; p++) {
            int v = out[p - from];
            if (v != 0) {
                clusterIds[act[p]] = v;
                matched++;
            }
        }
        return matched;
    }

    // ── Reference matching ──────────────────────────────────────────────────────
    /**
     * Annotates each cluster with the first matching reference. Members are
     * grouped by cluster in a single pass (only positive ids — i.e. seeds and
     * forward members — exactly as the original, which compared {@code clusterIds[i] ==
     * clusterId}), and each member is matched using the forward vector already
     * built during clustering rather than recomputing it from the store.
     */
    private int[] matchClustersToReferences(int numSeqs, int numClusters) {
        buildReferenceKmerVectors();

        final int[] clusterToRef = new int[numClusters];

        // Group positive members by cluster id (ascending index order within a cluster).
        int[] cnt = new int[numClusters];
        for (int i = 0; i < numSeqs; i++) {
            int c = clusterIds[i];
            if (c >= 3 && c < numClusters) {
                cnt[c]++;
            }
        }
        int[] memPtr = new int[numClusters + 1];
        for (int c = 0; c < numClusters; c++) {
            memPtr[c + 1] = memPtr[c] + cnt[c];
        }
        int[] members = new int[memPtr[numClusters]];
        int[] fill = new int[numClusters];
        for (int i = 0; i < numSeqs; i++) {
            int c = clusterIds[i];
            if (c >= 3 && c < numClusters) {
                members[memPtr[c] + fill[c]++] = i;
            }
        }

        IntStream.range(3, numClusters).parallel().forEach(clusterId -> {
            int mStart = memPtr[clusterId], mEnd = memPtr[clusterId + 1];
            for (int mi = mStart; mi < mEnd; mi++) {
                int seq = members[mi];
                int sFrom = fwdPtr[seq], sTo = fwdPtr[seq + 1];
                if ((sTo - sFrom) < MIN_COSINE_NNZ) {
                    continue;   // can't match anything
                }
                int matchedRef = findMatchingReference(refFwdIdx, refFwdVal, refFwdPtr, refFwdNorm, seq, sFrom, sTo);
                if (matchedRef < 0) {
                    matchedRef = findMatchingReference(refRevIdx, refRevVal, refRevPtr, refRevNorm, seq, sFrom, sTo);
                }
                if (matchedRef >= 0) {
                    clusterToRef[clusterId] = matchedRef + 1; // 1-indexed
                    break;
                }
            }
        });

        return clusterToRef;
    }

    /**
     * Finds the first reference whose vector is sufficiently similar to the
     * sequence's forward vector (CSR row {@code [seqFrom, seqTo)} of
     * {@code fwd*}). Mirrors the original
     * {@code computeSimilarity(refVector, seqVector, false)}: a=ref value,
     * b=seq value, full (non-early-exit) computation.
     *
     * @return reference index (0-based), or -1 if no match
     */
    private int findMatchingReference(int[] refIdx, int[] refVal, int[] refPtr, double[] refNorm,
            int seq, int seqFrom, int seqTo) {
        final int seqLen = (int) sequences[seq][1];
        for (int r = 0; r < numReferenceSeqs; r++) {
            int rFrom = refPtr[r], rTo = refPtr[r + 1];
            if ((rTo - rFrom) < MIN_COSINE_NNZ) {
                continue;   // short/empty ref can't match
            }
            // A library entry and a block are rarely the same length, so the floor is taken
            // at the shorter of the two — the noisier side is what sets the chance level.
            int scale = Math.min(seqLen, referenceSeqs[r].length());
            double floor = profileFloor(lengthClass(scale), 1);
            if (floor > 1.0) {
                continue;
            }
            double similarity = cosineMerge(refIdx, refVal, rFrom, rTo, refNorm[r],
                    fwdIdx, fwdVal, seqFrom, seqTo, fwdNorm[seq]);
            if (similarity > floor) {
                return r;
            }
        }
        return -1;
    }

    // ── Containment backend (the -contain mode) ─────────────────────────────────
    /**
     * The k a containment run uses: large enough to be specific for the typical
     * block length, k = log4(L*99), forced odd, clamped to [11,24]. Ported from
     * GeneDistance's containmentK. L is the MEDIAN block length (not the maximum),
     * so a run mixing a few very long blocks with many short ones does not pick a
     * k so large that the short blocks yield no k-mers at all.
     */
    private static int containmentK(int length) {
        int k = (int) Math.ceil(Math.log(Math.max(1, length) * 99.0) / Math.log(4));
        return Math.max(CONTAIN_K_MIN, Math.min(CONTAIN_K_MAX, k | 1));
    }

    /**
     * Smallest containment ratio that counts as evidence rather than coincidence,
     * for a query of {@code nQuery} distinct k-mers tested against a container of
     * {@code nContainer} distinct k-mers.
     *
     * <p>Null model: the container's set covers a fraction {@code p = nContainer /
     * 4^k} of k-mer space, so each of the query's k-mers falls inside it by chance
     * with probability p. The shared count is then Binomial(nQuery, p), and
     * {@code c = shared / nQuery} has mean p and standard deviation
     * {@code sqrt(p(1-p)/nQuery)}. We demand that c clear the null mean by
     * {@link #CHANCE_SIGMA} standard deviations. Note that the query's size sets the
     * VARIANCE, not the mean — a bigger query gives a tighter estimate of the same
     * chance level.
     *
     * <p>This replaces the original {@code c > 10 * chance} with
     * {@code chance = max(nQuery, nContainer) / 4^k}, which was wrong twice over.
     * The mean of the null depends only on the CONTAINER's coverage, so taking a max
     * with the query's size had no meaning. Worse, "10 * p" is not a probability
     * statement and runs off the end of the scale: it exceeds 1 once the container
     * covers more than a tenth of k-mer space (nContainer > 4^k/10, i.e. ~419k
     * distinct 11-mers, reached by roughly 440kb of complex sequence), while
     * {@code c <= 1} by construction. Past that point EVERY candidate was rejected
     * unconditionally, however well it matched. Since {@link #clusterContainment}
     * seeds longest-first, the blocks most able to contain a short fragment were
     * exactly the ones that went inert — the failure landed on this mode's whole
     * reason to exist.
     *
     * <p>A deviation floor degrades gracefully instead. It sits far below
     * MIN_CONTAINMENT for ordinary seeds, so MIN_CONTAINMENT stays the binding
     * constraint and the decisions there are unchanged; it climbs as the container
     * saturates, demanding near-perfect containment from a seed that holds most of
     * k-mer space; and it forbids everything only when the container really does
     * hold every k-mer, where rejection is correct because such a set matches
     * anything.
     *
     * <p>Caveat: the binomial null assumes the query's k-mers are independent draws.
     * Overlapping k-mers share k-1 bases and real DNA is composition-skewed, so the
     * true variance is somewhat larger and the true p somewhat higher than uniform
     * sampling implies. CHANCE_SIGMA absorbs this only crudely; the floor is a sanity
     * bound against coincidence, not a calibrated p-value.
     */
    private static double containmentFloor(int nQuery, int nContainer, double pow4k) {
        final double p = nContainer / pow4k;
        if (p >= 1.0) {
            return 1.0;   // container spans k-mer space: no observation can be evidence
        }
        return p + CHANCE_SIGMA * Math.sqrt(p * (1.0 - p) / nQuery);
    }

    /** Number of k-mers two sorted, de-duplicated sets have in common (linear merge). */
    private static int intersectionCount(long[] a, long[] b) {
        int i = 0, j = 0, n = 0;
        while (i < a.length && j < b.length) {
            long x = a[i], y = b[j];
            if (x < y) {
                i++;
            } else if (x > y) {
                j++;
            } else {
                n++;
                i++;
                j++;
            }
        }
        return n;
    }

    private static long[] sortDedup(long[] buf, int m) {
        if (m == 0) {
            return EMPTY_L;
        }
        long[] r = Arrays.copyOf(buf, m);
        Arrays.sort(r);
        int d = 0;
        for (int i = 0; i < r.length; i++) {
            if (i == 0 || r[i] != r[i - 1]) {
                r[d++] = r[i];
            }
        }
        return (d == r.length) ? r : Arrays.copyOf(r, d);
    }

    /**
     * The distinct forward k-mers of a window, sorted. A match in the SAME
     * orientation is measured against this set; a match in the opposite orientation
     * is measured against {@link #reverseComplementSet} of it. With odd k no k-mer
     * equals its own reverse complement, so the two sets are disjoint and the
     * winning one gives the strand of a match.
     *
     * A non-acgt symbol resets the rolling window (the k-mer is skipped), exactly
     * as the 4-mer accumulator does. Codes are packed 2 bits/base into a long,
     * first base in the high bits. May be empty.
     */
    private static long[] forwardSet(String seq, int offset, int length, int k) {
        if (length < k) {
            return EMPTY_L;
        }
        final long mask = (1L << (2 * k)) - 1;
        final long[] fbuf = new long[length - k + 1];
        int m = 0;
        long fwd = 0;
        int valid = 0;
        final int end = offset + length;
        for (int p = offset; p < end; p++) {
            final char ch = seq.charAt(p);
            final int b = (ch < 128) ? BASE_CODE[ch] : -1;
            if (b < 0) {
                valid = 0;
                fwd = 0;
                continue;
            }
            fwd = ((fwd << 2) | b) & mask;
            if (valid < k) {
                valid++;
            }
            if (valid == k) {
                fbuf[m] = fwd;
                m++;
            }
        }
        return sortDedup(fbuf, m);
    }

    /**
     * {@link #forwardSet} for a block of the store: resolve the owning segment
     * once and read it with int-local coordinates (the norm — a repeat block lies
     * inside one segment), falling back to per-base global access only for the
     * rare block that spans a boundary. Mirrors {@link #accumulateKmers}.
     */
    private static long[] forwardSetStore(SeqStore store, long offset, int length, int k) {
        if (length < k) {
            return EMPTY_L;
        }
        int seg = store.segmentOf(offset);
        long segStart = store.segmentStart(seg);
        String s = store.segment(seg);
        long local = offset - segStart;
        if (local + length <= s.length()) {
            return forwardSet(s, (int) local, length, k);
        }
        return forwardSetSpanning(store, offset, length, k);
    }

    private static long[] forwardSetSpanning(SeqStore store, long offset, int length, int k) {
        final long mask = (1L << (2 * k)) - 1;
        final long[] fbuf = new long[length - k + 1];
        int m = 0;
        long fwd = 0;
        int valid = 0;
        final long end = offset + length;
        for (long p = offset; p < end; p++) {
            final char ch = store.charAt(p);
            final int b = (ch < 128) ? BASE_CODE[ch] : -1;
            if (b < 0) {
                valid = 0;
                fwd = 0;
                continue;
            }
            fwd = ((fwd << 2) | b) & mask;
            if (valid < k) {
                valid++;
            }
            if (valid == k) {
                fbuf[m] = fwd;
                m++;
            }
        }
        return sortDedup(fbuf, m);
    }

    /** Reverse complement of one k-mer packed 2 bits/base, first base in the high bits. */
    private static long reverseComplement(long kmer, int k) {
        long r = 0;
        for (int i = 0; i < k; i++) {
            r = (r << 2) | (3 - (kmer & 3));
            kmer >>>= 2;
        }
        return r;
    }

    /**
     * The window's reverse-complement k-mer set, DERIVED from its forward set rather
     * than stored alongside it.
     *
     * <p>The two sets are not independent: the rolling loop that built the forward
     * set emitted, at every position, a k-mer and its reverse complement, and
     * reverse-complementation is a bijection on packed k-mers. So it maps duplicates
     * to duplicates and distinct to distinct, and re-sorting the mapped forward set
     * reproduces the reverse set exactly — this was checked against the previous
     * two-set builder over 400 blocks (complex, tandem, N-gapped and homopolymer, k
     * across [11,23]) with zero mismatches and zero size differences.
     *
     * <p>Keeping both cost 16 bytes per distinct k-mer per block, and the class doc
     * warned a large combined run wants -Xmx for it. Deriving the seed's set on
     * demand halves that to 8, and is not even a time trade: the old builder ran a
     * sortDedup per block for the reverse set, whereas a block pays this sort only
     * on the iteration where it is the seed, and each block is a seed at most once.
     */
    private static long[] reverseComplementSet(long[] fwd, int k) {
        if (fwd.length == 0) {
            return EMPTY_L;
        }
        final long[] r = new long[fwd.length];
        for (int i = 0; i < fwd.length; i++) {
            r[i] = reverseComplement(fwd[i], k);
        }
        Arrays.sort(r);   // already distinct: rc is injective, so no dedup pass is needed
        return r;
    }

    /** Builds cFwd for every block from the store, in parallel. */
    private void buildContainmentSets(int numSeqs) {
        final long[][] tF = new long[numSeqs][];
        final long[][] seqData = sequences;
        final SeqStore st = store;
        final int k = containK;
        IntStream.range(0, numSeqs).parallel().forEach(j -> {
            long offset = seqData[j][0];
            int length = (int) seqData[j][1];   // a repeat block always fits in an int
            tF[j] = forwardSetStore(st, offset, length, k);
        });
        cFwd = tF;
    }

    /**
     * Greedy seed-and-absorb clustering by CONTAINMENT. Same skeleton as
     * {@link #cluster} — sequential, deterministic outer loop; seed is the longest
     * remaining block ({@code active[0]}); compact ids from 3; a rejected singleton
     * reuses its id; the returned value is the exclusive upper bound. Only the
     * per-candidate test differs: a candidate joins the seed's cluster when a large
     * enough share of its k-mers occurs in the seed (see {@link #scanSeedContainment}).
     */
    private int clusterContainment(int numSeqs) {
        this.containK = containmentK((int) sequences[numSeqs / 2][1]);   // k from the median block length
        buildContainmentSets(numSeqs);

        final int[] active = new int[numSeqs];
        for (int j = 0; j < numSeqs; j++) {
            active[j] = j;
        }
        int activeLen = numSeqs;

        final int[] assignBuf = new int[numSeqs];
        int nextCluster = 2;

        while (activeLen > 0) {
            int seed = active[0];

            nextCluster++;
            clusterIds[seed] = nextCluster;

            int matched = scanSeedContainment(seed, active, activeLen, nextCluster, assignBuf);

            if (matched == 0) {
                nextCluster--;
                clusterIds[seed] = 0;
            }

            int w = 0;
            for (int t = 1; t < activeLen; t++) {
                int idx = active[t];
                if (clusterIds[idx] == 0) {
                    active[w++] = idx;
                }
            }
            activeLen = w;
        }
        return nextCluster + 1;
    }

    /**
     * Containment scan of one seed against its candidates. For each candidate the
     * share of its forward k-mers found in the seed is measured in both
     * orientations — {@code sff} against the seed's forward set, {@code sfr}
     * against the seed's reverse-complement set — and the larger decides both
     * membership and strand. A candidate is absorbed when it clears both
     * {@link #containmentFloor} and MIN_CONTAINMENT; it is assigned {@code +clusterId} for a forward
     * match and {@code -clusterId} for a reverse one, which is how the strand
     * reaches the report (the caller turns a negative id into a negative length).
     *
     * <p>Deterministic: the scan writes only private assignBuf slots and reads
     * read-only sets, so parallel and single-threaded give the same clusters; the
     * id writes are applied afterwards, in index order.
     */
    private int scanSeedContainment(int seedIdx, int[] active, int activeLen, int clusterId, int[] assignBuf) {
        final int from = 1, to = activeLen;
        if (to <= from) {
            return 0;
        }

        final long[] seedF = cFwd[seedIdx];
        final long[] seedR = reverseComplementSet(seedF, containK);   // derived once per seed
        final int nSeed = seedF.length;
        if (nSeed < MIN_SHARED_KMERS) {
            return 0;   // seed too small: no candidate can reach MIN_SHARED_KMERS
        }

        final int posId = clusterId, negId = -clusterId;
        final int[] act = active;
        final long[][] F = cFwd;
        final double pow4k = Math.pow(4, containK);
        final int[] out = assignBuf;

        IntStream stream = IntStream.range(from, to);
        if ((to - from) >= PARALLEL_SCAN_THRESHOLD) {
            stream = stream.parallel();
        }

        stream.forEach(p -> {
            final int j = act[p];
            final long[] candF = F[j];
            final int nCand = candF.length;
            int v = 0;
            if (nCand >= MIN_SHARED_KMERS) {
                final int sff = intersectionCount(candF, seedF);   // same orientation
                final int sfr = intersectionCount(candF, seedR);   // reverse-complement orientation
                final int shared = Math.max(sff, sfr);
                if (shared >= MIN_SHARED_KMERS) {
                    final double c = (double) shared / nCand;
                    if (c > containmentFloor(nCand, nSeed, pow4k) && c >= MIN_CONTAINMENT) {
                        v = (sff >= sfr) ? posId : negId;
                    }
                }
            }
            out[p - from] = v;
        });

        int matched = 0;
        for (int p = from; p < to; p++) {
            int v = out[p - from];
            if (v != 0) {
                clusterIds[act[p]] = v;
                matched++;
            }
        }
        return matched;
    }

    /**
     * Containment version of {@link #matchClustersToReferences}: annotates each
     * cluster with the first reference a member is contained in. References become
     * forward/reverse k-mer sets (short ones, length <= MIN_REF_LENGTH, are
     * dropped to empty). Members are grouped by cluster exactly as the default
     * matcher (positive ids, i.e. seeds and forward members), and the returned
     * array keeps the same shape: {@code int[numClusters]}, slot clusterId holds
     * refIndex+1 on a hit else 0.
     */
    private int[] matchClustersToReferencesContainment(int numSeqs, int numClusters) {
        final int n = numReferenceSeqs;
        final long[][] rF = new long[n][];
        final long[][] rR = new long[n][];
        final String[] refs = referenceSeqs;
        final int k = containK;
        IntStream.range(0, n).parallel().forEach(r -> {
            int length = refs[r].length();
            if (length <= MIN_REF_LENGTH) {
                rF[r] = EMPTY_L;
                rR[r] = EMPTY_L;
                return;
            }
            rF[r] = forwardSet(refs[r], 0, length, k);
            rR[r] = reverseComplementSet(rF[r], k);
        });

        final int[] clusterToRef = new int[numClusters];

        int[] cnt = new int[numClusters];
        for (int i = 0; i < numSeqs; i++) {
            int c = clusterIds[i];
            if (c >= 3 && c < numClusters) {
                cnt[c]++;
            }
        }
        int[] memPtr = new int[numClusters + 1];
        for (int c = 0; c < numClusters; c++) {
            memPtr[c + 1] = memPtr[c] + cnt[c];
        }
        int[] members = new int[memPtr[numClusters]];
        int[] fill = new int[numClusters];
        for (int i = 0; i < numSeqs; i++) {
            int c = clusterIds[i];
            if (c >= 3 && c < numClusters) {
                members[memPtr[c] + fill[c]++] = i;
            }
        }

        final double pow4k = Math.pow(4, k);
        IntStream.range(3, numClusters).parallel().forEach(clusterId -> {
            int mStart = memPtr[clusterId], mEnd = memPtr[clusterId + 1];
            for (int mi = mStart; mi < mEnd; mi++) {
                int seq = members[mi];
                long[] memF = cFwd[seq];
                if (memF.length < MIN_SHARED_KMERS) {
                    continue;
                }
                int matchedRef = -1;
                for (int r = 0; r < n; r++) {
                    long[] refF = rF[r];
                    if (refF.length < MIN_SHARED_KMERS) {
                        continue;
                    }
                    int shared = Math.max(intersectionCount(memF, refF), intersectionCount(memF, rR[r]));
                    if (shared < MIN_SHARED_KMERS) {
                        continue;
                    }
                    double c = (double) shared / memF.length;
                    if (c > containmentFloor(memF.length, refF.length, pow4k) && c >= MIN_CONTAINMENT) {
                        matchedRef = r;
                        break;
                    }
                }
                if (matchedRef >= 0) {
                    clusterToRef[clusterId] = matchedRef + 1;   // 1-indexed
                    break;
                }
            }
        });

        return clusterToRef;
    }
}