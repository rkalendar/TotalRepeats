import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

/**
 * Clusters DNA sequences based on k-mer frequency vector similarity.
 * Sequences are represented as 4-mer frequency vectors and clustered by comparing
 * pairwise ratio profiles. Both forward and reverse-complement orientations are considered.
 * Optionally matches clusters against a set of reference sequences.
 */
public final class SequencesClustering2 { // original code a little bit slow compared to SequenceClustering

    // ── Tuning constants ────────────────────────────────────────────────────────
    private static final int KMER_SIZE = 4;
    private static final int QUICK_POS_THRESHOLD = 7;
    private static final int QUICK_NEG_THRESHOLD = 12;
    private static final int MIN_VECTOR_PAIRS = 7;    // min unique-kmer-vector pairs per comparison
    private static final int MIN_KMER_HITS = 3;       // min kmer occurrences to keep a vector slot
    private static final int MIN_REF_LENGTH = 90;     // min length of a reference sequence
    private static final int MIN_SIMILARITY = 70;     // percent similarity threshold (0–100)
    private static final double MIN_DIM_SIMILARITY = 0.7; // per-dimension similarity floor

    /** Canonical 4-mer palette (256 entries covering all non-N 4-mers). */
    private static final String[] KMERS = {
        "ttaa", "acgt", "agct", "tgca", "tcga", "catg", "ctag", "gatc",
        "gtac", "ccgg", "ggcc", "tttc", "gaaa", "attt", "aaat", "tttg",
        "caaa", "tgaa", "ttct", "aaca", "ttca", "agaa", "tgtt", "tctt",
        "ttga", "aaga", "tcaa", "ttgt", "aaac", "gttt", "acaa", "aaag",
        "tgat", "cttt", "atca", "tcat", "atga", "catt", "acat", "attg",
        "caat", "atgt", "aatg", "attc", "gaat", "gatt", "aatc", "tatt",
        "aata", "atct", "ttcc", "agtt", "aact", "gaag", "cttc", "agat",
        "aagt", "ggaa", "actt", "ccaa", "tcca", "ttgg", "tgga", "ttgc",
        "tctc", "gcaa", "gatg", "catc", "tcct", "caac", "gttg", "agga",
        "caag", "tgta", "aagg", "taca", "agca", "tatc", "cctt", "gata",
        "tgct", "acca", "atgg", "tggt", "atac", "gtat", "cttg", "ccat",
        "cata", "gtga", "tatg", "tctg", "gctt", "ttat", "gcat", "gttc",
        "ggat", "agag", "ataa", "atcc", "aagc", "ggtt", "gaac", "caga",
        "aacc", "atgc", "ttac", "gtaa", "ctga", "tcag", "ctgt", "atag",
        "ctat", "agta", "aggt", "tact", "ctcc", "atta", "taat", "tgtc",
        "tcac", "tgag", "acct", "acag", "ctca", "ggag", "gtca", "agtg",
        "taga", "tcta", "cctc", "taaa", "tgac", "ttta", "gaca", "gagg",
        "actc", "gagt", "actg", "tggc", "cagt", "taac", "gtta", "gtgg",
        "ggga", "tagt", "ccac", "gctg", "acta", "cacc", "cagc", "tccc",
        "ctgc", "gcag", "ggta", "gtct", "tacc", "ccag", "cact", "gtag",
        "ctgg", "ggtg", "agac", "tgcc", "tggg", "agtc", "cagg", "gact",
        "ggca", "cctg", "aggg", "taag", "ctaa", "gcca", "tagc", "ccct",
        "ctta", "gctc", "gcta", "ttag", "gagc", "ggct", "ctac", "ttcg",
        "gcac", "tagg", "cgaa", "ccta", "gtgc", "gggt", "tcgt", "aggc",
        "cgtt", "acga", "gcct", "aacg", "atcg", "cgat", "gacc", "ggtc",
        "ggac", "accc", "ccca", "gtcc", "ctcg", "agcc", "ccga", "cgag",
        "acgg", "tcgc", "cgct", "cgga", "tccg", "cgtc", "cacg", "gacg",
        "tcgg", "agcg", "gcga", "cgcc", "gccc", "accg", "cgta", "gccg",
        "tacg", "cgtg", "ccgc", "cgac", "ggcg", "gggc", "cggt", "cggc",
        "gtcg", "ccgt", "cgca", "tgcg", "gcgg", "cggg", "acgc", "gcgt",
        "cccg"
    };

    /** Pre-built immutable lookup: kmer string → vector index. */
    private static final Map<String, Integer> KMER_INDEX;
    private static final int NUM_KMERS;

    static {
        Map<String, Integer> map = new HashMap<>(KMERS.length * 2);
        for (int i = 0; i < KMERS.length; i++) {
            map.put(KMERS[i], i);
        }
        KMER_INDEX = Collections.unmodifiableMap(map);
        NUM_KMERS = KMER_INDEX.size();
    }

    // ── Instance state ──────────────────────────────────────────────────────────
    private int[][] sequences;           // [i] = {offset, length}, sorted by length desc
    private int[] clusterIds;            // cluster assignment per sequence
    private int[] refClusterIds;         // cluster → reference-sequence mapping
    private int clusterCount;
    private final String[] referenceSeqs;
    private final int numReferenceSeqs;

    // ── Constructor ─────────────────────────────────────────────────────────────

    /**
     * @param seq        concatenated sequence data
     * @param refSeqs    reference sequences for cluster annotation (nullable)
     * @param offsets    flat array of [offset0, length0, offset1, length1, …]
     * @param useParallelClustering if true, uses the parallel clustering strategy
     */
    public SequencesClustering2(String seq, String[] refSeqs, int[] offsets, boolean useParallelClustering) {
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

        this.clusterCount = useParallelClustering
                ? clusterParallel(seq, numSeqs)
                : clusterSequential(seq, numSeqs);

        this.refClusterIds = new int[clusterCount];
        if (numReferenceSeqs > 0) {
            this.refClusterIds = matchClustersToReferences(seq, numSeqs, clusterCount);
        }
    }

    // ── Public API ──────────────────────────────────────────────────────────────

    public int getClusterCount()         { return clusterCount; }
    public int[] getClusterIds()         { return clusterIds; }
    public int[] getReferenceIds()       { return refClusterIds; }
    public int[][] getSequenceOffsets()  { return sequences; }

    // ── Parsing ─────────────────────────────────────────────────────────────────

    private static int[][] parseAndSortByLengthDesc(int[] offsets, int numSeqs) {
        int[][] seqs = new int[numSeqs][2];
        for (int j = 0; j < numSeqs; j++) {
            seqs[j][0] = offsets[j * 2];
            seqs[j][1] = offsets[j * 2 + 1];
        }
        Arrays.sort(seqs, (a, b) -> Integer.compare(b[1], a[1]));
        return seqs;
    }

    // ── K-mer vector construction ───────────────────────────────────────────────

    /**
     * Builds forward and reverse-complement k-mer frequency vectors for a set of
     * subsequences extracted from {@code seq}. Vectors are normalised in-place:
     * slots with ≤ MIN_KMER_HITS are zeroed; others become length / (1 + count).
     *
     * @return {forwardVectors, reverseVectors} each [numSeqs][NUM_KMERS]
     */
    private int[][][] buildKmerVectors(String seq, int numSeqs, int[][] seqData) {
        int[][] forward = new int[numSeqs][NUM_KMERS];
        int[][] reverse = new int[numSeqs][NUM_KMERS];

        IntStream.range(0, numSeqs).parallel().forEach(j -> {
            if (clusterIds[j] != 0) return;

            int offset = seqData[j][0];
            int length = seqData[j][1];

            for (int i = 0; i <= length - KMER_SIZE; i++) {
                String kmer = seq.substring(offset + i, offset + i + KMER_SIZE);

                Integer idx = KMER_INDEX.get(kmer);
                if (idx != null) {
                    forward[j][idx]++;
                }

                Integer rcIdx = KMER_INDEX.get(dna.ComplementDNA2(kmer));
                if (rcIdx != null) {
                    reverse[j][rcIdx]++;
                }
            }

            normalizeVector(forward[j], length);
            normalizeVector(reverse[j], length);
        });

        return new int[][][] { forward, reverse };
    }

    /** Builds k-mer vectors specifically for reference sequences. */
    private int[][][] buildReferenceKmerVectors() {
        int[][] forward = new int[numReferenceSeqs][NUM_KMERS];
        int[][] reverse = new int[numReferenceSeqs][NUM_KMERS];

        IntStream.range(0, numReferenceSeqs).parallel().forEach(j -> {
            int length = referenceSeqs[j].length();
            if (length <= MIN_REF_LENGTH) return;

            for (int i = 0; i <= length - KMER_SIZE; i++) {
                String kmer = referenceSeqs[j].substring(i, i + KMER_SIZE);

                Integer idx = KMER_INDEX.get(kmer);
                if (idx != null) {
                    forward[j][idx]++;
                }

                Integer rcIdx = KMER_INDEX.get(dna.ComplementDNA2(kmer));
                if (rcIdx != null) {
                    reverse[j][rcIdx]++;
                }
            }

            normalizeVector(forward[j], length);
            normalizeVector(reverse[j], length);
        });

        return new int[][][] { forward, reverse };
    }

    private static void normalizeVector(int[] vector, int length) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = (vector[i] > MIN_KMER_HITS) ? (length / (1 + vector[i])) : 0;
        }
    }

    // ── Similarity computation ──────────────────────────────────────────────────

    /**
     * Computes similarity between two k-mer vectors using pairwise ratio profiles.
     *
     * @param vecA          first vector
     * @param vecB          second vector
     * @param useEarlyExit  enable quick-accept / quick-reject heuristics
     * @return similarity score in 0–100, or -1 if insufficient vector pairs
     */
    private int computeSimilarity(int[] vecA, int[] vecB, boolean useEarlyExit) {
        int pairCount = 0;
        int failCount = 0;
        double totalSim = 0.0;

        for (int k = 0; k < NUM_KMERS - 1; k++) {
            if (vecA[k] <= 0 || vecB[k] <= 0) continue;

            int matchCount = 0;
            int mismatchCount = 0;
            double matchWeight = 0.0;

            for (int h = k + 1; h < NUM_KMERS; h++) {
                if (vecA[h] <= 0 || vecB[h] <= 0) continue;

                matchCount++;

                // Compute ratio similarity: how similarly the two vectors
                // express the ratio between dimensions k and h.
                double ratioA = computeMinRatio(vecA[k], vecA[h]);
                double ratioB = computeMinRatio(vecB[k], vecB[h]);
                double dimSim = Math.min(ratioA, ratioB) / Math.max(ratioA, ratioB);

                if (dimSim > MIN_DIM_SIMILARITY && dimSim <= 1.0) {
                    matchWeight += dimSim;
                } else {
                    mismatchCount++;
                }

                if (useEarlyExit) {
                    if ((mismatchCount > matchWeight && mismatchCount > QUICK_NEG_THRESHOLD)
                            || matchWeight > QUICK_POS_THRESHOLD) {
                        break;
                    }
                }
            }

            if (matchCount > 0) {
                double avgSim = matchWeight / matchCount;
                totalSim += (avgSim > MIN_DIM_SIMILARITY) ? avgSim : 0.0;
                pairCount++;

                if (useEarlyExit && matchWeight + matchWeight < matchCount) {
                    failCount++;
                }
            }

            if (useEarlyExit && failCount > QUICK_NEG_THRESHOLD) {
                return -1; // early rejection
            }
        }

        if (pairCount <= MIN_VECTOR_PAIRS) {
            return -1;
        }
        return (int) (100.0 * totalSim) / pairCount;
    }

    /** Returns 100 * min(a,b) / max(a,b), i.e. a normalised ratio in [0, 100]. */
    private static double computeMinRatio(int a, int b) {
        return (a > b) ? (100.0 * b) / a : (100.0 * a) / b;
    }

    // ── Sequential clustering ───────────────────────────────────────────────────

    private int clusterSequential(String seq, int numSeqs) {
        int[][][] vectors = buildKmerVectors(seq, numSeqs, sequences);
        int[][] fwd = vectors[0];
        int[][] rev = vectors[1];

        int nextCluster = 2;

        for (int i = 0; i < numSeqs; i++) {
            if (clusterIds[i] != 0) continue;

            nextCluster++;
            clusterIds[i] = nextCluster;
            int clusterSize = 0;

            // Forward direction
            clusterSize += assignMatchesSequential(fwd, fwd, i, numSeqs, nextCluster, /* negate */ false);

            // Reverse-complement direction
            clusterSize += assignMatchesSequential(fwd, rev, i, numSeqs, nextCluster, /* negate */ true);

            // Remove singleton clusters
            if (clusterSize == 0) {
                nextCluster--;
                clusterIds[i] = 0;
            }
        }
        return nextCluster;
    }

    /**
     * Tries to assign unclustered sequences to cluster {@code clusterId}
     * by comparing them against the seed vector {@code vecSeed[seedIdx]}.
     *
     * @param vecSeed   vectors for the seed (always forward)
     * @param vecTarget vectors for candidates (forward or reverse)
     * @param negate    if true, assigned ids are stored as negative (reverse match)
     * @return number of newly assigned sequences
     */
    private int assignMatchesSequential(int[][] vecSeed, int[][] vecTarget,
                                        int seedIdx, int numSeqs, int clusterId, boolean negate) {
        int assigned = 0;
        for (int j = seedIdx + 1; j < numSeqs; j++) {
            if (clusterIds[j] != 0) continue;

            int similarity = computeSimilarity(vecSeed[seedIdx], vecTarget[j], /* earlyExit */ true);
            if (similarity > MIN_SIMILARITY) {
                clusterIds[j] = negate ? -clusterId : clusterId;
                assigned++;
            }
        }
        return assigned;
    }

    // ── Parallel clustering ─────────────────────────────────────────────────────

    private int clusterParallel(String seq, int numSeqs) {
        int[][][] vectors = buildKmerVectors(seq, numSeqs, sequences);
        int[][] fwd = vectors[0];
        int[][] rev = vectors[1];

        AtomicInteger nextCluster = new AtomicInteger(3);

        IntStream.range(0, numSeqs).parallel().forEach(i -> {
            if (clusterIds[i] != 0) return;

            int localCluster = nextCluster.getAndIncrement();
            Set<Integer> members = Collections.newSetFromMap(new ConcurrentHashMap<>());

            clusterIds[i] = localCluster;
            members.add(i);

            // Forward direction
            assignMatchesParallel(fwd, fwd, i, numSeqs, localCluster, false, members);

            // Reverse-complement direction
            assignMatchesParallel(fwd, rev, i, numSeqs, localCluster, true, members);

            // Remove singleton clusters
            if (members.size() == 1) {
                clusterIds[i] = 0;
                // NOTE: this leaves a gap in cluster IDs, which is acceptable.
                // The original decrementAndGet() was racy; gaps are harmless.
            }
        });

        return nextCluster.get();
    }

    private void assignMatchesParallel(int[][] vecSeed, int[][] vecTarget,
                                       int seedIdx, int numSeqs, int clusterId,
                                       boolean negate, Set<Integer> members) {
        for (int j = seedIdx + 1; j < numSeqs; j++) {
            if (clusterIds[j] != 0) continue;

            int similarity = computeSimilarity(vecSeed[seedIdx], vecTarget[j], /* earlyExit */ false);
            if (similarity > MIN_SIMILARITY) {
                synchronized (clusterIds) {
                    if (clusterIds[j] == 0) {
                        clusterIds[j] = negate ? -clusterId : clusterId;
                        members.add(j);
                    }
                }
            }
        }
    }

    // ── Reference matching ──────────────────────────────────────────────────────

    private int[] matchClustersToReferences(String seq, int numSeqs, int numClusters) {
        int[][][] refVectors = buildReferenceKmerVectors();
        int[][] refFwd = refVectors[0];
        int[][] refRev = refVectors[1];

        int[] clusterToRef = new int[numClusters];

        IntStream.range(3, numClusters).parallel().forEach(clusterId -> {
            outer:
            for (int i = 0; i < numSeqs; i++) {
                if (clusterIds[i] != clusterId) continue;

                int[] seqVector = buildSingleSequenceVector(seq, i);

                // Try forward, then reverse
                for (int[][] refVec : new int[][][] { refFwd, refRev }) {
                    int matchedRef = findMatchingReference(refVec, seqVector);
                    if (matchedRef >= 0) {
                        clusterToRef[clusterId] = matchedRef + 1; // 1-indexed
                        break outer;
                    }
                }
            }
        });

        return clusterToRef;
    }

    /** Builds a normalised k-mer vector for a single sequence from the main data. */
    private int[] buildSingleSequenceVector(String seq, int seqIdx) {
        int[] vector = new int[NUM_KMERS];
        int offset = sequences[seqIdx][0];
        int length = sequences[seqIdx][1];

        for (int j = 0; j <= length - KMER_SIZE; j++) {
            String kmer = seq.substring(offset + j, offset + j + KMER_SIZE);
            Integer idx = KMER_INDEX.get(kmer);
            if (idx != null) {
                vector[idx]++;
            }
        }
        normalizeVector(vector, length);
        return vector;
    }

    /**
     * Finds the first reference whose vector is sufficiently similar to {@code seqVector}.
     *
     * @return reference index (0-based), or -1 if no match
     */
    private int findMatchingReference(int[][] refVectors, int[] seqVector) {
        for (int r = 0; r < numReferenceSeqs; r++) {
            int similarity = computeSimilarity(refVectors[r], seqVector, /* earlyExit */ false);
            if (similarity > MIN_SIMILARITY) {
                return r;
            }
        }
        return -1;
    }
}
