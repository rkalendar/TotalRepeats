import java.util.ArrayList;
import java.util.List;

/**
 * Single-sequence repeat masker using k-mer frequency analysis.
 *
 * Algorithm: 1. Count k-mers on both strands of the input sequence. 2. Retain
 * k-mers observed ≥ MIN_COUNT times (= repeated). 3. Mark per-base coverage by
 * repeated k-mers. 4. Extract contiguous blocks with sufficient coverage and
 * length.
 */
public class MaskingSequence {

    private long gapsLen = 0;
    private long repeatsLen = 0;

    /**
     * A k-mer must appear at least this many times (across both strands) to be
     * considered a repeat.
     */
    private static final int MIN_COUNT = 2;

    /**
     * A position must be covered by more than this many repeated k-mers before
     * it is counted as a "high-coverage" base that anchors a repeat block. (set
     * to 1 so that coverage > 1 means ≥ 2 overlapping k-mers)
     */
    private static final int HIGH_COV_THRESHOLD = 1;

    // ── Public API ──────────────────────────────────────────────────────────
    /**
     * Mask repeat regions in a single sequence.
     *
     * @param seq DNA sequence (A/C/G/T/N, any case)
     * @param ssrmsk per-base SSR mask; positions where ssrmsk[i] != 0 are
     * excluded from coverage marking
     * @param kmer k-mer length
     * @param minLenSeq minimum repeat-block length to report
     * @return flat array of [start₀, len₀, start₁, len₁, …]
     */
    public int[] mask(String seq, byte[] ssrmsk, int kmer, int minLenSeq) {
        /*
         * sens=true selects half-kmer sensitivity: blocks whose internal
         * high-coverage span exceeds kmer/2-1 are kept.  This prevents
         * very-short spurious repeats at the boundary of real repeats.
         */
        final boolean sens = true;
        final int top = sens ? (kmer / 2 - 1) : (kmer - 1);

        // Normalise once; both strands reuse the same byte-coded arrays.
        final byte[] fwd = normalise(seq);
        final byte[] rev = normalise(dna.ComplementDNA(seq));

        // Count gap (N) bases on the forward strand.
        gapsLen = 0;
        for (byte b : fwd) {
            if (b == 4) {
                gapsLen++;
            }
        }

        // Phase 1 – count k-mers on both strands.
        LongIntHashMap kmerMap = countKmers(fwd, rev, kmer);

        // Phase 2 – per-base coverage by repeated k-mers.
        int[] coverage = computeCoverage(kmerMap, fwd, rev, ssrmsk, kmer);

        // Phase 3 – build repeat blocks and apply length filter.
        return buildBlocks(coverage, kmer, minLenSeq, top);
    }

    public long gapsLength() {
        return gapsLen;
    }

    public long repeatLength() {
        return repeatsLen;
    }

    // ── Phase 1: k-mer counting ─────────────────────────────────────────────
    private LongIntHashMap countKmers(byte[] fwd, byte[] rev, int kmer) {
        // Pre-size so the map does not resize while filling.  The distinct
        // k-mer count is bounded by the number of windows on both strands and
        // by 4^kmer; withExpectedKeys() caps the capacity at MAX_CAPACITY.
        long expected = estimateDistinctKmers(fwd.length, kmer);
        LongIntHashMap map = LongIntHashMap.withExpectedKeys(expected);
        countStrand(map, fwd, kmer);
        countStrand(map, rev, kmer);
        return map;
    }

    /**
     * Upper bound on the number of distinct k-mers that {@link #countKmers}
     * will store: the number of k-mer windows across both strands, but never
     * more than the 4^kmer possible gap-free k-mers (alphabet A/C/G/T only).
     *
     * This is an upper bound — for highly repetitive input the true distinct
     * count is lower, so the map may be over-sized. If memory is tight, drop
     * the factor of 2 (use a single strand's window count) as a tighter, still
     * usually-sufficient estimate. The choice affects only speed/memory, never
     * correctness.
     */
    private static long estimateDistinctKmers(int seqLen, int kmer) {
        long windowsPerStrand = Math.max(0L, (long) seqLen - kmer + 1);
        long windows = 2L * windowsPerStrand;                  // forward + reverse
        long alphabetBound = (kmer < 32) ? (1L << (2 * kmer)) // 4^kmer
                : Long.MAX_VALUE;
        return Math.min(windows, alphabetBound);
    }

    /**
     * Counts all gap-free k-mers in {@code bases}, incrementing the map entry
     * up to MIN_COUNT. Capping at MIN_COUNT keeps the map small — we only need
     * to know whether the count reaches the threshold, not the exact value.
     *
     * Uses {@link LongIntHashMap#incrementCapped} so each k-mer resolves its
     * slot in a single probe sequence (insert-1-or-bump), instead of a separate
     * get followed by put.
     */
    private void countStrand(LongIntHashMap map, byte[] bases, int kmer) {
        final int len = bases.length;

        // Initialise the gap counter for the first (kmer-1) bases.
        int gapCount = 0;
        for (int i = 0; i < kmer - 1; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
        }

        for (int i = kmer - 1; i < len; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }

            if (gapCount == 0) {
                long code = encodeKmer(bases, i - kmer + 1, kmer);
                // Insert with count 1 if absent, else bump up to MIN_COUNT.
                // Once at MIN_COUNT the entry is a "confirmed repeat" and the
                // cap stops further increments — single hash probe per k-mer.
                map.incrementCapped(code, MIN_COUNT);
            }

            // Slide the window: drop the base leaving the left end.
            if (bases[i - kmer + 1] == 4) {
                gapCount--;
            }
        }
    }

    // ── Phase 2: coverage ───────────────────────────────────────────────────
    private int[] computeCoverage(LongIntHashMap kmerMap,
            byte[] fwd, byte[] rev,
            byte[] ssrmsk, int kmer) {
        final int len = fwd.length;
        final int[] coverage = new int[len];

        markStrandForward(kmerMap, fwd, ssrmsk, coverage, kmer);
        markStrandReverse(kmerMap, rev, ssrmsk, coverage, len, kmer);

        return coverage;
    }

    /**
     * For each gap-free k-mer on the forward strand whose count ≥ MIN_COUNT,
     * increment coverage for every base it spans (unless the base is
     * SSR-masked).
     */
    private void markStrandForward(LongIntHashMap kmerMap, byte[] bases, byte[] ssrmsk, int[] coverage, int kmer) {
        final int len = bases.length;
        int gapCount = 0;
        for (int i = 0; i < kmer - 1; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
        }

        for (int i = kmer - 1; i < len; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }

            int start = i - kmer + 1;
            if (gapCount == 0) {
                long code = encodeKmer(bases, start, kmer);
                if (kmerMap.get(code) >= MIN_COUNT) {
                    for (int j = start; j <= i; j++) {
                        if (ssrmsk[j] == 0) {
                            coverage[j]++;
                        }
                    }
                }
            }

            if (bases[start] == 4) {
                gapCount--;
            }
        }
    }

    /**
     * Same as {@link #markStrandForward} but for the reverse-complement strand.
     * A k-mer at position {@code [start, start+kmer-1]} on {@code rev} maps to
     * forward positions {@code [len-start-kmer, len-start-1]}.
     */
    private void markStrandReverse(LongIntHashMap kmerMap, byte[] bases, byte[] ssrmsk, int[] coverage, int len, int kmer) {
        int gapCount = 0;
        for (int i = 0; i < kmer - 1; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }
        }

        for (int i = kmer - 1; i < len; i++) {
            if (bases[i] == 4) {
                gapCount++;
            }

            int start = i - kmer + 1;
            if (gapCount == 0) {
                long code = encodeKmer(bases, start, kmer);
                if (kmerMap.get(code) >= MIN_COUNT) {
                    // Map reverse coordinates → forward coordinates.
                    int fwdStart = len - start - kmer;   // inclusive
                    int fwdEnd = len - start;           // exclusive
                    for (int j = fwdStart; j < fwdEnd; j++) {
                        if (ssrmsk[j] == 0) {
                            coverage[j]++;
                        }
                    }
                }
            }

            if (bases[start] == 4) {
                gapCount--;
            }
        }
    }

    // ── Phase 3: block building ─────────────────────────────────────────────
    /**
     * Converts the per-base coverage array into a list of repeat blocks.
     *
     * A run of positions is candidate if coverage[i] > HIGH_COV_THRESHOLD (i.e.
     * covered by at least 2 overlapping repeated k-mers). Within that run, the
     * block is retained only when at least 2 positions exceed {@code top} (the
     * half-kmer or full-kmer threshold), which guards against very short
     * spurious signals. Retained blocks are extended by (kmer-1) at the right
     * end and merged if they overlap.
     */
    private int[] buildBlocks(int[] coverage, int kmer, int minLenSeq, int top) {
        final int n = coverage.length;
        List<int[]> blocks = new ArrayList<>();
        int i = 0;

        while (i < n) {
            if (coverage[i] > HIGH_COV_THRESHOLD) {
                int start = i;
                int end = i;
                int highCount = 0;

                while (i < n && coverage[i] > HIGH_COV_THRESHOLD) {
                    if (coverage[i] > top) {
                        highCount++;
                    }
                    end = i++;
                }

                // Require at least 2 positions with "high" coverage to confirm.
                if (highCount > 1) {
                    int blockEnd = Math.min(n - 1, end + kmer - 1);
                    if (!blocks.isEmpty()
                            && blocks.get(blocks.size() - 1)[1] >= start) {
                        // Merge with the previous block.
                        blocks.get(blocks.size() - 1)[1]
                                = Math.max(blocks.get(blocks.size() - 1)[1], blockEnd);
                    } else {
                        blocks.add(new int[]{start, blockEnd});
                    }
                }
            } else {
                i++;
            }
        }

        // Filter by minimum length and accumulate total repeat length.
        long totalLen = 0L;
        List<int[]> filtered = new ArrayList<>();
        for (int[] bnd : blocks) {
            int length = bnd[1] - bnd[0] + 1;
            if (length > minLenSeq) {
                filtered.add(new int[]{bnd[0], length});
                totalLen += length;
            }
        }
        repeatsLen = totalLen;

        // Pack into a flat [start, len, start, len, ...] result array.
        int[] result = new int[filtered.size() * 2];
        int p = 0;
        for (int[] bl : filtered) {
            result[p++] = bl[0];
            result[p++] = bl[1];
        }
        return result;
    }

    // ── Encoding ────────────────────────────────────────────────────────────
    /**
     * Encodes a k-mer as a base-5 long integer.
     *
     * Alphabet mapping (from {@code tables.dx2}): A=0, C=1, G=2, T=3, N(gap)=4
     *
     * Using base 5 (one digit per base): - Clean: no ambiguity, no hash
     * collision from different k-mers. - Safe range: 5^27 = 7.45×10^18 <
     * Long.MAX_VALUE (9.22×10^18), so k-mers up to length 27 always produce
     * non-negative codes. - Codes of gap-containing k-mers are never inserted
     * (gapCount guard), so digit value 4 appears only in invalid/unreachable
     * codes.
     *
     * NOTE: Base-5 codes cannot equal EMPTY_KEY (Long.MIN_VALUE) or DELETED_KEY
     * (Long.MIN_VALUE+1) for k ≤ 27 because those sentinels are negative and
     * base-5 codes are non-negative in that range.
     *
     * Replacing the old base-10 encoding: - Base 10 wasted ~1.3 bits per base.
     * - Base 10 limited safe (non-negative) k-mer length to 19. - Base 5 raises
     * that limit to 27, sufficient for typical repeat masking.
     */
    private static long encodeKmer(byte[] b, int start, int kmerSize) {
        long r = 0L;
        for (int i = start; i < start + kmerSize; i++) {
            r = r * 5 + b[i];
        }
        return r;
    }

    /**
     * Converts a DNA string to a byte array using the dx2 lookup table.
     */
    private static byte[] normalise(String seq) {
        byte[] raw = seq.getBytes();
        for (int i = 0; i < raw.length; i++) {
            raw[i] = tables.dx2[raw[i]];
        }
        return raw;
    }
}
