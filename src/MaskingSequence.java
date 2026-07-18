import java.util.ArrayList;
import java.util.List;

/**
 * Single-sequence repeat masker using k-mer frequency analysis.
 *
 * Algorithm: 1. Count k-mers on both strands of the input sequence. 2. Retain
 * k-mers observed ≥ MIN_COUNT times (= repeated). 3. Mark per-base coverage by
 * repeated k-mers. 4. Extract contiguous blocks with sufficient coverage and
 * length.
 *
 * Uses LongCountSet: a value-less set that packs a "seen >= 2" flag into the
 * count is only ever capped at MIN_COUNT, so one byte per slot suffices and the
 * table is ~25% smaller and more cache-friendly than the int-valued map.
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
        if (kmer > LongCountSet.MAX_KMER) {
            throw new IllegalArgumentException(
                "LongCountSet requires kmer <= " + LongCountSet.MAX_KMER
              + " (it packs the count into the high bit of the base-5 key). "
              + "For longer k-mers use the byte/int map.");
        }
        /*
         * sens=true selects half-kmer sensitivity: blocks whose internal
         * high-coverage span exceeds kmer/2-1 are kept.  This prevents
         * very-short spurious repeats at the boundary of real repeats.
         */
        final boolean sens = true;
        final int top = sens ? (kmer / 2 - 1) : (kmer - 1);

        // Normalise once; both strands reuse the same byte-coded arrays.
        final byte[] fwd = normalise(seq);
        final byte[] rev = normalise(Dna.ComplementDNA(seq));

        // Count gap (N) bases on the forward strand.
        gapsLen = 0;
        for (byte b : fwd) {
            if (b == 4) {
                gapsLen++;
            }
        }

        // Phase 1 – count k-mers on both strands.
        LongCountSet kmerMap = countKmers(fwd, rev, kmer);

        // Drop unique k-mers before the read phase (smaller, cache-resident table).
        kmerMap.retainAtLeast(MIN_COUNT);

        // Phase 2 – per-base coverage by repeated k-mers.
        int[] coverage = computeCoverage(kmerMap, fwd, rev, ssrmsk, kmer);

        // Phase 3 – build repeat blocks and apply length filter.
        return buildBlocks(coverage, kmer, minLenSeq, top);
    }

    /**
     * Combined / multi-chunk masker.
     *
     * Counts k-mers across ALL chunks into ONE shared map, so a k-mer shared
     * between different chunks (i.e. different input sequences) reaches
     * MIN_COUNT and is masked in every chunk it occurs in — cross-sequence
     * repeats are detected, which {@link #mask} called per sequence cannot do.
     *
     * No k-mer window ever crosses a chunk boundary, so no repeat can span a
     * sequence junction. Blocks are returned in GLOBAL long coordinates. No chunk
     * is ever concatenated, so the total length may exceed the ~2.1 Gb
     * single-String / single-array limit.
     *
     * Provably equivalent to: concatenate every sequence into one, block the
     * junctions, and call {@link #mask} once — without building the merged
     * sequence.
     *
     * @param seqs      individual sequences (chunks); each &lt; 2.1 Gb
     * @param ssrmsk    per-chunk SSR mask (ssrmsk[c].length == seqs[c].length())
     * @param kmer      k-mer length
     * @param minLenSeq minimum repeat-block length to report
     * @return flat GLOBAL array [start0, len0, start1, len1, …] (long)
     */
    public long[] maskCombined(String[] seqs, byte[][] ssrmsk, int kmer, int minLenSeq) {
        if (kmer > LongCountSet.MAX_KMER) {
            throw new IllegalArgumentException(
                "LongCountSet requires kmer <= " + LongCountSet.MAX_KMER
              + " (it packs the count into the high bit of the base-5 key). "
              + "For longer k-mers use the byte/int map.");
        }
        final boolean sens = true;
        final int top = sens ? (kmer / 2 - 1) : (kmer - 1);

        // ── Phase 1: ONE shared k-mer map over every chunk, both strands. ──
        final long alphabetBound = (kmer < 32) ? (1L << (2 * kmer)) : Long.MAX_VALUE;
        long expected = 0L;
        for (String s : seqs) {
            long w = Math.max(0L, (long) s.length() - kmer + 1);
            expected += Math.min(2L * w, alphabetBound);
        }
        expected = Math.min(expected, alphabetBound);

        LongCountSet map = LongCountSet.withExpectedKeys(expected);
        for (String s : seqs) {
            countStrand(map, normalise(s), kmer);
            countStrand(map, normalise(Dna.ComplementDNA(s)), kmer);
        }

        // Counting is finished — drop unique k-mers once. The table shrinks
        // many-fold, so the per-base marking below hits a small, cache-resident
        // table and most of the memory is freed. Lossless: a count-1 entry here
        // is genuinely unique. (Must come AFTER all chunks are counted, because a
        // k-mer in the last chunk can confirm a repeat in the first one.)
        map.retainAtLeast(MIN_COUNT);

        // ── Phase 2: coverage + blocks per chunk against the SHARED map. ──
        ArrayList<long[]> parts = new ArrayList<>(seqs.length);
        long totalRep = 0L, totalGap = 0L, off = 0L;

        for (int c = 0; c < seqs.length; c++) {
            final byte[] fwd = normalise(seqs[c]);
            final byte[] rev = normalise(Dna.ComplementDNA(seqs[c]));
            final int len = fwd.length;

            for (byte b : fwd) {
                if (b == 4) {
                    totalGap++;
                }
            }

            int[] coverage = new int[len];
            markStrandForward(map, fwd, ssrmsk[c], coverage, kmer);
            markStrandReverse(map, rev, ssrmsk[c], coverage, len, kmer);

            int[] local = buildBlocks(coverage, kmer, minLenSeq, top); // sets repeatsLen for THIS chunk
            totalRep += repeatsLen;

            long[] g = new long[local.length];
            for (int j = 0; j + 1 < local.length; j += 2) {
                g[j] = off + local[j];   // local start → global start
                g[j + 1] = local[j + 1]; // length (kept as-is)
            }
            parts.add(g);
            off += len;
        }

        repeatsLen = totalRep;
        gapsLen = totalGap;

        // Flatten the per-chunk global arrays into one result (single allocation).
        int total = 0;
        for (long[] g : parts) {
            total += g.length;
        }
        long[] result = new long[total];
        int p = 0;
        for (long[] g : parts) {
            System.arraycopy(g, 0, result, p, g.length);
            p += g.length;
        }
        return result;
    }

    public long gapsLength() {
        return gapsLen;
    }

    public long repeatLength() {
        return repeatsLen;
    }

    // ── Phase 1: k-mer counting ─────────────────────────────────────────────
    private LongCountSet countKmers(byte[] fwd, byte[] rev, int kmer) {
        // Pre-size so the map does not resize while filling.  The distinct
        // k-mer count is bounded by the number of windows on both strands and
        // by 4^kmer; withExpectedKeys() caps the capacity at MAX_CAPACITY.
        long expected = estimateDistinctKmers(fwd.length, kmer);
        LongCountSet map = LongCountSet.withExpectedKeys(expected);
        countStrand(map, fwd, kmer);
        countStrand(map, rev, kmer);
        return map;
    }

    /**
     * Upper bound on the number of distinct k-mers that {@link #countKmers}
     * will store: the number of k-mer windows across both strands, but never
     * more than the 4^kmer possible gap-free k-mers (alphabet A/C/G/T only).
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
     */
    private void countStrand(LongCountSet map, byte[] bases, int kmer) {
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
                map.incrementCapped(code, MIN_COUNT);
            }

            // Slide the window: drop the base leaving the left end.
            if (bases[i - kmer + 1] == 4) {
                gapCount--;
            }
        }
    }

    // ── Phase 2: coverage ───────────────────────────────────────────────────
    private int[] computeCoverage(LongCountSet kmerMap,
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
    private void markStrandForward(LongCountSet kmerMap, byte[] bases, byte[] ssrmsk, int[] coverage, int kmer) {
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
    private void markStrandReverse(LongCountSet kmerMap, byte[] bases, byte[] ssrmsk, int[] coverage, int len, int kmer) {
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
     * Alphabet mapping (from {@code Tables.dx2}): A=0, C=1, G=2, T=3, N(gap)=4
     * Base 5 keeps codes non-negative and collision-free for k ≤ 27.
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
            raw[i] = Tables.dx2[raw[i]];
        }
        return raw;
    }
}
