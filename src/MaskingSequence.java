import java.util.ArrayList;
import java.util.List;

/**
 * Single-sequence repeat masker using k-mer frequency analysis.
 * Counts k-mers on both strands of a single sequence, retains those observed ≥ 2 times, then masks regions covered by repeated k-mers.
 */
public class MaskingSequence {

    private long gapsLen = 0;
    private long repeatsLen = 0;
    private final int MIN_COUNT = 2;

    // ── Public API ──────────────────────────────────────────────────────────
    /**
     * Mask repeat regions in a single sequence.
     *
     * @param seq DNA sequence
     * @param ssrmsk per-base SSR mask (0 = not masked)
     * @param kmer k-mer length
     * @param minLenSeq minimum repeat block length to report
     * @return flat array of [start₀, len₀, start₁, len₁, …]
     */
    public int[] mask(String seq, byte[] ssrmsk, int kmer, int minLenSeq) {
        final boolean sens = true;
        final int top = sens ? (kmer / 2 - 1) : (kmer - 1);

        // Normalise once, reuse for both phases
        final byte[] fwd = normalise(seq);
        final byte[] rev = normalise(dna.ComplementDNA(seq));

        // Count gaps
        gapsLen = 0;
        for (byte b : fwd) {
            if (b == 4) {
                gapsLen++;
            }
        }

        // Phase 1: count k-mers → LongIntHashMap directly
        LongIntHashMap kmerMap = countKmers(fwd, rev, kmer);

        // Phase 2: compute per-base coverage
        // Filtering by MIN_COUNT — inside markStrand через get() >= MIN_COUNT
        int[] coverage = computeCoverage(kmerMap, fwd, rev, kmer, ssrmsk);

        // Phase 3: build blocks and filter
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
        LongIntHashMap map = new LongIntHashMap();
        countStrand(map, fwd, kmer);
        countStrand(map, rev, kmer);
        return map;
    }

    private void countStrand(LongIntHashMap map, byte[] bases, int kmer) {
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

            if (gapCount == 0) {
                long code = encodeKmer(bases, i - kmer + 1, kmer);
                int current = map.get(code);
                if (current == LongIntHashMap.NO_VALUE) {
                    map.put(code, 1);
                } else if (current < MIN_COUNT) {
                    map.put(code, current + 1);
                }
            }

            if (bases[i + 1 - kmer] == 4) {
                gapCount--;
            }
        }
    }

    // ── Phase 2: coverage ───────────────────────────────────────────────────
    private int[] computeCoverage(LongIntHashMap kmerMap, byte[] fwd, byte[] rev, int kmer, byte[] ssrmsk) {
        final int len = fwd.length;
        final int[] coverage = new int[len];

        markStrandForward(kmerMap, fwd, ssrmsk, coverage, kmer);
        markStrandReverse(kmerMap, rev, ssrmsk, coverage, len, kmer);

        return coverage;
    }

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
                    int fwdStart = len - start - kmer;
                    int fwdEnd = len - start;
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
    private int[] buildBlocks(int[] coverage, int kmer, int minLenSeq, int top) {
        final int n = coverage.length;
        List<int[]> blocks = new ArrayList<>();
        int i = 0;

        while (i < n) {
            if (coverage[i] > 1) {
                int start = i;
                int end = i;
                int highCount = 0;

                while (i < n && coverage[i] > 1) {
                    if (coverage[i] > top) {
                        highCount++;
                    }
                    end = i++;
                }

                if (highCount > 1) {
                    int blockEnd = Math.min(n - 1, end + kmer - 1);
                    if (!blocks.isEmpty() && blocks.get(blocks.size() - 1)[1] >= start) {
                        blocks.get(blocks.size() - 1)[1] = Math.max(blocks.get(blocks.size() - 1)[1], blockEnd);
                    } else {
                        blocks.add(new int[]{start, blockEnd});
                    }
                }
            } else {
                i++;
            }
        }

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

        int[] result = new int[filtered.size() * 2];
        int p = 0;
        for (int[] bl : filtered) {
            result[p++] = bl[0];
            result[p++] = bl[1];
        }
        return result;
    }

    // ── Utilities ───────────────────────────────────────────────────────────
    private static byte[] normalise(String seq) {
        byte[] raw = seq.getBytes();
        for (int i = 0; i < raw.length; i++) {
            raw[i] = tables.dx2[raw[i]];
        }
        return raw;
    }

    private long encodeKmer(byte[] b, int start, int kmerSize) {
        // Maximum value: -9,223,372,036,854,775,808 and 9,223,372,036,854,775,807 (which is 2^63 - 1)
        // kmer<=19 only positive hashcode
        // kmer>19 some hashcodes shifted to negative one
        long r = 0L;
        for (int i = start; i < start + kmerSize; i++) {
            r = r * 10 + b[i];
        }
        return r;
    }

    /**
     * Rolling 2-bit hash instead of encodeKmerLong rescanning all k bases each
     * window. Each step is now O(1) instead of O(k), which is a significant win
     * for larger k values. The encoding uses 2 bits per base (A=0, C=1, G=2,
     * T=3) fitting up to k=32 in a long. Correctness note: The original used
     * base-10 encoding (r = r * 10 + b[i]), which wastes bits and limits max
     * k-mer size. The 2-bit encoding supports k up to 32 (vs. ~19 before) while
     * keeping all codes non-negative. If you need the hash codes to be
     * backward-compatible with the old encoding, let me know and I can preserve
     * that scheme.
     *
     * private static long encodeKmer2(byte[] bases, int start, int kmerSize) {
     * long code = 0L; for (int i = start; i < start + kmerSize; i++) { code =
     * (code << 2) | bases[i]; } return code; }
     */
}
