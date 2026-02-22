import java.util.ArrayList;
import java.util.List;

/**
 * Single-threaded repeat masker using k-mer frequency analysis.
 */
public class MaskingSequences {

    private long[] repeatsLen;
    private long[] gapsLen;
    private final int MIN_COUNT = 2;

    // ── Public API ──────────────────────────────────────────────────────────
    /**
     * Mask repeat regions in the given sequences.
     *
     * @param seq array of DNA sequences
     * @param ssrmsk per-base SSR mask (0 = not masked)
     * @param kmer k-mer length
     * @param minLenSeq minimum repeat block length to report
     * @param sens sensitivity flag (adjusts coverage threshold)
     * @return list of int[] per sequence; each array is [start₀, len₀, start₁,
     * len₁, …]
     */
    public ArrayList<int[]> mask(String[] seq, byte[][] ssrmsk, int kmer, int minLenSeq, boolean sens) {
        final int ns = seq.length;
        final int top = sens ? (kmer / 2 - 1) : (kmer - 1);

        repeatsLen = new long[ns];
        gapsLen = new long[ns];

        // Pre-normalise all sequences (forward + reverse complement)
        byte[][] fwd = new byte[ns][];
        byte[][] rev = new byte[ns][];

        // ── Phase 1: k-mer counting ────────────────────────────────────────
        LongIntHashMap map = new LongIntHashMap();
        for (int k = 0; k < ns; k++) {
            fwd[k] = normalise(seq[k]);
            rev[k] = normalise(dna.ComplementDNA(seq[k]));
            countKmers(map, k, fwd[k], rev[k], kmer);
        }

        // ── Phase 2: masking — map передаётся напрямую, без промежуточного Set
        ArrayList<int[]> results = new ArrayList<>(ns);
        for (int k = 0; k < ns; k++) {
            int[] coverage = computeCoverage(map, fwd[k], rev[k], kmer, ssrmsk[k]);
            int[] result = buildBlocks(coverage, kmer, minLenSeq, top, k);
            results.add(result);
        }

        return results;
    }

    public long[] gapsLength() {
        return gapsLen;
    }

    public long[] repeatLength() {
        return repeatsLen;
    }

    // ── Phase 1: k-mer counting ─────────────────────────────────────────────
    private void countKmers(LongIntHashMap map, int seqIdx, byte[] fwd, byte[] rev, int kmer) {
        // Count gaps
        long gaps = 0L;
        for (byte b : fwd) {
            if (b == 4) {
                gaps++;
            }
        }
        gapsLen[seqIdx] = gaps;

        // Count k-mers on both strands
        countStrand(map, fwd, kmer);
        countStrand(map, rev, kmer);
    }

    /**
     * Slide a k-mer window over one strand, incrementing counts for k-mers that
     * contain no gaps (base value 4). Counts are capped at MIN_COUNT.
     */
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

    // ── Phase 2: coverage & masking ─────────────────────────────────────────
    /**
     * LongIntHashMap: get(code) >= MIN_COUNT 
     */
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
                    // Map reverse-complement window → forward-strand coordinates
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

    // ── Block building ──────────────────────────────────────────────────────
    private int[] buildBlocks(int[] coverage, int kmer, int minLenSeq, int top, int seqIdx) {
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
        repeatsLen[seqIdx] = totalLen;

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

    private static long encodeKmer(byte[] bases, int start, int kmerSize) {
        long code = 0L;
        for (int i = start; i < start + kmerSize; i++) {
            code = code * 10 + bases[i];
        }
        return code;
    }
}
