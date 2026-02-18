import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Identifies and masks repeated/low-complexity regions in a DNA sequence using
 * k-mer counting on both the forward strand and its reverse complement.
 *
 * <p>Usage:
 * <pre>
 *   MaskingPairwiseAlignmentSequence masker = new MaskingPairwiseAlignmentSequence();
 *   int[] regions = masker.mask(sequence, kmerSize, minBlockLen);
 *   // regions is a flat [start, length, start, length, ...] array of masked blocks
 * </pre>
 */
public class MaskingPairwiseAlignmentSequence {

    private byte[] mask;
    private long noGapsLen = 0;
    private long repeatsLen = 0;

    // -------------------------------------------------------------------------
    // Public API
    // -------------------------------------------------------------------------

    /**
     * Identifies repeated blocks in {@code seq} and returns them as a flat
     * {@code [start0, len0, start1, len1, ...]} array.
     *
     * @param seq       the DNA sequence (A/C/G/T/N)
     * @param kmer      k-mer length used for repeat detection
     * @param minLenSeq minimum block length to include in the result
     * @return flat array of (start, length) pairs for each masked block
     */
    public int[] mask(String seq, int kmer, int minLenSeq) {
        mask = maskingL(seq, kmer, minLenSeq);
        int n = mask.length;

        List<int[]> blocks = new ArrayList<>();

        int i = 0;
        while (i < n) {
            if (mask[i] == 1) {
                int start = i;

                // Advance i to the end of this run of 1s.
                while (i < n && mask[i] == 1) {
                    i++;
                }
                int end = i - 1;                   // inclusive end of the run
                int runLen = end - start + 1;       // replaces redundant highCount

                if (runLen > kmer) {
                    int seqStart = start;
                    int seqEnd   = Math.min(n - 1, end + kmer - 1);

                    if (blocks.isEmpty() || blocks.get(blocks.size() - 1)[1] < seqStart) {
                        blocks.add(new int[]{seqStart, seqEnd});
                    } else {
                        // Extend the last block instead of adding a new one.
                        int[] prev = blocks.get(blocks.size() - 1);
                        prev[1] = Math.max(prev[1], seqEnd);
                    }
                }
            } else {
                i++;
            }
        }

        int totalLen = 0;
        List<int[]> filtered = new ArrayList<>();
        for (int[] block : blocks) {
            int start  = block[0];
            int end    = block[1];
            int length = end - start + 1;
            if (length > minLenSeq) {
                filtered.add(new int[]{start, length});
                totalLen += length;
            }
        }
        this.repeatsLen = totalLen;

        int[] result = new int[filtered.size() * 2];
        int idx = 0;
        for (int[] block : filtered) {
            result[idx++] = block[0];
            result[idx++] = block[1];
        }
        return result;
    }

    /** Returns the number of gap-free positions seen during the last {@link #mask} call.
     * @return  */
    public long noGapsLength() {
        return noGapsLen;
    }

    /** Returns the raw byte mask produced by the last {@link #mask} call.
     * @return  */
    public byte[] getByteMask() {
        return mask;
    }

    /** Returns the total length of detected repeat regions from the last {@link #mask} call.
     * @return  */
    public long repeatLength() {
        return repeatsLen;
    }

    // -------------------------------------------------------------------------
    // Core masking (fast 2-bit k-mer path)
    // -------------------------------------------------------------------------

    /**
     * Builds a per-position mask array where {@code 1} marks positions that
     * belong to a repeated region.
     *
     * <p>Implementation notes:
     * <ul>
     *   <li>Uses 2-bit encoding (A=0, C=1, G=2, T=3, N/gap=4) — no String substrings.</li>
     *   <li>{@code HashMap<Long,int[]>}: {@code entry[0]} = count (negated while being marked);
     *       {@code entry[1]} = rolling index into the position array {@code u[]}.</li>
     *   <li>Layout: {@code b[0..l-1]} forward, {@code b[l]=4} sentinel,
     *       {@code b[l+1..l+l]} reverse complement.</li>
     * </ul>
     */
    private byte[] maskingL(String seq, int kmer, int minLenBlock) {
        if (seq == null || kmer <= 0 || minLenBlock <= 0) {
            return new byte[0];
        }
        final int l = seq.length();
        if (l == 0 || l < kmer) {
            return new byte[l];
        }

        final int g = l + l + 1;
        final byte[] msk = new byte[g];

        // Build the combined forward + sentinel + reverse-complement code array.
        byte[] b = seq.getBytes();          // ASCII A/C/G/T/N
        for (int i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];       // 0..4
        }
        b = Arrays.copyOf(b, g);
        b[l] = 4;                           // sentinel (N/gap)
        for (int i = 1; i <= l; i++) {
            b[l + i] = tables.cdnat2[b[l - i]];
        }

        final long keyMask = (kmer == 32) ? -1L : ((1L << (2 * kmer)) - 1L);

        // ---- Pass 1: count all valid k-mers on both strands ----
        HashMap<Long, int[]> km = new HashMap<>(Math.min(l, 1 << 20));
        noGapsLen = kmer - 1;

        countKmers(b, 0, l, kmer, keyMask, km, true /* track noGapsLen */);
        countKmers(b, l + 1, l, kmer, keyMask, km, false);

        // ---- Pass 2: mark repeated k-mers and tally buckets ----
        int distinct = 0, total = 0;

        distinct += markRepeated(b, 0,     l, kmer, keyMask, km);
        distinct += markRepeated(b, l + 1, l, kmer, keyMask, km);

        // Recount total after marking (signs were flipped; restore temporarily).
        for (int[] v : km.values()) {
            if (v[0] < 0) {
                total += -v[0];
            }
        }

        if (total == 0 || distinct == 0) {
            return Arrays.copyOf(msk, l);
        }

        // ---- Pass 3: collect occurrence positions into bucket arrays ----
        int[] u   = new int[total + 1];          // u[0] = number of distinct buckets used
        int[][] x1 = new int[distinct + 1][2];   // x1[i] = {count, startIndexInU}
        int z = 1, t = 0;

        int[] tz;
        tz = collectPositions(b, 0,     l, kmer, keyMask, km, u, x1, z, t);
        t = tz[0]; z = tz[1];
        tz = collectPositions(b, l + 1, l, kmer, keyMask, km, u, x1, z, t);
        t = tz[0];

        if (t == 0) {
            return Arrays.copyOf(msk, l);
        }

        // ---- Pass 4: pairwise extension within each k-mer bucket ----
        for (int bi = 1; bi <= t; bi++) {
            int count = x1[bi][0];
            int start = x1[bi][1];
            for (int a = 1; a <= count; a++) {
                int x = u[start + a - 1];
                for (int bpos = a + 1; bpos <= count; bpos++) {
                    int y = u[start + bpos - 1];
                    if (msk[x] != 0 || msk[y] != 0) {
                        continue;
                    }
                    int h = extendBlock(b, l, g, x, y, kmer);
                    if (h > minLenBlock) {
                        for (int r = 0; r < h; r++) {
                            msk[x + r] = 1;
                            if (y < l) {
                                msk[y + r] = 1;
                            } else {
                                msk[l - y + l - r] = 1;  // map RC position back to forward
                            }
                        }
                    }
                }
            }
        }

        // ---- Pass 5: fill small internal gaps within masked blocks ----
        fillGaps(msk, l, minLenBlock);

        return Arrays.copyOf(msk, l);
    }

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    /**
     * Scans a segment of {@code b} and increments the k-mer count for each
     * valid (gap-free) k-mer window found.
     *
     * @param b           the encoded sequence array
     * @param offset      start index in {@code b} for this strand
     * @param len         number of positions to scan (= sequence length {@code l})
     * @param kmer        k-mer length
     * @param keyMask     bitmask for the rolling key
     * @param km          the k-mer count map to update
     * @param trackGaps   if {@code true}, increment {@link #noGapsLen} for each valid window
     */
    private void countKmers(byte[] b, int offset, int len, int kmer, long keyMask,
                             HashMap<Long, int[]> km, boolean trackGaps) {
        int valid = 0;
        long key  = 0;
        for (int i = 0; i < len; i++) {
            int c = b[offset + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    if (trackGaps) noGapsLen++;
                    int[] v = km.get(key);
                    if (v == null) {
                        km.put(key, new int[]{1, 0});
                    } else {
                        v[0]++;
                    }
                }
            } else {
                valid = 0;
                key   = 0;
            }
        }
    }

    /**
     * Scans a strand and negates counts of k-mers that appear more than once,
     * so they can be identified as repeated in later passes.
     *
     * @return the number of distinct repeated k-mers found on this strand
     */
    private int markRepeated(byte[] b, int offset, int len, int kmer, long keyMask,
                              HashMap<Long, int[]> km) {
        int valid    = 0;
        long key     = 0;
        int distinct = 0;
        for (int i = 0; i < len; i++) {
            int c = b[offset + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    int[] v = km.get(key);
                    if (v != null && v[0] > 1) {
                        distinct++;
                        v[0] = -v[0];   // negate to mark as "seen"
                    }
                }
            } else {
                valid = 0;
                key   = 0;
            }
        }
        return distinct;
    }

    /**
     * Collects the starting positions of each occurrence of every repeated k-mer
     * into the shared {@code u[]} position array, and registers bucket metadata
     * in {@code x1[]}.
     *
     * <p>The window-start position stored in {@code u[]} is {@code offset + i + 1 - kmer},
     * which correctly gives the first index of the k-mer regardless of strand.
     *
     * @return {@code int[]{t, z}} — updated bucket count and next free slot in {@code u[]}
     */
    private int[] collectPositions(byte[] b, int offset, int len, int kmer, long keyMask,
                                   HashMap<Long, int[]> km, int[] u, int[][] x1,
                                   int z, int t) {
        int valid = 0;
        long key  = 0;
        for (int i = 0; i < len; i++) {
            int c = b[offset + i];
            if (c < 4) {
                key = ((key << 2) | c) & keyMask;
                if (++valid >= kmer) {
                    // Window start position in b[]: one past the last base, minus kmer length.
                    int pos = offset + i + 1 - kmer;
                    int[] v = km.get(key);
                    if (v != null) {
                        if (v[0] < 0) {             // first occurrence seen → open bucket
                            v[0] = -v[0];
                            if (v[0] > 1) {
                                v[1] = z;
                                u[z] = pos;
                                t++;
                                u[0]     = t;
                                x1[t][0] = v[0];
                                x1[t][1] = z;
                                z       += v[0];
                            }
                        } else if (v[0] > 1) {      // subsequent occurrence
                            v[1]++;
                            u[v[1]] = pos;
                        }
                    }
                }
            } else {
                valid = 0;
                key   = 0;
            }
        }
        return new int[]{t, z};
    }

    /**
     * Extends a matched pair of k-mer seeds ({@code x} and {@code y}) as far
     * as possible, allowing up to 3 mismatches before stopping.
     *
     * @return the extended block length
     */
    private static int extendBlock(byte[] b, int l, int g, int x, int y, int kmer) {
        int h = kmer, e = 0, p = 2;
        for (;;) {
            // Bounds check: don't run off the end of either strand.
            if ((y > l  && y + h > g - 1) || (y <= l && y + h > l - 1)) break;
            if ((x < l  && x + h > l - 1) || (x >= l && x + h > g - 1)) break;
            // Stop at a sentinel on both sides simultaneously.
            if (b[x + h] == 4 && b[y + h] == 4) break;

            if (b[x + h] == b[y + h]) {
                if (p > 0) e = 0;   // reset mismatch counter on a match run
                p++;
            } else {
                if (e > 2) {
                    h -= 3;         // retract the last 3 positions
                    break;
                }
                e++;
                p = 0;
            }
            h++;
        }
        return h;
    }

    /**
     * Fills gaps of ≤ {@code minLenBlock} zeros that sit between two masked
     * (value = 1) regions, merging them into a single continuous block.
     */
    private static void fillGaps(byte[] msk, int l, int minLenBlock) {
        int gapLen  = 0;
        boolean inBlock = false;

        for (int i = 0; i < l; i++) {
            if (msk[i] > 0) {
                if (!inBlock) {
                    inBlock = true;
                    gapLen  = 0;
                } else if (gapLen > 0) {
                    if (gapLen <= minLenBlock) {
                        // Fill the small gap that preceded this position.
                        Arrays.fill(msk, i - gapLen, i, (byte) 1);
                    }
                    gapLen = 0;
                }
            } else {
                if (inBlock) gapLen++;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Utility methods (kept private — not part of the public contract)
    // -------------------------------------------------------------------------

    private byte[] arrayTrim(byte[] srcArray, int n) {
        return Arrays.copyOf(srcArray, n);
    }

    private byte[] arrayExtendByte(byte[] srcArray, int n) {
        return Arrays.copyOf(srcArray, srcArray.length + n);
    }
}
