/**
 * Read-only "virtual concatenation" of many sequences as one logical sequence,
 * addressed by a {@code long} global position — without ever building a single
 * >2.1 Gb String.
 *
 * Why this exists: a Java String (like any array) can hold at most
 * ~Integer.MAX_VALUE chars (~2.147e9).  Concatenating many chromosomes with
 * {@code String.join} throws "Requested string length exceeds VM limit" once the
 * total crosses that bound.  Each individual chromosome still fits in a String
 * (real chromosomes are < 2.1 Gb); only the concatenation was illegal.
 *
 * Coordinates: global positions are {@code long}.  Each segment is a normal
 * String; {@code segmentStart(i)} gives its global start.  A repeat block never
 * crosses a chromosome boundary, so the hot path resolves the owning segment
 * once and then indexes that segment with plain int-local coordinates (no
 * per-base binary search).
 */
public final class SeqStore {

    private final String[] segments;   // individual sequences (each ≤ ~2.1 Gb)
    private final long[]   starts;     // starts[i] = global start of segment i; starts[n] = total
    private final long     totalLength;

    public SeqStore(String[] segments) {
        this.segments = segments;
        this.starts   = new long[segments.length + 1];
        long acc = 0;
        for (int i = 0; i < segments.length; i++) {
            starts[i] = acc;
            acc += segments[i].length();
        }
        starts[segments.length] = acc;
        this.totalLength = acc;
    }

    /** Total length across all segments (the value the old {@code s.length()} held). */
    public long length() { return totalLength; }

    public int    segmentCount()      { return segments.length; }
    public String segment(int i)      { return segments[i]; }
    public long   segmentStart(int i) { return starts[i]; }

    /** Index of the segment containing global position {@code pos} (0 ≤ pos < length). */
    public int segmentOf(long pos) {
        int lo = 0, hi = segments.length - 1;
        while (lo < hi) {
            int mid = (lo + hi + 1) >>> 1;
            if (starts[mid] <= pos) lo = mid; else hi = mid - 1;
        }
        return lo;
    }

    /** Single base at a global position.  Use the resolve-once pattern for loops. */
    public char charAt(long pos) {
        int s = segmentOf(pos);
        return segments[s].charAt((int) (pos - starts[s]));
    }

    /**
     * Substring over global {@code [start, end)}; {@code end - start} must fit in
     * an int.  Fast path when the range is within one segment (always true for
     * repeat blocks).  Falls back to stitching only across boundaries.
     */
    public String substring(long start, long end) {
        int s = segmentOf(start);
        if (end <= starts[s + 1]) {                       // common: single segment
            return segments[s].substring((int) (start - starts[s]),
                                         (int) (end   - starts[s]));
        }
        int len = (int) (end - start);                    // rare: cross-boundary
        StringBuilder sb = new StringBuilder(len);
        long p = start;
        while (p < end) {
            int si = segmentOf(p);
            long segEnd = Math.min(end, starts[si + 1]);
            sb.append(segments[si], (int) (p - starts[si]), (int) (segEnd - starts[si]));
            p = segEnd;
        }
        return sb.toString();
    }
}
