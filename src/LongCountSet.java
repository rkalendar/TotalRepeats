/**
 * Open-addressing presence SET with a single "seen ≥ 2" flag per key, for
 * memory- and cache-optimal capped k-mer counting. NO value array at all.
 *
 * Packing: a base-5 k-mer code for k ≤ 21 satisfies code &lt; 5^21 &lt; 2^49, so
 * bit 49 of the slot is free. Each live slot stores
 *
 *     slot = code | (seenTwice ? FLAG : 0)        FLAG = 1L &lt;&lt; 49
 *
 * so the count is encoded in one bit and the int[]/byte[] values array is gone.
 * At MAX_CAPACITY this is one long[] of 8 GiB total (vs 9 GiB byte-valued, 12 GiB
 * int-valued), and each probe touches a SINGLE cache line (key+flag in one word)
 * instead of two. The distinct-k-mer ceiling is unchanged: still 1&lt;&lt;30 slots
 * (~805 M keys at load 0.75) — packing saves memory/latency, not the slot cap.
 *
 * REQUIRES kmer ≤ 21. For longer k-mers (base-5 code can reach 2^62) use the
 * byte- or int-valued map, which keep the count in a separate array.
 *
 * Drop-in API match with LongByteHashMap (incrementCapped / get / retainAtLeast /
 * withExpectedKeys), so MaskingSequence only needs the type renamed. Counts
 * saturate at 2: get() returns -1 (absent), 1 (seen once) or 2 (seen ≥ twice).
 *
 * Sentinels (bit 63 set, never collide with a packed slot whose bit 63 is clear):
 *   EMPTY_SLOT   = Long.MIN_VALUE
 *   DELETED_SLOT = Long.MIN_VALUE + 1
 */
public class LongCountSet {

    public static final int  NO_VALUE = -1;
    /** Largest k-mer length whose base-5 code stays below the FLAG bit. */
    public static final int  MAX_KMER = 21;

    private static final long EMPTY_SLOT   = Long.MIN_VALUE;
    private static final long DELETED_SLOT = Long.MIN_VALUE + 1;
    private static final int  MAX_CAPACITY = 1 << 30;

    private static final long FLAG      = 1L << 49;   // "seen ≥ 2" bit (free for k ≤ 21)
    private static final long CODE_MASK = FLAG - 1;   // bits 0..48 hold the base-5 code

    /** 2^64 / golden ratio, odd — Fibonacci-hashing multiplier. */
    private static final long GOLDEN = 0x9E3779B97F4A7C15L;

    private long[] keys;            // packed slots: code | FLAG?
    private int    size;
    private int    capacity;
    private int    threshold;
    private int    shift;           // = 64 - log2(capacity)

    // ---------- Constructors ----------

    public LongCountSet() {
        this(1 << 20);
    }

    public LongCountSet(int initialCapacity) {
        capacity  = nextPowerOfTwo(initialCapacity);
        keys      = new long[capacity];
        threshold = (int) (capacity * 0.75);
        shift     = 64 - Integer.numberOfTrailingZeros(capacity);
        java.util.Arrays.fill(keys, EMPTY_SLOT);
    }

    // ---------- Factories ----------

    public static LongCountSet withExpectedKeys(long expectedKeys) {
        return new LongCountSet(capacityForExpectedKeys(expectedKeys));
    }

    public static LongCountSet forSequenceLength(long sequenceLength) {
        return withExpectedKeys(sequenceLength);
    }

    /** True if this set can hold base-5 codes for the given k-mer length. */
    public static boolean fits(int kmer) {
        return kmer <= MAX_KMER;
    }

    private static int capacityForExpectedKeys(long expectedKeys) {
        if (expectedKeys < 1) return 2;
        long needed = (long) Math.ceil(expectedKeys / 0.75);
        if (needed >= MAX_CAPACITY) return MAX_CAPACITY;
        int cap = Integer.highestOneBit((int) needed);
        if (cap < needed) cap <<= 1;
        return Math.max(2, Math.min(cap, MAX_CAPACITY));
    }

    // ---------- Public API ----------

    /**
     * Observes {@code code}: inserts it with count 1 if absent, otherwise sets
     * the "seen ≥ 2" flag (saturating count at 2). {@code cap} ≥ 2 enables the
     * flag; the set has no states above 2.
     *
     * Precondition: {@code 0 ≤ code < 2^49} (guaranteed for base-5 codes with
     * kmer ≤ MAX_KMER). Codes touching the FLAG bit would corrupt the count, so
     * the caller must restrict kmer ≤ 21 (cheap to assert once, not per call).
     *
     * @return the count after the operation (1 or 2)
     */
    public int incrementCapped(long code, int cap) {
        ensureCapacity();

        int i            = index(code);
        int firstDeleted = -1;

        while (keys[i] != EMPTY_SLOT) {
            long slot = keys[i];
            if (slot == DELETED_SLOT) {
                if (firstDeleted == -1) firstDeleted = i;
            } else if ((slot & CODE_MASK) == code) {
                // Present — mark "seen ≥ 2" (saturating) when cap allows.
                if (cap >= 2 && (slot & FLAG) == 0) {
                    keys[i] = slot | FLAG;
                }
                return (keys[i] & FLAG) != 0 ? 2 : 1;
            }
            i = next(i);
        }

        if (firstDeleted != -1) i = firstDeleted;
        keys[i] = code;   // flag 0 → count 1
        size++;
        return 1;
    }

    /**
     * Returns the count for the key: NO_VALUE (-1) if absent, 1 if seen once,
     * 2 if seen ≥ twice. A present slot never collides with NO_VALUE.
     */
    public int get(long code) {
        int i = index(code);
        while (keys[i] != EMPTY_SLOT) {
            long slot = keys[i];
            // DELETED check first: (DELETED_SLOT & CODE_MASK) could otherwise
            // spuriously equal a small code.
            if (slot != DELETED_SLOT && (slot & CODE_MASK) == code) {
                return (slot & FLAG) != 0 ? 2 : 1;
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    public boolean containsKey(long code) {
        return get(code) != NO_VALUE;
    }

    /** Removes the key. Returns the count it had, or NO_VALUE if absent. */
    public int remove(long code) {
        int i = index(code);
        while (keys[i] != EMPTY_SLOT) {
            long slot = keys[i];
            if (slot != DELETED_SLOT && (slot & CODE_MASK) == code) {
                keys[i] = DELETED_SLOT;
                size--;
                return (slot & FLAG) != 0 ? 2 : 1;
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    /**
     * Drops every key with count &lt; {@code min} and rehashes once. For
     * {@code min == 2} (the masking threshold) this keeps only the flagged
     * (repeated) k-mers and removes all unique ones — the table shrinks
     * many-fold so the read phase becomes cache-resident. {@code min ≤ 1} keeps
     * everything; {@code min > 2} is treated as 2 (the set saturates there).
     *
     * @return number of entries removed
     */
    public int retainAtLeast(int min) {
        if (min <= 1) return 0;            // every live entry has count ≥ 1
        int removed = pruneUnflagged();    // drop count-1 (unflagged) keys
        if (removed > 0) rehash();
        return removed;
    }

    public int size()       { return size; }
    public boolean isEmpty() { return size == 0; }

    public void clear() {
        java.util.Arrays.fill(keys, EMPTY_SLOT);
        size = 0;
    }

    /**
     * Merges this set into {@code target}, carrying the "seen ≥ 2" flag. Used to
     * merge per-thread sets after parallel counting.
     */
    public void mergeInto(LongCountSet target) {
        for (int i = 0; i < capacity; i++) {
            long slot = keys[i];
            if (slot == EMPTY_SLOT || slot == DELETED_SLOT) continue;
            long code = slot & CODE_MASK;
            target.incrementCapped(code, 2);                 // ensure present
            if ((slot & FLAG) != 0) target.incrementCapped(code, 2); // bump to ≥2
        }
    }

    // ---------- Private helpers ----------

    /** Grows, or (at MAX_CAPACITY) prunes unique k-mers only — never repeats. */
    private void ensureCapacity() {
        if (size < threshold) return;

        if (capacity < MAX_CAPACITY) {
            resize();
            return;
        }

        int removed = pruneUnflagged();   // only count-1 keys; repeats are kept
        if (removed > 0) rehash();
        if (size >= threshold) {
            throw new IllegalStateException(
                    "Cannot free space: distinct repeated k-mers exceed table capacity "
                    + "(MAX_CAPACITY=" + MAX_CAPACITY + "). size=" + size);
        }
    }

    /** Fibonacci hashing on the CODE (not the packed slot). */
    private int index(long code) {
        return (int) ((code * GOLDEN) >>> shift);
    }

    private int next(int i) {
        return (i + 1) & (capacity - 1);
    }

    private void resize() {
        final int    oldCap  = capacity;
        final long[] oldKeys = keys;

        final int newCapacity = oldCap << 1;
        final long[] newKeys  = new long[newCapacity];
        java.util.Arrays.fill(newKeys, EMPTY_SLOT);

        capacity  = newCapacity;
        threshold = (int) (newCapacity * 0.75);
        shift     = 64 - Integer.numberOfTrailingZeros(newCapacity);
        keys      = newKeys;

        for (int i = 0; i < oldCap; i++) {
            long slot = oldKeys[i];
            if (slot == EMPTY_SLOT || slot == DELETED_SLOT) continue;
            int j = index(slot & CODE_MASK);
            while (newKeys[j] != EMPTY_SLOT) {
                j = next(j);
            }
            newKeys[j] = slot;   // move packed slot as-is (flag preserved)
        }
    }

    private void rehash() {
        long[] newKeys = new long[capacity];
        java.util.Arrays.fill(newKeys, EMPTY_SLOT);

        for (int i = 0; i < capacity; i++) {
            long slot = keys[i];
            if (slot == EMPTY_SLOT || slot == DELETED_SLOT) continue;
            int j = index(slot & CODE_MASK);
            while (newKeys[j] != EMPTY_SLOT) {
                j = next(j);
            }
            newKeys[j] = slot;
        }
        keys = newKeys;
    }

    /** Tombstones every live slot without the "seen ≥ 2" flag (count == 1). */
    private int pruneUnflagged() {
        int removed = 0;
        for (int i = 0; i < capacity; i++) {
            long slot = keys[i];
            if (slot != EMPTY_SLOT && slot != DELETED_SLOT && (slot & FLAG) == 0) {
                keys[i] = DELETED_SLOT;
                size--;
                removed++;
            }
        }
        return removed;
    }

    private static int nextPowerOfTwo(int n) {
        if (n <= 1) return 2;
        if (n > MAX_CAPACITY) {
            throw new IllegalArgumentException(
                    "initialCapacity exceeds maximum: " + MAX_CAPACITY);
        }
        int p = Integer.highestOneBit(n);
        return (p == n) ? n : p << 1;
    }

    @Override
    public String toString() {
        StringBuilder sb    = new StringBuilder("{");
        boolean       first = true;
        for (int i = 0; i < capacity; i++) {
            long slot = keys[i];
            if (slot != EMPTY_SLOT && slot != DELETED_SLOT) {
                if (!first) sb.append(", ");
                sb.append(slot & CODE_MASK).append('=').append((slot & FLAG) != 0 ? 2 : 1);
                first = false;
            }
        }
        return sb.append('}').toString();
    }
}
