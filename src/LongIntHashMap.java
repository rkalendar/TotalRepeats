/**
 * Hash table with open addressing.
 * Key: long, Value: int (positive)
 * No boxing, no objects — 5-7× faster than standard HashMap.
 *
 * Sentinel values:
 *   EMPTY_KEY   = Long.MIN_VALUE       (never a real key)
 *   DELETED_KEY = Long.MIN_VALUE + 1   (tombstone after remove)
 *
 * Changes vs the original (both behaviour-preserving for correctness):
 *   • Multiplicative (Fibonacci) hashing in index(): better dispersion of the
 *     structured base-5 k-mer codes → shorter probe chains → fewer cache misses.
 *     Only the bucket distribution changes, never get/put semantics.
 *   • retainAtLeast(min): public one-shot compaction to drop low-count keys
 *     (e.g. unique k-mers) after counting and before the read-heavy phase.
 */
public class LongIntHashMap {

    public static final int NO_VALUE = -1;
    private static final long EMPTY_KEY   = Long.MIN_VALUE;
    private static final long DELETED_KEY = Long.MIN_VALUE + 1;
    private static final int  MAX_CAPACITY = 1 << 30; // 1_073_741_824

    /** 2^64 / golden ratio, odd — the standard Fibonacci-hashing multiplier. */
    private static final long GOLDEN = 0x9E3779B97F4A7C15L;

    private int pruneThreshold = 1; // start by deleting entries with value == 1

    private long[] keys;
    private int[]  values;
    private int    size;
    private int    capacity;
    private int    threshold;
    private int    shift;           // = 64 - log2(capacity); top `log2(capacity)` bits index the table

    // ---------- Constructors ----------

    public LongIntHashMap() {
        this(1 << 20); // ~1 M slots
    }

    public LongIntHashMap(int initialCapacity) {
        capacity  = nextPowerOfTwo(initialCapacity);
        keys      = new long[capacity];
        values    = new int[capacity];
        threshold = (int) (capacity * 0.75);
        shift     = 64 - Integer.numberOfTrailingZeros(capacity);
        java.util.Arrays.fill(keys, EMPTY_KEY);
    }

    // ---------- Factories ----------

    /**
     * Pre-sizes the map for an expected number of distinct keys so the table
     * doesn't resize while filling. Capacity = smallest power of two whose 0.75
     * load threshold covers expectedKeys, capped at MAX_CAPACITY (~12 GiB arrays).
     */
    public static LongIntHashMap withExpectedKeys(long expectedKeys) {
        return new LongIntHashMap(capacityForExpectedKeys(expectedKeys));
    }

    /**
     * Pre-sizes for a sequence of {@code sequenceLength} bases. Distinct k-mers
     * can't exceed the number of positions (≈ sequence length), so that is the
     * upper bound. Multi-Gb inputs saturate at MAX_CAPACITY and rely on prune/cap.
     */
    public static LongIntHashMap forSequenceLength(long sequenceLength) {
        return withExpectedKeys(sequenceLength);
    }

    /** Smallest power-of-two capacity (≤ MAX_CAPACITY) whose 0.75 threshold ≥ expectedKeys. */
    private static int capacityForExpectedKeys(long expectedKeys) {
        if (expectedKeys < 1) return 2;
        long needed = (long) Math.ceil(expectedKeys / 0.75);  // slots required at load 0.75
        if (needed >= MAX_CAPACITY) return MAX_CAPACITY;
        int cap = Integer.highestOneBit((int) needed);
        if (cap < needed) cap <<= 1;                          // round up to next power of two
        return Math.max(2, Math.min(cap, MAX_CAPACITY));
    }

    // ---------- Public API ----------

    /**
     * Insert or update a key→value mapping. value must be ≥ 1.
     * Keys equal to the internal sentinel values are silently ignored.
     */
    public void put(long key, int value) {
        if (key == EMPTY_KEY || key == DELETED_KEY) return;
        if (value < 1) return;

        ensureCapacity();

        int  i            = index(key);
        int  firstDeleted = -1;   // first tombstone slot seen during probe

        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                // Key already present — update in place (no size change).
                values[i] = value;
                return;
            }
            if (keys[i] == DELETED_KEY && firstDeleted == -1) {
                firstDeleted = i; // remember but keep scanning for the key
            }
            i = next(i);
        }

        // Key is absent.  Prefer reusing a tombstone slot to keep chains short.
        if (firstDeleted != -1) {
            i = firstDeleted;
        }
        keys[i]   = key;
        values[i] = value;
        size++;
    }

    /**
     * Increments the value for {@code key} by 1, but never above {@code cap};
     * inserts the key with value 1 if it is absent.  This is equivalent to a
     * get-then-put pair but resolves the slot in a single probe sequence, which
     * matters in tight counting loops.  {@code cap} should be ≥ 1.  Sentinel
     * keys are silently ignored.
     *
     * @return the value stored for the key after the operation, or NO_VALUE if
     *         the key was a sentinel and was therefore ignored
     */
    public int incrementCapped(long key, int cap) {
        if (key == EMPTY_KEY || key == DELETED_KEY) return NO_VALUE;

        ensureCapacity(); // may resize/rehash; index() is computed afterwards

        int i            = index(key);
        int firstDeleted = -1;   // first tombstone slot seen during probe

        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                // Key present — bump the count unless it already hit the cap.
                if (values[i] < cap) {
                    values[i]++;
                }
                return values[i];
            }
            if (keys[i] == DELETED_KEY && firstDeleted == -1) {
                firstDeleted = i; // remember but keep scanning for the key
            }
            i = next(i);
        }

        // Key absent — insert with count 1, reusing a tombstone if seen.
        if (firstDeleted != -1) {
            i = firstDeleted;
        }
        keys[i]   = key;
        values[i] = 1;
        size++;
        return 1;
    }

    /**
     * Returns the value for the key, or NO_VALUE (-1) if absent.
     * Probes past tombstones correctly.
     */
    public int get(long key) {
        int i = index(key);
        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                return values[i];
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    /** Returns true if the key exists in the map. */
    public boolean containsKey(long key) {
        return get(key) != NO_VALUE;
    }

    /**
     * Removes the key. Returns the deleted value, or NO_VALUE if absent.
     */
    public int remove(long key) {
        int i = index(key);
        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                int old = values[i];
                keys[i]   = DELETED_KEY; // tombstone
                values[i] = 0;
                size--;
                return old;
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    /**
     * Drops every live entry whose value is &lt; {@code min} (tombstones it) and
     * rehashes ONCE to remove the resulting tombstone chains.
     *
     * Intended use: call after the counting phase is complete and before the
     * read-heavy phase. In k-mer repeat masking the unique (count == 1) k-mers
     * are the large majority, so retainAtLeast(MIN_COUNT) shrinks the table
     * many-fold — the subsequent get() calls then hit a small, cache-resident
     * table and most of the memory is freed. Lossless for repeat detection,
     * because once counting is finished a count-1 entry is genuinely unique.
     *
     * @return number of entries removed
     */
    public int retainAtLeast(int min) {
        int removed = 0;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY && values[i] < min) {
                keys[i]   = DELETED_KEY;
                values[i] = 0;
                size--;
                removed++;
            }
        }
        if (removed > 0) {
            rehash();
        }
        return removed;
    }

    public int size()      { return size; }
    public boolean isEmpty() { return size == 0; }

    public void clear() {
        java.util.Arrays.fill(keys, EMPTY_KEY);
        java.util.Arrays.fill(values, 0);
        size           = 0;
        pruneThreshold = 1;
    }

    /**
     * Merges all entries of this map into {@code target}, summing values and
     * capping at {@code maxValue}.  Used to merge per-thread maps.
     */
    public void mergeInto(LongIntHashMap target, int maxValue) {
        for (int i = 0; i < capacity; i++) {
            long k = keys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;

            int existing = target.get(k);
            if (existing == NO_VALUE) {
                target.put(k, Math.min(values[i], maxValue));
            } else {
                target.put(k, Math.min(existing + values[i], maxValue));
            }
        }
    }

    // ---------- Private helpers ----------

    /** Grows or prunes the table if we are at/above the load threshold. */
    private void ensureCapacity() {
        if (size < threshold) return;

        if (capacity < MAX_CAPACITY) {
            resize();
            return;
        }

        // At max capacity: prune low-frequency entries instead.
        int removed = pruneBelowOrEqual(pruneThreshold);
        if (removed == 0) {
            pruneThreshold++;
            removed = pruneBelowOrEqual(pruneThreshold);
        } else {
            pruneThreshold = 1; // reset after a successful cleanup
        }
        if (removed > 0) {
            rehash();
        }
        if (size >= threshold) {
            throw new IllegalStateException(
                    "Cannot free space. size=" + size
                    + ", pruneThreshold=" + pruneThreshold);
        }
    }

    /** Fibonacci hashing: take the top log2(capacity) bits of key·GOLDEN. */
    private int index(long key) {
        return (int) ((key * GOLDEN) >>> shift);
    }

    private int next(int i) {
        return (i + 1) & (capacity - 1);
    }

    /** Doubles the table capacity. */
    private void resize() {
        final int    oldCap    = capacity;
        final long[] oldKeys   = keys;
        final int[]  oldValues = values;

        final int newCapacity = oldCap << 1;
        final long[] newKeys   = new long[newCapacity];
        final int[]  newValues = new int[newCapacity];
        java.util.Arrays.fill(newKeys, EMPTY_KEY);

        // Publish the new geometry first so index()/next() use the new shift/cap.
        capacity  = newCapacity;
        threshold = (int) (newCapacity * 0.75);
        shift     = 64 - Integer.numberOfTrailingZeros(newCapacity);
        keys      = newKeys;
        values    = newValues;

        for (int i = 0; i < oldCap; i++) {
            long k = oldKeys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;

            int j = index(k);
            while (newKeys[j] != EMPTY_KEY) {
                j = next(j);
            }
            newKeys[j]   = k;
            newValues[j] = oldValues[i];
        }
    }

    /**
     * Rebuilds the table at the same capacity to eliminate tombstone chains.
     * Called after a prune to restore O(1) average probe length.
     */
    private void rehash() {
        long[] newKeys   = new long[capacity];
        int[]  newValues = new int[capacity];
        java.util.Arrays.fill(newKeys, EMPTY_KEY);

        for (int i = 0; i < capacity; i++) {
            long k = keys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;

            int j = index(k);
            while (newKeys[j] != EMPTY_KEY) {
                j = next(j);
            }
            newKeys[j]   = k;
            newValues[j] = values[i];
        }

        keys   = newKeys;
        values = newValues;
    }

    /**
     * Marks entries with value ≤ {@code minValue} as deleted.
     *
     * @return number of entries removed
     */
    private int pruneBelowOrEqual(int minValue) {
        int removed = 0;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY
                    && values[i] <= minValue) {
                keys[i]   = DELETED_KEY;
                values[i] = 0;
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

    // ---------- toString ----------

    @Override
    public String toString() {
        StringBuilder sb    = new StringBuilder("{");
        boolean       first = true;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY) {
                if (!first) sb.append(", ");
                sb.append(keys[i]).append('=').append(values[i]);
                first = false;
            }
        }
        return sb.append('}').toString();
    }
}
