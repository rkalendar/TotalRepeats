/**
 * Hash table with open addressing, specialised for capped COUNTING.
 * Key: long, Value: unsigned byte (0..255).
 *
 * Same layout/algorithm as LongIntHashMap, but the value array is byte[] instead
 * of int[]: in k-mer repeat masking the value is only a count capped at MIN_COUNT
 * (=2), so one byte per slot is plenty. At MAX_CAPACITY this is 1 GiB of values
 * instead of 4 GiB (keys stay 8 GiB), i.e. ~9 GiB total instead of ~12 GiB, and
 * the denser value array improves cache locality during the read phase.
 *
 * Values are read UNSIGNED (value & 0xFF), so the usable count range is 0..255.
 * {@code cap} passed to incrementCapped is clamped to 255 to avoid byte wrap.
 *
 * Sentinel values:
 *   EMPTY_KEY   = Long.MIN_VALUE       (never a real key)
 *   DELETED_KEY = Long.MIN_VALUE + 1   (tombstone after remove)
 */
public class LongByteHashMap {

    public static final int NO_VALUE = -1;
    private static final long EMPTY_KEY   = Long.MIN_VALUE;
    private static final long DELETED_KEY = Long.MIN_VALUE + 1;
    private static final int  MAX_CAPACITY = 1 << 30; // 1_073_741_824
    private static final int  MAX_VALUE    = 0xFF;     // unsigned byte ceiling

    /** 2^64 / golden ratio, odd — the standard Fibonacci-hashing multiplier. */
    private static final long GOLDEN = 0x9E3779B97F4A7C15L;

    private int pruneThreshold = 1; // start by deleting entries with value == 1

    private long[] keys;
    private byte[] values;          // unsigned count 0..255
    private int    size;
    private int    capacity;
    private int    threshold;
    private int    shift;           // = 64 - log2(capacity)

    // ---------- Constructors ----------

    public LongByteHashMap() {
        this(1 << 20); // ~1 M slots
    }

    public LongByteHashMap(int initialCapacity) {
        capacity  = nextPowerOfTwo(initialCapacity);
        keys      = new long[capacity];
        values    = new byte[capacity];
        threshold = (int) (capacity * 0.75);
        shift     = 64 - Integer.numberOfTrailingZeros(capacity);
        java.util.Arrays.fill(keys, EMPTY_KEY);
    }

    // ---------- Factories ----------

    /**
     * Pre-sizes the map for an expected number of distinct keys so the table
     * doesn't resize while filling. Capacity = smallest power of two whose 0.75
     * load threshold covers expectedKeys, capped at MAX_CAPACITY (~9 GiB arrays).
     */
    public static LongByteHashMap withExpectedKeys(long expectedKeys) {
        return new LongByteHashMap(capacityForExpectedKeys(expectedKeys));
    }

    /**
     * Pre-sizes for a sequence of {@code sequenceLength} bases. Distinct k-mers
     * can't exceed the number of positions (≈ sequence length), so that is the
     * upper bound. Multi-Gb inputs saturate at MAX_CAPACITY and rely on prune/cap.
     */
    public static LongByteHashMap forSequenceLength(long sequenceLength) {
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
     * Insert or update a key→value mapping. value is clamped to [1, 255].
     * Keys equal to the internal sentinel values are silently ignored.
     */
    public void put(long key, int value) {
        if (key == EMPTY_KEY || key == DELETED_KEY) return;
        if (value < 1) return;
        final byte v = (byte) Math.min(value, MAX_VALUE);

        ensureCapacity();

        int i            = index(key);
        int firstDeleted = -1;

        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                values[i] = v;        // update in place (no size change)
                return;
            }
            if (keys[i] == DELETED_KEY && firstDeleted == -1) {
                firstDeleted = i;
            }
            i = next(i);
        }

        if (firstDeleted != -1) {
            i = firstDeleted;
        }
        keys[i]   = key;
        values[i] = v;
        size++;
    }

    /**
     * Increments the count for {@code key} by 1, never above {@code cap} (and
     * never above 255); inserts the key with count 1 if absent. Resolves the
     * slot in a single probe sequence. Sentinel keys are silently ignored.
     *
     * @return the count stored after the operation, or NO_VALUE if the key was a
     *         sentinel and was therefore ignored
     */
    public int incrementCapped(long key, int cap) {
        if (key == EMPTY_KEY || key == DELETED_KEY) return NO_VALUE;
        final int eff = Math.min(cap, MAX_VALUE);   // guard against byte wrap

        ensureCapacity();

        int i            = index(key);
        int firstDeleted = -1;

        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                int v = values[i] & 0xFF;
                if (v < eff) {
                    v++;
                    values[i] = (byte) v;
                }
                return v;
            }
            if (keys[i] == DELETED_KEY && firstDeleted == -1) {
                firstDeleted = i;
            }
            i = next(i);
        }

        if (firstDeleted != -1) {
            i = firstDeleted;
        }
        keys[i]   = key;
        values[i] = 1;
        size++;
        return 1;
    }

    /**
     * Returns the (unsigned) count for the key, or NO_VALUE (-1) if absent.
     * A present count is always 1..255, so it never collides with NO_VALUE.
     */
    public int get(long key) {
        int i = index(key);
        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                return values[i] & 0xFF;
            }
            i = next(i);
        }
        return NO_VALUE;
    }

    /** Returns true if the key exists in the map. */
    public boolean containsKey(long key) {
        return get(key) != NO_VALUE;
    }

    /** Removes the key. Returns the deleted count, or NO_VALUE if absent. */
    public int remove(long key) {
        int i = index(key);
        while (keys[i] != EMPTY_KEY) {
            if (keys[i] == key) {
                int old = values[i] & 0xFF;
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
     * Drops every live entry whose count is &lt; {@code min} (tombstones it) and
     * rehashes ONCE to remove the resulting tombstone chains. Call after counting
     * and before the read-heavy phase: unique (count==1) k-mers are the large
     * majority, so retainAtLeast(MIN_COUNT) shrinks the table many-fold, making
     * subsequent get() calls cache-resident. Lossless once counting is finished.
     *
     * @return number of entries removed
     */
    public int retainAtLeast(int min) {
        int removed = 0;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY && (values[i] & 0xFF) < min) {
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

    public int size()       { return size; }
    public boolean isEmpty() { return size == 0; }

    public void clear() {
        java.util.Arrays.fill(keys, EMPTY_KEY);
        java.util.Arrays.fill(values, (byte) 0);
        size           = 0;
        pruneThreshold = 1;
    }

    /**
     * Merges all entries of this map into {@code target}, summing counts and
     * capping at {@code maxValue}. Used to merge per-thread maps.
     */
    public void mergeInto(LongByteHashMap target, int maxValue) {
        for (int i = 0; i < capacity; i++) {
            long k = keys[i];
            if (k == EMPTY_KEY || k == DELETED_KEY) continue;

            int mine     = values[i] & 0xFF;
            int existing = target.get(k);
            if (existing == NO_VALUE) {
                target.put(k, Math.min(mine, maxValue));
            } else {
                target.put(k, Math.min(existing + mine, maxValue));
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
            pruneThreshold = 1;
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
        final byte[] oldValues = values;

        final int newCapacity  = oldCap << 1;
        final long[] newKeys    = new long[newCapacity];
        final byte[] newValues  = new byte[newCapacity];
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

    /** Rebuilds the table at the same capacity to eliminate tombstone chains. */
    private void rehash() {
        long[] newKeys   = new long[capacity];
        byte[] newValues = new byte[capacity];
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
     * Marks entries with count ≤ {@code minValue} as deleted.
     *
     * @return number of entries removed
     */
    private int pruneBelowOrEqual(int minValue) {
        int removed = 0;
        for (int i = 0; i < capacity; i++) {
            if (keys[i] != EMPTY_KEY && keys[i] != DELETED_KEY
                    && (values[i] & 0xFF) <= minValue) {
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
                sb.append(keys[i]).append('=').append(values[i] & 0xFF);
                first = false;
            }
        }
        return sb.append('}').toString();
    }
}
