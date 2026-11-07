/**
 * The value-less presence set used for capped k-mer counting, where the
 * "seen at least twice" flag is packed into bit 49 of the slot itself.
 *
 * <p>The packing is what makes the k &le; 21 limit real, so the tests cover both
 * the saturating count and the {@link LongCountSet#fits} guard.
 */
final class LongCountSetTest {

    private LongCountSetTest() {
    }

    static void run() {
        TestRunner.suite("LongCountSet");

        // ── the k-mer bound the packing depends on ──
        TestRunner.eq("MAX_KMER", 21, LongCountSet.MAX_KMER);
        TestRunner.ok("fits(21)", LongCountSet.fits(21));
        TestRunner.ok("fits(19) (the default k)", LongCountSet.fits(19));
        TestRunner.ok("does not fit 22", !LongCountSet.fits(22));

        LongCountSet s = new LongCountSet(16);

        // ── absent ──
        TestRunner.eq("get on empty set", LongCountSet.NO_VALUE, s.get(123));
        TestRunner.ok("containsKey false when absent", !s.containsKey(123));
        TestRunner.ok("new set is empty", s.isEmpty());

        // ── counts saturate at 2, which is all the masker needs ──
        TestRunner.eq("first sighting", 1, s.incrementCapped(123, 2));
        TestRunner.eq("second sighting", 2, s.incrementCapped(123, 2));
        TestRunner.eq("third sighting stays at 2", 2, s.incrementCapped(123, 2));
        TestRunner.eq("get reports the capped count", 2, s.get(123));
        TestRunner.eq("size counts distinct codes", 1, s.size());

        s.incrementCapped(456, 2);
        TestRunner.eq("a once-seen code reports 1", 1, s.get(456));
        TestRunner.eq("two distinct codes", 2, s.size());

        // ── remove and probe past the tombstone ──
        TestRunner.eq("remove returns the count", 2, s.remove(123));
        TestRunner.eq("get after remove", LongCountSet.NO_VALUE, s.get(123));
        TestRunner.eq("other code still reachable", 1, s.get(456));

        // ── retainAtLeast(2) is the unique-k-mer prune of the masking phase ──
        LongCountSet p = new LongCountSet(64);
        p.incrementCapped(1, 2);                 // seen once
        p.incrementCapped(2, 2);                 // seen once
        p.incrementCapped(3, 2);
        p.incrementCapped(3, 2);                 // seen twice
        p.incrementCapped(4, 2);
        p.incrementCapped(4, 2);                 // seen twice
        int removed = p.retainAtLeast(2);
        TestRunner.eq("prune drops the singletons", 2, removed);
        TestRunner.eq("size after prune", 2, p.size());
        TestRunner.eq("singleton is gone", LongCountSet.NO_VALUE, p.get(1));
        TestRunner.eq("repeated code retained", 2, p.get(3));
        TestRunner.eq("repeated code retained (2)", 2, p.get(4));

        // ── growth past the initial capacity ──
        LongCountSet big = new LongCountSet(16);
        final int n = 5000;
        for (int i = 1; i <= n; i++) {
            big.incrementCapped(i, 2);
        }
        TestRunner.eq("size after growth", n, big.size());
        boolean allPresent = true;
        for (int i = 1; i <= n; i++) {
            if (big.get(i) != 1) {
                allPresent = false;
                break;
            }
        }
        TestRunner.ok("every code survives rehashing", allPresent);

        // ── large codes: a base-5 k-mer code for k=21 stays below the flag bit ──
        long bigCode = 1L << 47;
        LongCountSet wide = new LongCountSet(16);
        wide.incrementCapped(bigCode, 2);
        TestRunner.eq("code near the flag bit is stored", 1, wide.get(bigCode));
        wide.incrementCapped(bigCode, 2);
        TestRunner.eq("code near the flag bit counts to 2", 2, wide.get(bigCode));

        // ── mergeInto, as used to combine per-thread sets ──
        LongCountSet a = new LongCountSet(16);
        LongCountSet b = new LongCountSet(16);
        a.incrementCapped(10, 2);                // once in a
        b.incrementCapped(10, 2);                // once in b -> twice combined
        b.incrementCapped(11, 2);
        b.mergeInto(a);
        TestRunner.eq("merge lifts a shared code to 2", 2, a.get(10));
        TestRunner.eq("merge adds the unshared code", 1, a.get(11));

        // ── clear ──
        a.clear();
        TestRunner.ok("clear empties the set", a.isEmpty());
        TestRunner.eq("get after clear", LongCountSet.NO_VALUE, a.get(10));
    }
}
