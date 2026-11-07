/**
 * The open-addressing long-&gt;int map used for k-mer counting: insertion,
 * tombstones, growth, and the prune/merge operations the masking phase relies on.
 */
final class LongIntHashMapTest {

    private LongIntHashMapTest() {
    }

    static void run() {
        TestRunner.suite("LongIntHashMap");

        LongIntHashMap m = new LongIntHashMap(16);

        // ── absent keys ──
        TestRunner.eq("get on empty map", LongIntHashMap.NO_VALUE, m.get(42));
        TestRunner.ok("containsKey false when absent", !m.containsKey(42));
        TestRunner.ok("new map is empty", m.isEmpty());

        // ── put / get ──
        m.put(42, 7);
        TestRunner.eq("get after put", 7, m.get(42));
        TestRunner.eq("size after put", 1, m.size());
        TestRunner.ok("containsKey true when present", m.containsKey(42));

        // update in place must not change size
        m.put(42, 9);
        TestRunner.eq("put overwrites", 9, m.get(42));
        TestRunner.eq("overwrite keeps size", 1, m.size());

        // ── values below 1 are rejected by contract ──
        m.put(43, 0);
        TestRunner.eq("value 0 is ignored", LongIntHashMap.NO_VALUE, m.get(43));

        // ── negative keys are ordinary keys ──
        m.put(-5, 3);
        TestRunner.eq("negative key", 3, m.get(-5));

        // ── remove leaves a tombstone that later probes must see past ──
        TestRunner.eq("remove returns old value", 9, m.remove(42));
        TestRunner.eq("get after remove", LongIntHashMap.NO_VALUE, m.get(42));
        TestRunner.eq("remove absent key", LongIntHashMap.NO_VALUE, m.remove(999));
        TestRunner.eq("size after remove", 1, m.size());
        TestRunner.eq("probing past a tombstone still finds -5", 3, m.get(-5));

        // ── incrementCapped: insert at 1, then saturate ──
        LongIntHashMap c = new LongIntHashMap(16);
        TestRunner.eq("first increment inserts 1", 1, c.incrementCapped(100, 3));
        TestRunner.eq("second increment", 2, c.incrementCapped(100, 3));
        TestRunner.eq("third increment reaches cap", 3, c.incrementCapped(100, 3));
        TestRunner.eq("increment saturates at cap", 3, c.incrementCapped(100, 3));
        TestRunner.eq("counted key size", 1, c.size());

        // ── growth: insert well past the 0.75 load factor of the initial table ──
        LongIntHashMap big = new LongIntHashMap(16);
        final int n = 5000;
        for (int i = 1; i <= n; i++) {
            big.put(i, i);
        }
        TestRunner.eq("size after growth", n, big.size());
        boolean allFound = true;
        for (int i = 1; i <= n; i++) {
            if (big.get(i) != i) {
                allFound = false;
                break;
            }
        }
        TestRunner.ok("every key survives rehashing", allFound);
        TestRunner.eq("absent key after growth", LongIntHashMap.NO_VALUE, big.get(n + 1));

        // ── retainAtLeast drops low counts, which is how unique k-mers are pruned ──
        LongIntHashMap p = new LongIntHashMap(64);
        p.put(1, 1);
        p.put(2, 1);
        p.put(3, 2);
        p.put(4, 5);
        int removed = p.retainAtLeast(2);
        TestRunner.eq("retainAtLeast removes count-1 entries", 2, removed);
        TestRunner.eq("size after prune", 2, p.size());
        TestRunner.eq("pruned key is gone", LongIntHashMap.NO_VALUE, p.get(1));
        TestRunner.eq("retained key keeps its value", 2, p.get(3));
        TestRunner.eq("retained key keeps its value (2)", 5, p.get(4));

        // ── mergeInto sums and caps, as used to combine per-thread maps ──
        LongIntHashMap a = new LongIntHashMap(16);
        LongIntHashMap b = new LongIntHashMap(16);
        a.put(10, 2);
        a.put(11, 1);
        b.put(10, 3);
        b.put(12, 4);
        b.mergeInto(a, 100);
        TestRunner.eq("merge sums shared key", 5, a.get(10));
        TestRunner.eq("merge keeps unshared key of target", 1, a.get(11));
        TestRunner.eq("merge adds unshared key of source", 4, a.get(12));

        LongIntHashMap capA = new LongIntHashMap(16);
        LongIntHashMap capB = new LongIntHashMap(16);
        capA.put(1, 8);
        capB.put(1, 9);
        capB.mergeInto(capA, 10);
        TestRunner.eq("merge respects the cap", 10, capA.get(1));

        // ── clear ──
        p.clear();
        TestRunner.ok("clear empties the map", p.isEmpty());
        TestRunner.eq("get after clear", LongIntHashMap.NO_VALUE, p.get(3));

        // ── pre-sizing factories ──
        LongIntHashMap sized = LongIntHashMap.withExpectedKeys(1000);
        sized.put(7, 1);
        TestRunner.eq("withExpectedKeys works", 1, sized.get(7));
        LongIntHashMap forSeq = LongIntHashMap.forSequenceLength(10_000);
        forSeq.put(7, 2);
        TestRunner.eq("forSequenceLength works", 2, forSeq.get(7));
    }
}
