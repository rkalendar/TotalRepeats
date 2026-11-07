import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Minimal, dependency-free test harness for TotalRepeats.
 *
 * <p>TotalRepeats ships as a single JAR with no external dependencies, and the
 * NetBeans metadata that would carry a JUnit configuration is not tracked in the
 * repository. The test suite therefore uses nothing but the JDK: compile the
 * sources together with this directory and run {@code TestRunner}.
 *
 * <pre>
 *   ./run-tests.sh          (Linux/macOS)
 *   run-tests.bat           (Windows)
 * </pre>
 *
 * <p>The process exits with status 1 if any check fails, so the suite can gate a
 * build or a CI job. Test classes live in the default package, as the code under
 * test does, which lets them reach package-private members directly.
 */
public final class TestRunner {

    private static int passed;
    private static final List<String> FAILURES = new ArrayList<>();
    private static String suite = "(none)";

    private TestRunner() {
    }

    // ── Suite bookkeeping ───────────────────────────────────────────────────

    /** Starts a named group of checks; the name is printed with any failure. */
    public static void suite(String name) {
        suite = name;
        System.out.println("── " + name);
    }

    private static void fail(String what, String detail) {
        FAILURES.add(suite + " :: " + what + "\n      " + detail);
        System.out.println("   FAIL  " + what + "\n         " + detail);
    }

    private static void pass(String what) {
        passed++;
        System.out.println("   ok    " + what);
    }

    // ── Assertions ──────────────────────────────────────────────────────────

    public static void ok(String what, boolean condition) {
        if (condition) {
            pass(what);
        } else {
            fail(what, "expected true");
        }
    }

    public static void eq(String what, long expected, long actual) {
        if (expected == actual) {
            pass(what);
        } else {
            fail(what, "expected " + expected + " but was " + actual);
        }
    }

    public static void eq(String what, String expected, String actual) {
        if (expected == null ? actual == null : expected.equals(actual)) {
            pass(what);
        } else {
            fail(what, "expected \"" + expected + "\" but was \"" + actual + "\"");
        }
    }

    public static void eq(String what, int[] expected, int[] actual) {
        if (Arrays.equals(expected, actual)) {
            pass(what);
        } else {
            fail(what, "expected " + Arrays.toString(expected)
                    + " but was " + Arrays.toString(actual));
        }
    }

    public static void eq(String what, long[] expected, long[] actual) {
        if (Arrays.equals(expected, actual)) {
            pass(what);
        } else {
            fail(what, "expected " + Arrays.toString(expected)
                    + " but was " + Arrays.toString(actual));
        }
    }

    /** Passes when {@code body} throws {@code expected} (or a subclass). */
    public static void throwsEx(String what, Class<? extends Throwable> expected, Runnable body) {
        try {
            body.run();
        } catch (Throwable t) {
            if (expected.isInstance(t)) {
                pass(what);
            } else {
                fail(what, "expected " + expected.getSimpleName()
                        + " but threw " + t.getClass().getSimpleName() + ": " + t.getMessage());
            }
            return;
        }
        fail(what, "expected " + expected.getSimpleName() + " but nothing was thrown");
    }

    /** Passes when {@code body} completes without throwing. */
    public static void noThrow(String what, Runnable body) {
        try {
            body.run();
            pass(what);
        } catch (Throwable t) {
            fail(what, "threw " + t.getClass().getSimpleName() + ": " + t.getMessage());
        }
    }

    // ── Entry point ─────────────────────────────────────────────────────────

    /** A suite body, allowed to throw so tests can call code that declares exceptions. */
    private interface Suite {
        void run() throws Exception;
    }

    /**
     * Runs one suite. An exception escaping a suite is recorded as a failure and
     * the remaining suites still run — otherwise a single unexpected throw (which
     * is itself a defect worth reporting) would hide every later result.
     */
    private static void guard(String name, Suite body) {
        try {
            body.run();
        } catch (Throwable t) {
            suite = name;
            StringBuilder sb = new StringBuilder(t.getClass().getName());
            if (t.getMessage() != null) {
                sb.append(": ").append(t.getMessage());
            }
            StackTraceElement[] st = t.getStackTrace();
            for (int i = 0; i < st.length && i < 3; i++) {
                sb.append("\n        at ").append(st[i]);
            }
            fail("suite aborted by an unexpected exception", sb.toString());
        }
    }

    public static void main(String[] args) {
        long t0 = System.nanoTime();
        System.out.println("TotalRepeats test suite");
        System.out.println();

        guard("Dna", DnaTest::run);
        guard("FastaReader", FastaReaderTest::run);
        guard("SeqStore", SeqStoreTest::run);
        guard("LongIntHashMap", LongIntHashMapTest::run);
        guard("LongCountSet", LongCountSetTest::run);
        guard("MaskResult", MaskResultTest::run);
        guard("PangenomeAnalysis", PangenomeAnalysisTest::run);
        guard("CLI", CliTest::run);

        long ms = (System.nanoTime() - t0) / 1_000_000;
        System.out.println();
        if (FAILURES.isEmpty()) {
            System.out.println("PASSED  " + passed + " checks in " + ms + " ms");
        } else {
            System.out.println("FAILED  " + FAILURES.size() + " of "
                    + (passed + FAILURES.size()) + " checks in " + ms + " ms");
            System.out.println();
            for (String f : FAILURES) {
                System.out.println("  " + f);
            }
            System.exit(1);
        }
    }
}
