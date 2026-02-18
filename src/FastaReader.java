import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Reads sequences from FASTA-formatted input (files, byte arrays, or strings).
 * <p>
 * By default, nucleotide characters are normalized to lowercase and certain
 * IUPAC equivalences are applied (e.g. U→t, I→c). Use {@link #withoutNormalization()}
 * to preserve the original casing (soft-masking awareness).
 * </p>
 *
 * <pre>
 *   FastaReader reader = FastaReader.fromPath(path);
 *   List&lt;String&gt; names     = reader.getNames();
 *   List&lt;String&gt; sequences = reader.getSequences();
 * </pre>
 */
public final class FastaReader {

    // ── Static nucleotide lookup table (shared, immutable) ──────────────

    /** Maps ASCII code points to their normalized lowercase nucleotide value, or 0 if invalid. */
    private static final byte[] NUCLEOTIDE_TABLE = buildNucleotideTable();

    private static byte[] buildNucleotideTable() {
        byte[] table = new byte[128];

        // Uppercase → normalized lowercase
        table['A'] = 'a';
        table['B'] = 'b';
        table['C'] = 'c';
        table['D'] = 'd';
        table['G'] = 'g';
        table['H'] = 'h';
        table['I'] = 'c';   // Inosine → c
        table['K'] = 'k';
        table['M'] = 'm';
        table['N'] = 'n';
        table['R'] = 'r';
        table['S'] = 's';
        table['T'] = 't';
        table['U'] = 't';   // Uracil  → t
        table['V'] = 'v';
        table['W'] = 'w';
        table['Y'] = 'y';

        // Lowercase → normalized lowercase (identity, except i/u)
        table['a'] = 'a';
        table['b'] = 'b';
        table['c'] = 'c';
        table['d'] = 'd';
        table['g'] = 'g';
        table['h'] = 'h';
        table['i'] = 'c';   // inosine → c
        table['k'] = 'k';
        table['m'] = 'm';
        table['n'] = 'n';
        table['r'] = 'r';
        table['s'] = 's';
        table['t'] = 't';
        table['u'] = 't';   // uracil  → t
        table['v'] = 'v';
        table['w'] = 'w';
        table['y'] = 'y';

        return table;
    }

    // ── Instance state ──────────────────────────────────────────────────

    private final List<String> names;
    private final List<String> sequences;
    private final int totalResidues;

    private FastaReader(List<String> names, List<String> sequences, int totalResidues) {
        this.names = names;
        this.sequences = sequences;
        this.totalResidues = totalResidues;
    }

    // ── Public factory methods ──────────────────────────────────────────

    /** Parse a FASTA file with nucleotide normalization. */
    public static FastaReader fromPath(Path fastaPath) throws IOException {
        try (BufferedReader br = Files.newBufferedReader(fastaPath, StandardCharsets.US_ASCII)) {
            return parse(br, true);
        }
    }

    /** Parse a FASTA file, preserving original casing (e.g. for soft-masking). */
    public static FastaReader fromPathRaw(Path fastaPath) throws IOException {
        try (BufferedReader br = Files.newBufferedReader(fastaPath, StandardCharsets.US_ASCII)) {
            return parse(br, false);
        }
    }

    /** Parse FASTA content from a byte array with normalization. */
    public static FastaReader fromBytes(byte[] data) {
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(new ByteArrayInputStream(data), StandardCharsets.US_ASCII))) {
            return parse(br, true);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /** Parse FASTA content from a byte array without normalization. */
    public static FastaReader fromBytesRaw(byte[] data) {
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(new ByteArrayInputStream(data), StandardCharsets.US_ASCII))) {
            return parse(br, false);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /** Parse FASTA content from a String with normalization. */
    public static FastaReader fromString(String fasta) {
        try (BufferedReader br = new BufferedReader(new StringReader(fasta))) {
            return parse(br, true);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    // ── Accessors ───────────────────────────────────────────────────────

    /** Returns an unmodifiable list of sequence names (header lines without '>'). */
    public List<String> getNames() {
        return Collections.unmodifiableList(names);
    }

    /** Returns an unmodifiable list of sequences. */
    public List<String> getSequences() {
        return Collections.unmodifiableList(sequences);
    }

    /** Number of sequences read. */
    public int getSequenceCount() {
        return sequences.size();
    }

    /** Total number of valid residues across all sequences. */
    public int getTotalResidues() {
        return totalResidues;
    }

    /** True if no sequences were found. */
    public boolean isEmpty() {
        return sequences.isEmpty();
    }

    // ── Core parser ─────────────────────────────────────────────────────

    /**
     * Parses FASTA records from the given reader.
     *
     * @param br        buffered reader positioned at the start of FASTA content
     * @param normalize if true, map every residue through {@link #NUCLEOTIDE_TABLE};
     *                  if false, keep the original character (but still filter to valid nucleotides)
     */
    private static FastaReader parse(BufferedReader br, boolean normalize) throws IOException {
        List<String> names = new ArrayList<>();
        List<String> seqs = new ArrayList<>();
        int residueCount = 0;

        String line;
        String currentName = null;
        StringBuilder currentSeq = new StringBuilder(1 << 20); // ~1 MB initial capacity

        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                continue;
            }

            if (line.charAt(0) == '>') {
                if (currentName != null) {
                    seqs.add(currentSeq.toString());
                    currentSeq.setLength(0);
                }
                currentName = line.substring(1).trim();
                names.add(currentName);
                continue;
            }

            for (int i = 0, len = line.length(); i < len; i++) {
                char ch = line.charAt(i);
                if (ch < 128 && NUCLEOTIDE_TABLE[ch] != 0) {
                    currentSeq.append(normalize ? (char) NUCLEOTIDE_TABLE[ch] : ch);
                    residueCount++;
                }
            }
        }

        // Flush the last sequence
        if (currentName != null) {
            seqs.add(currentSeq.toString());
        }

        return new FastaReader(names, seqs, residueCount);
    }
}
