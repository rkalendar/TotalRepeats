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

    /**
     * Name of the single record read from input with no header line at all;
     * a file read by path is named after the file instead.
     */
    public static final String PLAIN_NAME = "sequence";

    /**
     * Share of the letters of header-less input that must be nucleotide codes
     * for it to be read as a plain sequence. Prose falls well short — about a
     * third of English letters (e, o, l, p, f, ...) are no IUPAC code — so a
     * stray text file in an input folder is not analysed as DNA.
     */
    private static final double PLAIN_MIN_SHARE = 0.9;

    /** Parse a FASTA file with nucleotide normalization. */
    public static FastaReader fromPath(Path fastaPath) throws IOException {
        try (BufferedReader br = Files.newBufferedReader(fastaPath, StandardCharsets.US_ASCII)) {
            return parse(br, true, fastaPath.getFileName().toString());
        }
    }

    /** Parse a FASTA file, preserving original casing (e.g. for soft-masking). */
    public static FastaReader fromPathRaw(Path fastaPath) throws IOException {
        try (BufferedReader br = Files.newBufferedReader(fastaPath, StandardCharsets.US_ASCII)) {
            return parse(br, false, fastaPath.getFileName().toString());
        }
    }

    /** Parse FASTA content from a byte array with normalization. */
    public static FastaReader fromBytes(byte[] data) {
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(new ByteArrayInputStream(data), StandardCharsets.US_ASCII))) {
            return parse(br, true, PLAIN_NAME);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /** Parse FASTA content from a byte array without normalization. */
    public static FastaReader fromBytesRaw(byte[] data) {
        try (BufferedReader br = new BufferedReader(
                new InputStreamReader(new ByteArrayInputStream(data), StandardCharsets.US_ASCII))) {
            return parse(br, false, PLAIN_NAME);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /** Parse FASTA content from a String with normalization. */
    public static FastaReader fromString(String fasta) {
        try (BufferedReader br = new BufferedReader(new StringReader(fasta))) {
            return parse(br, true, PLAIN_NAME);
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
     * <p>Text before the first header belongs to no record and is ignored; it
     * used to be prepended to the first record, shifting all its coordinates.
     * Input with no header at all is read as one plain sequence named
     * {@code plainName}, provided it looks like one (see {@link #PLAIN_MIN_SHARE}).
     *
     * @param br        buffered reader positioned at the start of FASTA content
     * @param normalize if true, map every residue through {@link #NUCLEOTIDE_TABLE};
     *                  if false, keep the original character (but still filter to valid nucleotides)
     * @param plainName name of the record read from header-less input
     */
    private static FastaReader parse(BufferedReader br, boolean normalize, String plainName) throws IOException {
        List<String> names = new ArrayList<>();
        List<String> seqs = new ArrayList<>();
        int residueCount = 0;
        long plainLetters = 0;   // letters before the first header, nucleotide codes or not

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
                } else {
                    residueCount -= currentSeq.length();   // text before the first header
                }
                currentSeq.setLength(0);
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
                if (currentName == null && Character.isLetter(ch)) {
                    plainLetters++;
                }
            }
        }

        // Flush the last sequence
        if (currentName != null) {
            seqs.add(currentSeq.toString());
        } else if (currentSeq.length() > 0 && currentSeq.length() >= PLAIN_MIN_SHARE * plainLetters) {
            names.add(plainName);
            seqs.add(currentSeq.toString());
        } else {
            residueCount = 0;
        }

        return new FastaReader(names, seqs, residueCount);
    }
}
