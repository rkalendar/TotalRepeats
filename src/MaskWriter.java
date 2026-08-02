import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Writes the soft-masked FASTA files, in which repeat and short-tandem-repeat
 * positions are lowercased and everything else is left upper case.
 *
 * <p>Extracted from {@code TotalRepeatsSearching}. Only the writers that are
 * free of side effects on the analysis state live here; the two that also
 * update the owner ({@code SaveMask}, which back-fills the gap length, and
 * {@code lowercaseAndSaveMsk}, which repoints the current file path) stay with
 * it, so that state remains owned in one place.
 *
 * <p>A writer is built per file, capturing the path and sequence set in use at
 * that moment, exactly as the previous field reads did.
 */
final class MaskWriter {

    private final String filePath;
    private final String[] seq;
    private final String[] sname;
    private final int nseq;
    /** Repeat coverage of the current sequence, quoted in the FASTA header. */
    private final double repeatslen;

    MaskWriter(String filePath, String[] seq, String[] sname, int nseq, double repeatslen) {
        this.filePath = filePath;
        this.seq = seq;
        this.sname = sname;
        this.nseq = nseq;
        this.repeatslen = repeatslen;
    }

    void writeMask(int n, int[] m, int[] ssr) throws IOException {
        String maskedfile = filePath + "_" + (n + 1) + ".msk";
        if (nseq == 1) {
            maskedfile = filePath + ".msk";
        }
        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            System.out.println("Saving masked file: " + maskedfile);

            byte[] c = seq[n].toUpperCase().getBytes();
// UPPER letter to lower for repeats            
            for (int j = 0; j < m.length; j += 2) {
                for (int i = m[j]; i < m[j] + m[j + 1]; i++) {
                    if (c[i] > 64 && c[i] < 90) {
                        c[i] = (byte) (c[i] + 32);
                    }
                }
            }
//  SSR masking
            for (int j = 0; j < ssr.length; j += 2) {
                for (int i = ssr[j]; i < ssr[j] + ssr[j + 1]; i++) {
                    if (c[i] > 64 && c[i] < 90) {
                        c[i] = (byte) (c[i] + 32);
                    }
                }
            }

            fileWriter.write(">" + sname[n] + " TotalRepeats: Sequence coverage by repeats = " + String.format("%.2f", repeatslen) + "%\n");
//            fileWriter.write(new String(c));                   
            String seqStr = new String(c);
            for (int i = 0; i < seqStr.length(); i += 70) {
                int end = Math.min(i + 70, seqStr.length());
                fileWriter.write(seqStr.substring(i, end));
                fileWriter.write("\n");
            }

        }
    }
    void writeMaskTo(String maskedfile, int n, int[] m, int[] ssr) throws IOException {
        if (maskedfile.length() == 0) {
            maskedfile = filePath + ".msk";
        } else {
            maskedfile = maskedfile + ".msk";
        }
        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            System.out.println("Saving masked file: " + maskedfile);

            byte[] c = seq[n].toUpperCase().getBytes();
// UPPER letter to lower for repeats            
            for (int j = 0; j < m.length; j += 2) {
                for (int i = m[j]; i < m[j] + m[j + 1]; i++) {
                    if (c[i] > 64 && c[i] < 90) {
                        c[i] = (byte) (c[i] + 32);
                    }
                }
            }
//  SSR masking
            for (int j = 0; j < ssr.length; j += 2) {
                for (int i = ssr[j]; i < ssr[j] + ssr[j + 1]; i++) {
                    if (c[i] > 64 && c[i] < 90) {
                        c[i] = (byte) (c[i] + 32);
                    }
                }
            }

            fileWriter.write(">" + sname[n] + " TotalRepeats: Sequence coverage by repeats = " + String.format("%.2f", repeatslen) + "%\n");
//            fileWriter.write(new String(c));                   
            String seqStr = new String(c);
            for (int i = 0; i < seqStr.length(); i += 70) {
                int end = Math.min(i + 70, seqStr.length());
                fileWriter.write(seqStr.substring(i, end));
                fileWriter.write("\n");
            }

        }
    }
    void writeBytes(int n, byte[] c, int x1, int x2) throws IOException {
        String maskedfile = (nseq == 1) ? (filePath + ".msk") : (filePath + "_" + (n + 1) + ".msk");

        System.out.println("Saving masked file: " + maskedfile);
        Path out = Path.of(maskedfile);
        try (OutputStream os = Files.newOutputStream(out)) {
            os.write(('>' + sname[n] + " TotalRepeats mask\n").getBytes(StandardCharsets.US_ASCII));
            final int width = 70;
            int i = x1;
            while (i < x2) {
                int end = Math.min(i + width, x2);
                os.write(c, i, end - i);
                os.write('\n');
                i = end;
            }
        }
    }
    static void softMask(byte[] ci, int[] blocks, int l) {
        for (int j = 0; j + 1 < blocks.length; j += 2) {
            int from = Math.max(0, blocks[j]);
            int to = Math.min(l, blocks[j] + Math.abs(blocks[j + 1]));
            for (int p = from; p < to; p++) {
                byte bch = ci[p];
                if (bch >= 'A' && bch < 'Z') {
                    ci[p] = (byte) (bch + 32);
                }
            }
        }
    }}
