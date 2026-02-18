import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.NumberFormat;
import java.util.Locale;

public class FileMasksComparator {

    public FileMasksComparator() {
    }

    private long gapslen = 0;
    private long repeatslen = 0;

    public long GapsLength() {
        return gapslen;
    }

    public long RepeatLength() {
        return repeatslen;
    }

    public StringBuilder AnalysisFiles(String inputFile1, String inputFile2) {
        Path p1 = Path.of(inputFile1);
        Path p2 = Path.of(inputFile2);

        StringBuilder sb = new StringBuilder();
        Locale locale = Locale.US;
        NumberFormat intFmt = NumberFormat.getIntegerInstance(locale);
        NumberFormat pctFmt = NumberFormat.getNumberInstance(locale);
        pctFmt.setMinimumFractionDigits(2);
        pctFmt.setMaximumFractionDigits(2);

        long sameLower = 0;
        long diffCase1 = 0;
        long diffCase2 = 0;
        long mask1 = 0;
        long mask2 = 0;
        long totalCompared = 0;

        String name1 = "";
        String name2 = "";
        long len1 = 0;
        long len2 = 0;

        try (BufferedReader br1 = Files.newBufferedReader(p1, StandardCharsets.US_ASCII); BufferedReader br2 = Files.newBufferedReader(p2, StandardCharsets.US_ASCII)) {
            FastaLetterIterator it1 = new FastaLetterIterator(br1);
            FastaLetterIterator it2 = new FastaLetterIterator(br2);
            name1 = it1.headerName();
            name2 = it2.headerName();

            int ch1, ch2;
            while (true) {
                ch1 = it1.nextLetter();
                ch2 = it2.nextLetter();
                if (ch1 == -1 || ch2 == -1) {
                    break;
                }

                boolean isLower1 = (ch1 >= 'a' && ch1 < 'z');
                boolean isLower2 = (ch2 >= 'a' && ch2 < 'z');

                if (isLower1) {
                    mask1++;
                }
                if (isLower2) {
                    mask2++;
                }

                if (isLower1 && isLower2) {
                    sameLower++;
                } else if (isLower1) {
                    diffCase1++;
                } else if (isLower2) {
                    diffCase2++;
                }
                totalCompared++;
            }

            len1 = it1.lettersSeen() + countRemaining(it1);
            len2 = it2.lettersSeen() + countRemaining(it2);

            if (len1 == len2) {
                sb.append("Sequence length (nt): ").append(intFmt.format(totalCompared)).append("\n");
            } else {
                sb.append("Overlapping sequence length (nt): ").append(intFmt.format(totalCompared))
                        .append(" (").append(intFmt.format(len1)).append("/")
                        .append(intFmt.format(len2)).append(")\n");
            }
            sb.append("(1) ").append(inputFile1).append(" : ").append(name1).append("\n");
            sb.append("(2) ").append(inputFile2).append(" : ").append(name2).append("\n");

            double p1masked = (len1 == 0) ? 0.0 : (mask1 * 100.0) / len1;
            double p2masked = (len2 == 0) ? 0.0 : (mask2 * 100.0) / len2;
            double poverlap = (Math.max(mask1, mask2) == 0) ? 0.0 : (sameLower * 100.0) / Math.max(mask1, mask2);
            double pseqcov = (totalCompared == 0) ? 0.0 : (sameLower * 100.0) / totalCompared;
            double pdiff1 = (mask1 == 0) ? 0.0 : (diffCase1 * 100.0) / mask1;
            double pdiff2 = (mask2 == 0) ? 0.0 : (diffCase2 * 100.0) / mask2;

            sb.append("(1) Masked (nt): ").append(intFmt.format(mask1)).append(" : ").append(pctFmt.format(p1masked)).append("%\n");
            sb.append("(2) Masked (nt): ").append(intFmt.format(mask2)).append(" : ").append(pctFmt.format(p2masked)).append("%\n");
            sb.append("Overlapping masking (nt): ").append(intFmt.format(sameLower)).append(" : ").append(pctFmt.format(poverlap)).append("%\n");
            sb.append("Overlapping coverage at sequence level = ").append(pctFmt.format(pseqcov)).append("%\n");
            sb.append("(1) Different mask not overlap (nt): ").append(intFmt.format(diffCase1)).append(" : ").append(pctFmt.format(pdiff1)).append("%\n");
            sb.append("(2) Different mask not overlap (nt): ").append(intFmt.format(diffCase2)).append(" : ").append(pctFmt.format(pdiff2)).append("%\n");

            sb.append(totalCompared).append("\t")
                    .append(mask1).append("\t")
                    .append(mask2).append("\t")
                    .append(sameLower).append("\t")
                    .append(diffCase1).append("\t")
                    .append(diffCase2).append("\n");

            this.repeatslen = sameLower;
            this.gapslen = diffCase1 + diffCase2;

        } catch (IOException e) {
            sb.append("ERROR: ").append(e.getMessage()).append("\n");
        }

        return sb;
    }

    private static long countRemaining(FastaLetterIterator it) throws IOException {
        long c = 0;
        int ch;
        while ((ch = it.nextLetter()) != -1) {
            c++;
        }
        return c;
    }

    private static final class FastaLetterIterator {

        private final BufferedReader br;
        private String header = "";
        private String line = null;
        private int idx = 0;
        private boolean headerParsed = false;
        private long seen = 0L;

        FastaLetterIterator(BufferedReader br) throws IOException {
            this.br = br;
            parseHeader();
        }

        String headerName() {
            return header;
        }

        long lettersSeen() {
            return seen;
        }

        private void parseHeader() throws IOException {
            String first;
            while ((first = br.readLine()) != null) {
                if (first.isEmpty()) {
                    continue;
                }
                if (first.charAt(0) == '>') {
                    header = first.substring(1).trim();
                    break;
                } else {
                    line = first;
                    idx = 0;
                    header = "";
                    break;
                }
            }
            headerParsed = true;
        }

        int nextLetter() throws IOException {
            while (true) {
                if (line == null || idx >= line.length()) {
                    line = br.readLine();
                    idx = 0;
                    if (line == null) {
                        return -1;
                    }
                    if (!line.isEmpty() && line.charAt(0) == '>' && headerParsed) {
                        return -1;
                    }
                    continue;
                }
                char ch = line.charAt(idx++);
                if ((ch >= 'A' && ch < 'Z') || (ch >= 'a' && ch < 'z')) {
                    seen++;
                    return ch;
                }

            }
        }
    }
}
