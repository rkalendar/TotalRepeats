import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Writes the repeat annotation table and the pangenome report for one run.
 *
 * <p>Extracted from {@code TotalRepeatsSearching}, where these writers read
 * nineteen mutable instance fields of the enclosing class. Those fields were in
 * effect undeclared parameters: they were reassigned as a run walked its input
 * files, so the output of a writer depended on when it happened to be called.
 * Here the same values arrive through the constructor, which makes the contract
 * explicit and the writers reusable.
 *
 * <p>The dependencies are grouped in two records to keep that contract readable:
 * {@link Inputs} carries the sequences, names and clusters being reported, and
 * {@link Stats} the run parameters and coverage figures that go into the header.
 * A writer is built per report, so it captures the statistics as they stood at
 * that moment — matching the previous behaviour exactly.
 *
 * <p>The drawing of coordinates is unchanged from the original writers, including
 * the deliberate int/long split: {@link #writeTable} addresses one sequence with
 * int coordinates, and {@link #writeTableLong} the whole virtual concatenation
 * with long ones, so a combined run above ~2.1 Gb is still reportable.
 */
final class AnnotationWriter {

    /** The data being reported. */
    record Inputs(String[] seq, String[] sname, String[] filesPath, String filePath,
            int nseq, ArrayList<int[]> bb, int[] refclust, String[] refsname) {
    }

    /** Run parameters and coverage statistics, as printed in the report header. */
    record Stats(int kmerln, int minlenseq, int flanks, int gap, boolean seqShow,
            double repeatslen, double ssrglobal, double gapslen, double gaps,
            long maskduration, long startTime, boolean pangenome) {
    }

    // Unpacked into fields whose names match the originals, so the writer bodies
    // below are the same code that used to live in TotalRepeatsSearching.
    private final String[] seq;
    private final String[] sname;
    private final String[] filesPath;
    private final String filePath;
    private final int nseq;
    private final ArrayList<int[]> bb;
    private final int[] refclust;
    private final String[] refsname;

    private final int kmerln;
    private final int minlenseq;
    private final int flanks;
    private final int gap;
    private final boolean SeqShow;
    private final double repeatslen;
    private final double ssrglobal;
    private final double gapslen;
    private final double gaps;
    private final long maskduration;
    private final long startTime;
    private final boolean pangenome;

    AnnotationWriter(Inputs in, Stats st) {
        this.seq = in.seq();
        this.sname = in.sname();
        this.filesPath = in.filesPath();
        this.filePath = in.filePath();
        this.nseq = in.nseq();
        this.bb = in.bb();
        this.refclust = in.refclust();
        this.refsname = in.refsname();

        this.kmerln = st.kmerln();
        this.minlenseq = st.minlenseq();
        this.flanks = st.flanks();
        this.gap = st.gap();
        this.SeqShow = st.seqShow();
        this.repeatslen = st.repeatslen();
        this.ssrglobal = st.ssrglobal();
        this.gapslen = st.gapslen();
        this.gaps = st.gaps();
        this.maskduration = st.maskduration();
        this.startTime = st.startTime();
        this.pangenome = st.pangenome();
    }

    void writeTable(String reportfile, int n, int l, int[] h) throws IOException {
        String b = sname[n];
        long duration = (System.nanoTime() - startTime) / 1000000000;

        if (reportfile.length() == 0) {
            reportfile = filePath + "_" + (n + 1) + ".gff";
            if (h.length > 0) {
                reportfile = filePath + ".gff";
            }
        } else {
            reportfile = reportfile + ".gff";
        }

        try (FileWriter fileWriter = new FileWriter(reportfile); BufferedWriter bufferedWriter = new BufferedWriter(fileWriter)) {
            System.out.println("Saving report file: " + reportfile);
            StringBuilder sr = new StringBuilder();
            sr.append("#TotalRepeats (2024-2026) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi) https://github.com/rkalendar/TotalRepeats\n");
            sr.append("#kmer=").append(kmerln).append("\n").append("#Minimal repeat block size=").append(minlenseq).append("\n");
            sr.append("#Sequence length (bp)=").append(l).append("\n");
            sr.append("#Sequence coverage by repeats=").append(String.format("%.2f", repeatslen)).append("%\n");
            sr.append("#Short tandem repeat (STR) sequence coverage=").append(String.format("%.2f", ssrglobal)).append("%\n");
            sr.append("#Sequence gap (bp)=").append((int) gapslen).append(" (").append(String.format("%.4f", gaps)).append("%)\n");
            sr.append("#Masking time taken: ").append(maskduration).append(" seconds\n");
            sr.append("#Total duration: ").append(duration).append(" seconds\n");
            sr.append("#Repeats search for: ");

            if (h.length > 0) {
                for (String filesPath1 : filesPath) {
                    sr.append(filesPath1).append("\n");
                }
            } else {
                sr.append(b).append("\n");
            }

            if (SeqShow) {
                sr.append("\nSeqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence\n");
            } else {
                sr.append("\nSeqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\n");
            }

            bufferedWriter.write(sr.toString());
            int k = 0;
            for (int i = 0; i < bb.size(); i++) {
                int[] z7 = bb.get(i);
                k++;

                String gf = "CRP";
                if (k > 1) {
                    if (refclust != null && k < refclust.length) {
                        if (refclust[k] > 0) {
                            gf = refclust[k] + ":" + refsname[refclust[k] - 1];
                        }
                    }
                }

                for (int j = 0; j < z7.length - 1; j += 2) {

                    for (int w = 0; w < h.length; w++) {
                        if (h[w] > z7[j]) {
                            b = sname[w];
                            break;
                        }
                    }

                    String s0 = "";
                    int x = z7[j] + Math.abs(z7[j + 1]) - 1;
                    if (SeqShow) {
                        if (x > l) {
                            s0 = seq[n].substring(z7[j]);
                        } else {
                            s0 = seq[n].substring(z7[j], x);
                        }
                        if (flanks > 0) {
                            String s1 = "";
                            String s2 = "";
                            if (z7[j] - flanks > 0) {
                                s1 = seq[n].substring(z7[j] - flanks, z7[j]).toUpperCase();
                            } else {
                                if (z7[j] > 1) {
                                    s1 = seq[n].substring(1, z7[1] - 1).toUpperCase();
                                }
                            }
                            if (x + flanks < l) {
                                s2 = seq[n].substring(x, x + flanks).toUpperCase();
                            } else {
                                if (l - x > 0) {
                                    s2 = seq[n].substring(x, l).toUpperCase();
                                }
                            }
                            s0 = s1 + s0 + s2;
                        }
                        if (z7[j + 1] < 0) {
                            s0 = Dna.ComplementDNA2(s0);
                        }
                    }
                    sr = new StringBuilder();
                    String type = (k == 1) ? "STR" : (k == 2) ? "UCRP" : gf;
                    String strand = (z7[j + 1] > 0 || k <= 2) ? "+" : "-";
                    int end = (strand.equals("+")) ? z7[j + 1] : -z7[j + 1];
                    sr.append(b)
                            .append("\t").append(type)
                            .append("\t").append(k)
                            .append("\t").append(z7[j] + 1)
                            .append("\t").append(x + 1)
                            .append("\t").append(end)
                            .append("\t").append(strand)
                            .append("\t").append(s0).append("\n");
                    bufferedWriter.write(sr.toString());
                }
            }
        }
    }
    //  Long-coordinate combined writers (used only by RunCombine).
    //  These mirror the int versions above but address the whole virtual
    //  concatenation with long global coordinates and read sequence through
    //  the SeqStore, so no >2.1 Gb String is ever built. The per-file
    //  reports continue to use the unchanged int writers on local slices.
    // ===================================================================
    void writeTableLong(String reportfile, long l, long[] h, ArrayList<long[]> bbL, SeqStore store) throws IOException {
        String b = (sname != null && sname.length > 0) ? sname[0] : "";
        long duration = (System.nanoTime() - startTime) / 1000000000;

        if (reportfile.length() == 0) {
            reportfile = filePath + "_1.gff";
            if (h.length > 0) {
                reportfile = filePath + ".gff";
            }
        } else {
            reportfile = reportfile + ".gff";
        }

        try (FileWriter fileWriter = new FileWriter(reportfile); BufferedWriter bufferedWriter = new BufferedWriter(fileWriter)) {
            System.out.println("Saving report file: " + reportfile);
            StringBuilder sr = new StringBuilder();
            sr.append("#TotalRepeats (2024-2026) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi) https://github.com/rkalendar/TotalRepeats\n");
            sr.append("#kmer=").append(kmerln).append("\n").append("#Minimal repeat block size=").append(minlenseq).append("\n");
            sr.append("#Sequence length (bp)=").append(l).append("\n");
            sr.append("#Sequence coverage by repeats=").append(String.format("%.2f", repeatslen)).append("%\n");
            sr.append("#Short tandem repeat (STR) sequence coverage=").append(String.format("%.2f", ssrglobal)).append("%\n");
            sr.append("#Sequence gap (bp)=").append((long) gapslen).append(" (").append(String.format("%.4f", gaps)).append("%)\n");
            sr.append("#Masking time taken: ").append(maskduration).append(" seconds\n");
            sr.append("#Total duration: ").append(duration).append(" seconds\n");
            sr.append("#Repeats search for: ");

            if (h.length > 0) {
                for (String filesPath1 : filesPath) {
                    sr.append(filesPath1).append("\n");
                }
            } else {
                sr.append(b).append("\n");
            }

            if (SeqShow) {
                sr.append("\nSeqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence\n");
            } else {
                sr.append("\nSeqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\n");
            }

            bufferedWriter.write(sr.toString());
            int k = 0;
            for (int i = 0; i < bbL.size(); i++) {
                long[] z7 = bbL.get(i);
                k++;

                String gf = "CRP";
                if (k > 1) {
                    if (refclust != null && k < refclust.length) {
                        if (refclust[k] > 0) {
                            gf = refclust[k] + ":" + refsname[refclust[k] - 1];
                        }
                    }
                }

                for (int j = 0; j < z7.length - 1; j += 2) {

                    for (int w = 0; w < h.length; w++) {
                        if (h[w] > z7[j]) {
                            b = sname[w];
                            break;
                        }
                    }

                    String s0 = "";
                    long x = z7[j] + Math.abs(z7[j + 1]) - 1;
                    if (SeqShow) {
                        if (x > l) {
                            s0 = store.substring(z7[j], l);
                        } else {
                            s0 = store.substring(z7[j], x);
                        }
                        if (flanks > 0) {
                            String s1 = "";
                            String s2 = "";
                            if (z7[j] - flanks > 0) {
                                s1 = store.substring(z7[j] - flanks, z7[j]).toUpperCase();
                            } else {
                                if (z7[j] > 1) {
                                    long e = z7[1] - 1;     // faithful to the int version's fixed index
                                    if (e > 1) {
                                        s1 = store.substring(1, e).toUpperCase();
                                    }
                                }
                            }
                            if (x + flanks < l) {
                                s2 = store.substring(x, x + flanks).toUpperCase();
                            } else {
                                if (l - x > 0) {
                                    s2 = store.substring(x, l).toUpperCase();
                                }
                            }
                            s0 = s1 + s0 + s2;
                        }
                        if (z7[j + 1] < 0) {
                            s0 = Dna.ComplementDNA2(s0);
                        }
                    }
                    sr = new StringBuilder();
                    String type = (k == 1) ? "STR" : (k == 2) ? "UCRP" : gf;
                    String strand = (z7[j + 1] > 0 || k <= 2) ? "+" : "-";
                    long end = (strand.equals("+")) ? z7[j + 1] : -z7[j + 1];
                    sr.append(b)
                            .append("\t").append(type)
                            .append("\t").append(k)
                            .append("\t").append(z7[j] + 1)
                            .append("\t").append(x + 1)
                            .append("\t").append(end)
                            .append("\t").append(strand)
                            .append("\t").append(s0).append("\n");
                    bufferedWriter.write(sr.toString());
                }
            }
        }
    }
    /**
     * Pangenomic analysis for the combined run. PangenomeAnalysis is now
     * long-based, so the combined clusters/boundaries are passed straight
     * through in global long coordinates and the report is produced at any size
     * (no int-limit restriction). {@code l} is accepted for call-site
     * uniformity across the combine modes.
     */
    void writePangenome(String reportBase, long[] seqslen, ArrayList<long[]> bbL, long l) throws IOException {
        if (!pangenome || bbL == null || nseq < 2 || seqslen == null || seqslen.length < nseq) {
            return;
        }
        String base = (reportBase == null || reportBase.isEmpty()) ? filePath : reportBase;
        PangenomeAnalysis pa = new PangenomeAnalysis(bbL, seqslen, sname, refclust, refsname);
        pa.write(base);
    }
}
