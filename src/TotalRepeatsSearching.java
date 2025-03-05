
import java.util.Arrays;
import java.awt.BasicStroke;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;

public final class TotalRepeatsSearching {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
        nseq = seq.length;
    }

    public void SetFileName(String a) {
        filePath = a;
    }

    public void SetFileNames(String[] n) {
        filesPath = n;
    }

    public void SetRepeatLen(int kmerln, int minlenseq, int gap) {
        this.kmerln = kmerln;
        if (this.kmerln < 9) {
            this.kmerln = 9;
        }
        this.gap = gap;
        if (this.gap < kmerln) {
            this.gap = kmerln;
        }
        this.minlenseq = minlenseq;
        if (this.minlenseq < kmerln) {
            this.minlenseq = kmerln;
        }
    }

    public void SetFlanks(int i) {
        flanks = i;
    }

    public void SetMaskedPicture(boolean i) {
        MaskedPicture = i;
    }

    public void SetShowSeq(boolean i) {
        SeqShow = i;
    }

    public void SetMasked(boolean i) {
        MaskedShow = i;
    }

    public void SetGFF(boolean i) {
        GFFShow = i;
    }

    public void SetImage(int w, int h) {
        if (w > 0) {
            iwidth = w;
        }
        if (h > 0) {
            iheight = h;
        }
    }

    public void Run2(int k, int nkmer, boolean sensitve) throws IOException {
        startTime = System.nanoTime();
        int[] elemts = new int[nseq];
        int sz = 0;
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            sz = sz + l;
            elemts[i] = sz;
        }

        String s = String.join("n", seq);
        int l = s.length();
        seq = new String[]{s};
        s = "";

        repeatslen = 0;
        gapslen = 0;
        bb = new ArrayList<>();

        LowComplexitySequence m1 = new LowComplexitySequence();
        m1.FindAllSSRs(seq[0], telomers);
        byte[] ssrmsk = m1.MapBytes();
        int[] ssr = m1.IntBlocks();
        repeatslen = m1.GetTotalRepeats();

        bb.add(ssr);
        MaskingSequence ms = new MaskingSequence();

        int[] u = ms.Mask(seq[0], ssrmsk, kmerln, minlenseq, sensitve);
        repeatslen += ms.RepeatLength();
        gapslen = ms.GapsLength();
        repeatslen = (repeatslen * 100) / (l - gapslen);
        gaps = (gapslen * 100) / l;

        System.out.println("Target sequence length = " + l + " nt");
        System.out.println("Sequence coverage by repeats =" + String.format("%.2f", repeatslen) + "%");
        System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)\n");

        if (u.length > 1) {
            if (MaskedShow) {
                MaskSave(0, u, ssr);
            }
            if (MaskedPicture) {
                PictureMasking(u);
            } else {
                System.out.println("Clustering started...");
                ClusteringMasking(nkmer, seq[0], u, sim);
            }
        }

        if (bb != null) {
            if (GFFShow) {
                GffSave(0, l, elemts);
            }
            PictureSave(k, 0, iwidth, iheight, elemts);
        }

    }

    public void Run(int k, int nkmer, boolean sensitve) throws IOException {
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            repeatslen = 0;
            gapslen = 0;
            bb = new ArrayList<>();

            if (l > minlenseq) {
                startTime = System.nanoTime();

                LowComplexitySequence m1 = new LowComplexitySequence();
                m1.FindAllSSRs(seq[i], telomers);
                byte[] ssrmsk = m1.MapBytes();
                int[] ssr = m1.IntBlocks();
                repeatslen = m1.GetTotalRepeats();

                bb.add(ssr);

                MaskingSequence ms = new MaskingSequence();

                int[] u = ms.Mask(seq[i], ssrmsk, kmerln, minlenseq, sensitve);
                repeatslen += ms.RepeatLength();
                gapslen = ms.GapsLength();
                repeatslen = (repeatslen * 100) / (l - gapslen);
                gaps = (gapslen * 100) / l;

                System.out.println("Target sequence length = " + l + " nt");
                System.out.println("Sequence coverage by repeats =" + String.format("%.2f", repeatslen) + "%");
                System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)\n");

                if (u.length > 1) {
                    if (MaskedShow) {
                        MaskSave(i, u, ssr);
                    }
                    if (MaskedPicture) {
                        PictureMasking(u);
                    } else {
                        System.out.println("Clustering started...");
                        ClusteringMasking(nkmer, seq[i], u, sim);
                    }
                }

                if (bb != null) {
                    if (GFFShow) {
                        GffSave(0, l, new int[0]);
                    }
                    PictureSave(k, i, iwidth, iheight, new int[0]);
                }
            }
        }
    }

    private int ClusteringMasking(int nkmer, String seq, int[] z2, int sim) {
        int[][] d;   // d[j][0] = x1; d[j][1] = length;
        int[] q;     // cluster ID for each block
        int ncl;

        SequencesClustering sc = new SequencesClustering(seq, nkmer, z2, sim);
        d = sc.ResultArray();
        q = sc.Result();
        ncl = sc.getNcl();

        if (ncl < 1) {
            return -1;
        }

        for (int i = 0; i < q.length; i++) {
            if (q[i] == 0) {
                q[i] = 2;
            }
        }

        for (int j = 2; j <= ncl; j++) {
            ArrayList<Integer> z = new ArrayList<>();
            z.add(0);
            for (int i = 0; i < q.length; i++) {
                if (q[i] == j) {
                    z.add(d[i][0]);
                    z.add(d[i][1]);
                }
                if (-q[i] == j) {
                    z.add(d[i][0]);
                    z.add(-d[i][1]);
                }
            }

            bb.add(z.stream().mapToInt(Integer::intValue).toArray());
        }
        return bb.size();
    }

    private int PictureMasking(int[] x1) {
        int n = x1.length / 2;
        int[][] d = new int[n][2];
        int k = (n / 50 < 10) ? n : n / 50;

        for (int j = 0; j < n; j++) {
            int p = j * 2;
            d[j][0] = x1[p];
            d[j][1] = x1[p + 1] - x1[p];
        }
        Arrays.sort(d, (int[] a, int[] b) -> {
            return Integer.compare(b[1], a[1]);
        });

        for (int j = 0; j < n; j += k) {
            int[] k7 = new int[k + k + 1];
            k7[0] = k + k;
            int u = 0;
            int i = 0;
            int t = j + k;
            if (t > n) {
                t = n;
            }
            for (i = j; i < t; i++) {
                k7[++u] = d[i][0];
                k7[++u] = d[i][1];
            }
            bb.add(k7);
        }
        return bb.size();
    }

    private void MaskSave(int n, int[] m, int[] ssr) throws IOException {
        String maskedfile = filePath + "_" + (n + 1) + ".msk";
        if (nseq == 1) {
            maskedfile = filePath + ".msk";
        }
        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            System.out.println("Saving masked file: " + maskedfile);

            byte[] c = seq[n].getBytes();
// small to UPPER letter for repeats            
            for (int j = 0; j < m.length; j += 2) {
                for (int i = m[j]; i < m[j] + m[j + 1]; i++) {
                    c[i] = (byte) (c[i] - 32);
                }
            }
//  SSR masking
            for (int j = 1; j < ssr.length; j += 2) {
                for (int i = ssr[j]; i < ssr[j] + ssr[j + 1]; i++) {
                    if (c[i] > 96) {
                        c[i] = (byte) (c[i] - 32);
                    }
                }
            }

            fileWriter.write(">" + sname[n] + " Sequence coverage by repeats = " + String.format("%.2f", repeatslen) + "%\n");
            fileWriter.write(new String(c));
        }
    }

    private void GffSave(int n, int l, int[] h) throws IOException {
        String b = sname[n];
        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Time taken: " + duration + " seconds\n");

        String reportfile = filePath + "_" + (n + 1) + ".gff";
        if (h.length > 0) {
            reportfile = filePath + ".gff";
        }

        try (FileWriter fileWriter = new FileWriter(reportfile); BufferedWriter bufferedWriter = new BufferedWriter(fileWriter)) {
            System.out.println("Saving report file: " + reportfile);
            StringBuilder sr = new StringBuilder();
            sr.append("kmer=").append(kmerln).append("\n").append("Minimal repeat block size=").append(minlenseq).append("\n");
            sr.append("Sequence length (bp)=").append(l).append("\n");
            sr.append("Sequence coverage by repeats=").append(String.format("%.2f", repeatslen)).append("%\n");
            sr.append("Sequence gap (bp)=").append((int) gapslen).append(" (").append(String.format("%.4f", gaps)).append("%)\n");
            sr.append("Time taken: ").append(duration).append(" seconds\n\n");
            sr.append("Repeats search for:\n");

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

                for (int j = 1; j < z7.length; j += 2) {
                    
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
                            s0 = dna.ComplementDNA2(s0);
                        }
                    }
                    sr = new StringBuilder();
                    String type = (k == 1) ? "SSR" : (k == 2) ? "UCRP" : "CRP";
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

    private void PictureSave(int k, int n, int dw, int dh, int[] fastanames) throws IOException {
        int maxClusters = 500;
        int maxImageDimension = 120000;
        int minImageWidth = 4000;
        int minImageHeight = 100;
        int stepPadding = 50;
        int defaultStep = 10;

        // Adjust number of clusters `b`
        int b = Math.min(bb.size(), maxClusters); // Maximum of 1000 clusters

        // Adjust `z` (step between clusters) based on `b`
        int z = calculateClusterStep(b, defaultStep);

        // Calculate width and height
        int l = seq[n].length();
        int width = calculateWidth(k, l, dw, maxImageDimension, minImageWidth);
        int height = calculateHeight(b, z, dh, maxImageDimension, minImageHeight, stepPadding);

        // Calculate dot size
        float dotSize = calculateDotSize(b);

        // Attempt to save the image
        try {
            saveImage(n, l, b, z, width, height, dotSize, fastanames);
        } catch (IOException e) {
            System.out.println("Error saving the picture: " + e.getMessage());
        }
    }

    private int calculateClusterStep(int b, int defaultStep) {
        if (b > 400) {
            return 6;
        }
        if (b > 250) {
            return 7;
        }
        if (b > 100) {
            return 8;
        }
        if (b > 50) {
            return 9;
        }
        return defaultStep;
    }

    private int calculateWidth(int k, int l, int dw, double maxImageDimension, double minImageWidth) {
        double width = k * Math.sqrt(l);
        if (dw > 0) {
            width = dw;
        }
        return (int) Math.max(Math.min(width, maxImageDimension), minImageWidth);
    }

    private int calculateHeight(int b, int z, int dh, int maxImageDimension, int minImageHeight, int stepPadding) {
        int height = b * z + stepPadding; // Adding some padding
        if (dh > 0) {
            height = dh;
        }
        return Math.max(Math.min(height, maxImageDimension), minImageHeight);
    }

    private float calculateDotSize(int b) {
        float dotSize = 12 - (b / 100.0f);
        return Math.max(7.0f, Math.min(dotSize, 7.0f));
    }

    private void saveImage(int n, int l, int b, int z, int width, int height, float dotSize, int[] fastanames) throws IOException {
        double nucleotidesPerPixel = (double) width / l;
        String pngfile = filePath + "_" + (n + 1) + ".png";
        if (nseq == 1) {
            pngfile = filePath + ".png";
        }

        System.out.println("Saving picture " + (width + 100) + "x" + (height + 200) + " : " + pngfile);

        BufferedImage image = new BufferedImage(width + 100, height + 200, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        g2d.setStroke(new BasicStroke(dotSize));
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, width + 100, height + 200);
        g2d.setColor(Color.BLACK);
        drawLinesAndLabels(g2d, l, width);

        if (fastanames.length > 0) {
            int x1 = 0;
            for (int i = 0; i < fastanames.length; i++) {
                x1 = (int) (x1 * nucleotidesPerPixel);
                g2d.drawLine(x1 + 50, 1, x1 + 50, 20);//(x1, y, x2, y)
                g2d.drawString(sname[i], x1 + 55, 25);
                x1 = fastanames[i];
            }
        } else {
            g2d.drawLine(50, 1, 50, 20);
            g2d.drawString(sname[n], 55, 25);
        }

        drawClusters(g2d, b, z, nucleotidesPerPixel);
        g2d.dispose();
        File outputFile = new File(pngfile);
        ImageIO.write(image, "png", outputFile);
    }

    private void drawLinesAndLabels(Graphics2D g2d, int l, int width) {
        g2d.drawLine(50, 55, width + 50, 55); // top line (x1, y, x2, y)
        int f = 10;
        int w = width / f;
        int d = l / f;
        for (int i = 0; i <= f; i++) {
            g2d.drawLine(i * w + 50, 45, i * w + 50, 55);
            int v = 1 + i * d;
            if (v > l) {
                v = l;
            }
            g2d.drawString(String.format("%,d", v), 55 + i * w, 50);
        }
    }

    private void drawClusters(Graphics2D g2d, int b, int z, double w1) {
        // Color DarkRed = new Color(153, 0, 0); //https://teaching.csse.uwa.edu.au/units/CITS1001/colorinfo.html
        Color DarkGreen = new Color(0, 102, 0);
        Color Brown = new Color(102, 51, 0);
        for (int i = 0; i < b; i++) {
            int[] z7 = bb.get(i);

            // Gray lines at height 22
            for (int j = 1; j < z7.length; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = (z7[j + 1] > 0) ? 50 + (int) ((z7[j] + z7[j + 1]) * w1) : 50 + (int) ((z7[j] - z7[j + 1]) * w1);
                g2d.setColor(Brown);
                g2d.drawLine(x1, 63, x2, 63); // draw dark gray line (x1, y, x2, y)
            }

            //  int y = (i > 1) ? 120 + (i * z) : 80 + (i * 20);
            int y = (i > 10) ? 190 + (i * z) : 80 + (i * 20);

            for (int j = 1; j < z7.length; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 50;
                if (z7[j + 1] > 0) {
                    x2 = x2 + (int) ((z7[j] + z7[j + 1]) * w1);
                    if (i == 0) {
                        g2d.setColor(DarkGreen);
                    } else {
                        g2d.setColor(Color.BLUE);
                    }
                } else {
                    x2 = x2 + (int) ((z7[j] - z7[j + 1]) * w1);
                    g2d.setColor(Color.RED);
                }
                g2d.drawLine(x1, y, x2, y); // draw blue line
            }

        }
    }

    private double gaps = 0;
    private double gapslen = 0;
    private double repeatslen = 0;
    private long startTime;
    private int nseq = 0;
    private int iwidth = 0;
    private int iheight = 0;
    private int minlenseq = 50;     // Minimal repeat block size
    private int kmerln = 21;       // kmer=12-21
    private int flanks = 20;
    private int gap = 21;         // gap between repeat blocks, gap=kmer
    private final int sim = 80;
    private final int telomers = 14; // Kmax=11 ->SSR  Kmax=14 -> telomers
    private boolean SeqShow;
    private boolean MaskedShow;
    private boolean MaskedPicture;
    private boolean GFFShow;
    private String filePath;
    private String[] filesPath;
    private String[] seq;
    private String[] sname;
    private ArrayList<int[]> bb;
}
