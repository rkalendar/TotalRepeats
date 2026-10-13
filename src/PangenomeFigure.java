import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Draws the pangenome summary figure of a combined run as a single SVG page
 * with four panels:
 *
 * <ol>
 *   <li><b>Rarefaction</b> — expected pangenome and core size against the
 *       number of genomes sampled, with the Heaps' law exponent;</li>
 *   <li><b>Family frequency spectrum</b> — how many families occur in exactly
 *       k genomes (the pangenome U-curve);</li>
 *   <li><b>Per-genome composition</b> — core / accessory / unique families in
 *       every genome, in the order of the genome tree;</li>
 *   <li><b>Presence/absence map</b> — families × genomes, with the UPGMA tree
 *       of the genomes above the columns. Families are sorted by frequency and
 *       pattern; when there are more families than drawable rows, consecutive
 *       families are binned and the cell opacity shows the fraction present.</li>
 * </ol>
 *
 * <p>Everything is computed by {@link PangenomeAnalysis}; this class only
 * draws. The output is plain SVG with no external resources, so it opens in any
 * browser or vector editor, and all numbers are written with {@link Locale#ROOT}.
 */
final class PangenomeFigure {

    // Okabe–Ito colours: distinguishable under the common forms of colour blindness.
    static final String CORE = "#009E73";
    static final String ACCESSORY = "#0072B2";
    static final String UNIQUE = "#E69F00";
    private static final String PAN = "#222222";
    private static final String GRID = "#DDDDDD";
    private static final String ABSENT = "#F0F0F0";

    private static final int W = 1400;
    private static final int H = 880;
    private static final int GUTTER = 40;
    private static final int ROW1_Y = 80, ROW1_H = 320;
    private static final int ROW2_Y = 455, ROW2_H = 385;
    private static final int PANEL_W = (W - 3 * GUTTER) / 2;

    /** Upper bound on heat-map cells, which bounds the file size. */
    private static final int MAX_CELLS = 100_000;
    private static final int MIN_BINS = 200;

    private final PangenomeAnalysis pa;

    PangenomeFigure(PangenomeAnalysis pa) {
        this.pa = pa;
    }

    void write(String file) throws IOException {
        String svg = render();
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome figure " + W + "x" + H + " : " + file);
            w.write(svg);
        }
    }

    String render() {
        StringBuilder sb = new StringBuilder(1 << 16);
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        sb.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"").append(W).append("\" height=\"").append(H)
                .append("\" viewBox=\"0 0 ").append(W).append(' ').append(H).append("\">\n");
        sb.append("  <rect x=\"0\" y=\"0\" width=\"").append(W).append("\" height=\"").append(H).append("\" fill=\"#FFFFFF\"/>\n");
        sb.append("  <style><![CDATA[\n")
                .append("    text { font-family: Arial, Helvetica, sans-serif; fill:#000; }\n")
                .append("    .title { font-size:18px; font-weight:bold; }\n")
                .append("    .panel { font-size:14px; font-weight:bold; }\n")
                .append("    .lab { font-size:11px; }\n")
                .append("    .small { font-size:9px; }\n")
                .append("    .axis { stroke:#000; stroke-width:1; fill:none; }\n")
                .append("    .grid { stroke:").append(GRID).append("; stroke-width:1; }\n")
                .append("    .tree { stroke:#000; stroke-width:1; fill:none; }\n")
                .append("  ]]></style>\n");

        int N = pa.genomeCount();
        int F = pa.families().size();
        sb.append("  <text class=\"title\" x=\"").append(GUTTER).append("\" y=\"30\">TotalRepeats pangenome: ")
                .append(N).append(" genomes, ").append(F).append(" repeat families (core ").append(pa.coreCount())
                .append(", accessory ").append(pa.accessoryCount()).append(", unique ").append(pa.uniqueCount())
                .append(")</text>\n");
        legend(sb, GUTTER, 48);

        int x1 = GUTTER, x2 = 2 * GUTTER + PANEL_W;
        rarefaction(sb, x1, ROW1_Y, PANEL_W, ROW1_H);
        spectrum(sb, x2, ROW1_Y, PANEL_W, ROW1_H);
        perGenome(sb, x1, ROW2_Y, PANEL_W, ROW2_H);
        heatmap(sb, x2, ROW2_Y, PANEL_W, ROW2_H);

        sb.append("</svg>\n");
        return sb.toString();
    }

    // ── Panels ──────────────────────────────────────────────────────────────
    private void legend(StringBuilder sb, int x, int y) {
        String[] lab = {"core (all genomes)", "accessory (some)", "unique (one genome)"};
        String[] col = {CORE, ACCESSORY, UNIQUE};
        int cx = x;
        for (int i = 0; i < 3; i++) {
            rect(sb, cx, y - 10, 12, 12, col[i], 1.0);
            text(sb, "lab", cx + 16, y, lab[i], null);
            cx += 16 + 7 * lab[i].length() + 24;
        }
        String note = "Heaps' law: alpha <= 1 open, > 1 closed pangenome";
        text(sb, "lab", W - GUTTER, y, note, "end");
    }

    /** Panel 1: expected pangenome / core size against genomes sampled. */
    private void rarefaction(StringBuilder sb, int x, int y, int w, int h) {
        text(sb, "panel", x, y - 8, "Rarefaction: expected pan-genome and core size", null);
        int N = pa.genomeCount();
        double[] pan = pa.expectedPan();
        double[] core = pa.expectedCore();

        int px0 = x + 60, px1 = x + w - 20, py0 = y + 15, py1 = y + h - 45;
        double ymax = 0;
        for (int n = 1; n <= N; n++) {
            ymax = Math.max(ymax, pan[n]);
        }
        double step = niceStep(ymax, 5);
        ymax = Math.max(step, Math.ceil(ymax / step) * step);
        yAxis(sb, px0, px1, py0, py1, ymax, step, "families");

        // x axis: n = 1..N
        line(sb, "axis", px0, py1, px1, py1);
        int every = Math.max(1, (int) Math.ceil(N / 20.0));
        for (int n = 1; n <= N; n++) {
            double xx = xOf(n, N, px0, px1);
            if ((n - 1) % every == 0 || n == N) {
                line(sb, "axis", xx, py1, xx, py1 + 4);
                text(sb, "lab", xx, py1 + 16, Integer.toString(n), "middle");
            }
        }
        text(sb, "lab", (px0 + px1) / 2.0, py1 + 32, "genomes sampled", "middle");

        polyline(sb, pan, N, px0, px1, py0, py1, ymax, PAN, 2.0);
        polyline(sb, core, N, px0, px1, py0, py1, ymax, CORE, 2.0);

        // in-plot legend and the Heaps' law estimate
        double lx = px0 + 12, ly = py0 + 14;
        line(sb, "", lx, ly, lx + 24, ly, PAN, 2.0);
        text(sb, "lab", lx + 30, ly + 4, "pan-genome (families seen)", null);
        line(sb, "", lx, ly + 18, lx + 24, ly + 18, CORE, 2.0);
        text(sb, "lab", lx + 30, ly + 22, "core (families in every sampled genome)", null);

        String heaps;
        double a = pa.heapsAlpha();
        if (Double.isNaN(a)) {
            heaps = "Heaps' law alpha: n/a (needs >= 3 genomes)";
        } else {
            heaps = "Heaps' law alpha = " + PangenomeAnalysis.f2(a) + " (" + (a > 1 ? "closed" : "open") + ")";
            if (!Double.isNaN(pa.heapsGamma())) {
                heaps += ", gamma = " + PangenomeAnalysis.f2(pa.heapsGamma());
            }
        }
        text(sb, "lab", px1 - 6, py1 - 8, heaps, "end");
    }

    /** Panel 2: number of families present in exactly k genomes. */
    private void spectrum(StringBuilder sb, int x, int y, int w, int h) {
        text(sb, "panel", x, y - 8, "Family frequency spectrum", null);
        int N = pa.genomeCount();
        int[] spec = pa.spectrum();

        int px0 = x + 60, px1 = x + w - 20, py0 = y + 15, py1 = y + h - 45;
        double ymax = 0;
        for (int k = 1; k <= N; k++) {
            ymax = Math.max(ymax, spec[k]);
        }
        double step = niceStep(ymax, 5);
        ymax = Math.max(step, Math.ceil(ymax / step) * step);
        yAxis(sb, px0, px1, py0, py1, ymax, step, "families");
        line(sb, "axis", px0, py1, px1, py1);

        double slot = (double) (px1 - px0) / N;
        double bw = Math.max(1.0, slot * 0.8);
        int every = Math.max(1, (int) Math.ceil(N / 25.0));
        for (int k = 1; k <= N; k++) {
            double bx = px0 + (k - 1) * slot + (slot - bw) / 2;
            double bh = (py1 - py0) * spec[k] / ymax;
            String col = (k == N) ? CORE : (k == 1 && N > 1) ? UNIQUE : ACCESSORY;
            if (spec[k] > 0) {
                rect(sb, bx, py1 - bh, bw, bh, col, 1.0);
                if (slot >= 18) {
                    text(sb, "small", bx + bw / 2, py1 - bh - 3, Integer.toString(spec[k]), "middle");
                }
            }
            if ((k - 1) % every == 0 || k == N) {
                text(sb, "lab", bx + bw / 2, py1 + 16, Integer.toString(k), "middle");
            }
        }
        text(sb, "lab", (px0 + px1) / 2.0, py1 + 32, "present in k genomes", "middle");
    }

    /** Panel 3: core / accessory / unique family counts per genome, in tree order. */
    private void perGenome(StringBuilder sb, int x, int y, int w, int h) {
        text(sb, "panel", x, y - 8, "Families per genome (genomes in tree order)", null);
        int N = pa.genomeCount();
        int[] order = pa.leafOrder();
        int[][] cc = pa.genomeCategoryCounts();
        String[] names = pa.genomeNames();

        int labelH = 95;
        int px0 = x + 60, px1 = x + w - 20, py0 = y + 15, py1 = y + h - labelH;
        double ymax = 0;
        for (int g = 0; g < N; g++) {
            ymax = Math.max(ymax, cc[g][0]);
        }
        double step = niceStep(ymax, 5);
        ymax = Math.max(step, Math.ceil(ymax / step) * step);
        yAxis(sb, px0, px1, py0, py1, ymax, step, "families");
        line(sb, "axis", px0, py1, px1, py1);

        double slot = (double) (px1 - px0) / N;
        double bw = Math.max(1.0, slot * 0.8);
        double scale = (py1 - py0) / ymax;
        for (int c = 0; c < N; c++) {
            int g = order[c];
            double bx = px0 + c * slot + (slot - bw) / 2;
            double top = py1;
            String[] col = {CORE, ACCESSORY, UNIQUE};
            for (int part = 1; part <= 3; part++) {
                double bh = cc[g][part] * scale;
                if (bh > 0) {
                    rect(sb, bx, top - bh, bw, bh, col[part - 1], 1.0);
                    top -= bh;
                }
            }
            columnLabel(sb, bx + bw / 2, py1 + 6, slot, names[g]);
        }
    }

    /** Panel 4: UPGMA tree over a families × genomes presence map. */
    private void heatmap(StringBuilder sb, int x, int y, int w, int h) {
        text(sb, "panel", x, y - 8, "Presence/absence of families (rows) in genomes (columns), UPGMA on Jaccard distance", null);
        int N = pa.genomeCount();
        int[] order = pa.leafOrder();
        String[] names = pa.genomeNames();
        List<PangenomeAnalysis.Family> fams = new ArrayList<>(pa.families());
        int F = fams.size();

        int treeH = 80, labelH = 95;
        int px0 = x + 30, px1 = x + w - 15;
        int treeTop = y + 10, treeBottom = treeTop + treeH;
        int hy0 = treeBottom + 6, hy1 = y + h - labelH;
        double cw = (double) (px1 - px0) / N;

        // column centre per tree position
        double[] leafX = new double[N];
        int[] colOf = new int[N];
        for (int c = 0; c < N; c++) {
            leafX[order[c]] = px0 + (c + 0.5) * cw;
            colOf[order[c]] = c;
        }
        drawTree(sb, pa.tree(), leafX, treeTop, treeBottom);

        // background for absent cells, and the category band on the left
        rect(sb, px0, hy0, px1 - px0, hy1 - hy0, ABSENT, 1.0);
        if (F == 0) {
            text(sb, "lab", (px0 + px1) / 2.0, (hy0 + hy1) / 2.0, "no repeat families", "middle");
        } else {
            // Sort rows: frequency desc, then the presence pattern in column order
            // (so shared blocks line up), then size.
            final String[] patInTreeOrder = new String[F];
            for (int i = 0; i < F; i++) {
                PangenomeAnalysis.Family f = fams.get(i);
                char[] p = new char[N];
                for (int g = 0; g < N; g++) {
                    p[colOf[g]] = f.count[g] > 0 ? '1' : '0';
                }
                patInTreeOrder[i] = new String(p);
            }
            Integer[] idx = new Integer[F];
            for (int i = 0; i < F; i++) {
                idx[i] = i;
            }
            java.util.Arrays.sort(idx, Comparator
                    .comparingInt((Integer i) -> -fams.get(i).frequency)
                    .thenComparing(i -> patInTreeOrder[i], Comparator.reverseOrder())
                    .thenComparingLong(i -> -fams.get(i).totalBp)
                    .thenComparingInt(i -> fams.get(i).clusterId));

            int targetBins = Math.max(MIN_BINS, MAX_CELLS / Math.max(1, N));
            int rowsPerBin = Math.max(1, (int) Math.ceil((double) F / targetBins));
            int bins = (F + rowsPerBin - 1) / rowsPerBin;
            double rh = (double) (hy1 - hy0) / bins;

            // category band (left of the map): majority category of each bin
            String prevCol = null;
            double runStart = 0;
            for (int b = 0; b <= bins; b++) {
                String col = (b < bins) ? binCategoryColour(fams, idx, b * rowsPerBin, Math.min(F, (b + 1) * rowsPerBin)) : null;
                if (b == bins || !col.equals(prevCol)) {
                    if (prevCol != null) {
                        rect(sb, px0 - 10, hy0 + runStart * rh, 7, (b - runStart) * rh, prevCol, 1.0);
                    }
                    prevCol = col;
                    runStart = b;
                }
            }

            // cells: per column, runs of identical (colour, opacity) bins are merged into one rect
            for (int c = 0; c < N; c++) {
                int g = order[c];
                double cx = px0 + c * cw;
                String prev = null;
                int prevB = 0;
                for (int b = 0; b <= bins; b++) {
                    String key = null;
                    if (b < bins) {
                        int from = b * rowsPerBin, to = Math.min(F, (b + 1) * rowsPerBin);
                        int present = 0, core = 0, acc = 0, uniq = 0;
                        for (int r = from; r < to; r++) {
                            PangenomeAnalysis.Family f = fams.get(idx[r]);
                            if (f.count[g] > 0) {
                                present++;
                                switch (f.category) {
                                    case "core" -> core++;
                                    case "accessory" -> acc++;
                                    default -> uniq++;
                                }
                            }
                        }
                        if (present > 0) {
                            String col = (core >= acc && core >= uniq) ? CORE : (acc >= uniq) ? ACCESSORY : UNIQUE;
                            double op = Math.round(100.0 * present / (to - from)) / 100.0;
                            key = col + "|" + op;
                        }
                    }
                    if (b == bins || (key == null ? prev != null : !key.equals(prev))) {
                        if (prev != null) {
                            int bar = prev.indexOf('|');
                            rect(sb, cx, hy0 + prevB * rh, cw, (b - prevB) * rh, prev.substring(0, bar),
                                    Double.parseDouble(prev.substring(bar + 1)));
                        }
                        prev = key;
                        prevB = b;
                    }
                }
            }
            String rows = (rowsPerBin == 1) ? F + " families, one row each"
                    : F + " families in " + bins + " rows of " + rowsPerBin + "; cell opacity = fraction present";
            text(sb, "small", px1, hy1 + labelH - 6, rows + "; left band = category", "end");
        }

        for (int c = 0; c < N; c++) {
            columnLabel(sb, px0 + (c + 0.5) * cw, hy1 + 6, cw, names[order[c]]);
        }
    }

    private static String binCategoryColour(List<PangenomeAnalysis.Family> fams, Integer[] idx, int from, int to) {
        int core = 0, acc = 0, uniq = 0;
        for (int r = from; r < to; r++) {
            switch (fams.get(idx[r]).category) {
                case "core" -> core++;
                case "accessory" -> acc++;
                default -> uniq++;
            }
        }
        return (core >= acc && core >= uniq) ? CORE : (acc >= uniq) ? ACCESSORY : UNIQUE;
    }

    /** Draws the UPGMA tree with the leaves at {@code bottom} and the root at {@code top}. */
    private void drawTree(StringBuilder sb, PangenomeAnalysis.Node root, double[] leafX, int top, int bottom) {
        if (root == null) {
            return;
        }
        double maxH = root.height > 0 ? root.height : 1.0;
        treeNode(sb, root, leafX, top, bottom, maxH);
    }

    private double treeNode(StringBuilder sb, PangenomeAnalysis.Node t, double[] leafX, int top, int bottom, double maxH) {
        if (t.leaf >= 0) {
            return leafX[t.leaf];
        }
        double lx = treeNode(sb, t.left, leafX, top, bottom, maxH);
        double rx = treeNode(sb, t.right, leafX, top, bottom, maxH);
        double yy = bottom - (bottom - top) * (t.height / maxH);
        double ly = bottom - (bottom - top) * (t.left.height / maxH);
        double ry = bottom - (bottom - top) * (t.right.height / maxH);
        sb.append("  <path class=\"tree\" d=\"M").append(c(lx)).append(' ').append(c(ly))
                .append(" V").append(c(yy)).append(" H").append(c(rx)).append(" V").append(c(ry)).append("\"/>\n");
        return (lx + rx) / 2.0;
    }

    // ── Drawing helpers ─────────────────────────────────────────────────────
    private static double xOf(int n, int N, int px0, int px1) {
        return (N <= 1) ? (px0 + px1) / 2.0 : px0 + (double) (n - 1) / (N - 1) * (px1 - px0);
    }

    private static void polyline(StringBuilder sb, double[] v, int N, int px0, int px1, int py0, int py1,
            double ymax, String colour, double width) {
        StringBuilder pts = new StringBuilder();
        for (int n = 1; n <= N; n++) {
            double xx = xOf(n, N, px0, px1);
            double yy = py1 - (py1 - py0) * (v[n] / ymax);
            pts.append(c(xx)).append(',').append(c(yy)).append(' ');
        }
        sb.append("  <polyline points=\"").append(pts.toString().trim()).append("\" fill=\"none\" stroke=\"")
                .append(colour).append("\" stroke-width=\"").append(c(width)).append("\"/>\n");
        for (int n = 1; n <= N; n++) {
            double xx = xOf(n, N, px0, px1);
            double yy = py1 - (py1 - py0) * (v[n] / ymax);
            sb.append("  <circle cx=\"").append(c(xx)).append("\" cy=\"").append(c(yy))
                    .append("\" r=\"3\" fill=\"").append(colour).append("\"/>\n");
        }
    }

    private static void yAxis(StringBuilder sb, int px0, int px1, int py0, int py1, double ymax, double step, String label) {
        line(sb, "axis", px0, py0, px0, py1);
        for (double v = 0; v <= ymax + 1e-9; v += step) {
            double yy = py1 - (py1 - py0) * (v / ymax);
            if (v > 0) {
                line(sb, "grid", px0 + 1, yy, px1, yy);
            }
            line(sb, "axis", px0 - 4, yy, px0, yy);
            text(sb, "lab", px0 - 7, yy + 4, num(v), "end");
        }
        sb.append("  <text class=\"lab\" transform=\"translate(").append(px0 - 46).append(' ')
                .append((py0 + py1) / 2).append(") rotate(-90)\" text-anchor=\"middle\">").append(esc(label)).append("</text>\n");
    }

    /** Rotated column label, or none when the columns are too narrow to read. */
    private static void columnLabel(StringBuilder sb, double cx, double top, double slot, String name) {
        if (slot < 6) {
            return;
        }
        double size = Math.min(11, Math.max(6, slot * 0.9));
        String s = (name.length() > 28) ? name.substring(0, 27) + "…" : name;
        sb.append("  <text style=\"font-size:").append(c(size)).append("px\" text-anchor=\"end\" transform=\"translate(")
                .append(c(cx + size / 3)).append(' ').append(c(top + 4)).append(") rotate(-60)\">").append(esc(s)).append("</text>\n");
    }

    static double niceStep(double range, int targetTicks) {
        if (range <= 0) {
            return 1;
        }
        double raw = range / targetTicks;
        double mag = Math.pow(10, Math.floor(Math.log10(raw)));
        double r = raw / mag;
        double step = (r < 1.5) ? 1 : (r < 3.5) ? 2 : (r < 7.5) ? 5 : 10;
        return step * mag;
    }

    private static void rect(StringBuilder sb, double x, double y, double w, double h, String fill, double opacity) {
        sb.append("  <rect x=\"").append(c(x)).append("\" y=\"").append(c(y)).append("\" width=\"").append(c(w))
                .append("\" height=\"").append(c(h)).append("\" fill=\"").append(fill).append('"');
        if (opacity < 1.0) {
            sb.append(" fill-opacity=\"").append(c(opacity)).append('"');
        }
        sb.append("/>\n");
    }

    private static void line(StringBuilder sb, String cls, double x1, double y1, double x2, double y2) {
        sb.append("  <line class=\"").append(cls).append("\" x1=\"").append(c(x1)).append("\" y1=\"").append(c(y1))
                .append("\" x2=\"").append(c(x2)).append("\" y2=\"").append(c(y2)).append("\"/>\n");
    }

    private static void line(StringBuilder sb, String cls, double x1, double y1, double x2, double y2, String stroke, double width) {
        sb.append("  <line x1=\"").append(c(x1)).append("\" y1=\"").append(c(y1))
                .append("\" x2=\"").append(c(x2)).append("\" y2=\"").append(c(y2)).append("\" stroke=\"").append(stroke)
                .append("\" stroke-width=\"").append(c(width)).append("\"/>\n");
    }

    private static void text(StringBuilder sb, String cls, double x, double y, String s, String anchor) {
        sb.append("  <text class=\"").append(cls).append("\" x=\"").append(c(x)).append("\" y=\"").append(c(y)).append('"');
        if (anchor != null) {
            sb.append(" text-anchor=\"").append(anchor).append('"');
        }
        sb.append('>').append(esc(s)).append("</text>\n");
    }

    /** Coordinate formatting: one decimal, point separator, no trailing ".0". */
    private static String c(double v) {
        if (v == Math.rint(v)) {
            return Long.toString((long) v);
        }
        return String.format(Locale.ROOT, "%.2f", v);
    }

    private static String num(double v) {
        if (v == Math.rint(v)) {
            return Long.toString((long) v);
        }
        return String.format(Locale.ROOT, "%.1f", v);
    }

    private static String esc(String s) {
        if (s == null) {
            return "";
        }
        StringBuilder out = new StringBuilder(s.length() + 8);
        for (int i = 0; i < s.length(); i++) {
            char ch = s.charAt(i);
            switch (ch) {
                case '&' -> out.append("&amp;");
                case '<' -> out.append("&lt;");
                case '>' -> out.append("&gt;");
                case '"' -> out.append("&quot;");
                case '\'' -> out.append("&apos;");
                default -> out.append(ch);
            }
        }
        return out.toString();
    }
}
