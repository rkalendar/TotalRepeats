import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageTypeSpecifier;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.metadata.IIOMetadataNode;
import javax.imageio.stream.ImageOutputStream;

/**
 * Draws the repeat-landscape figures — the scalable SVG and the raster PNG —
 * for one analysed sequence or for a whole combined run.
 *
 * <p>Extracted from {@code TotalRepeatsSearching}, which previously held the
 * orchestration, the report writers and the two renderers in a single class and
 * passed data between them through mutable instance fields. The renderers read
 * only four of those fields, so they are taken here as constructor arguments and
 * the drawing code itself is unchanged.
 *
 * <p>Both output formats derive every coordinate from the same helpers —
 * {@link #clusterRowY} for the vertical position of a cluster row, and the
 * layout constants below for the canvas — so an SVG and a PNG produced from one
 * run cannot disagree about where a row sits. They did disagree before these
 * values were shared.
 *
 * <p>Instances are cheap and short-lived: the caller builds one per figure, so
 * the cluster list is read at the moment of drawing rather than cached.
 */
final class RepeatFigure {

    // ── Canvas limits, formerly repeated in each of the three entry points ──
    private static final int MAX_CLUSTERS = 500;
    private static final int MAX_IMAGE_DIMENSION = 120000;
    private static final int MIN_IMAGE_WIDTH = 4000;
    private static final int MIN_IMAGE_HEIGHT = 100;
    private static final int STEP_PADDING = 20;

    /** Repeat blocks per cluster, in the coordinates of the sequence being drawn. */
    private final ArrayList<int[]> bb;
    /** Sequence names, used for the ruler labels. */
    private final String[] sname;
    /** Base path of the analysed file; the extension is appended per format. */
    private final String filePath;
    /** Number of sequences in the run, which selects the output naming. */
    private final int nseq;

    RepeatFigure(ArrayList<int[]> bb, String[] sname, String filePath, int nseq) {
        this.bb = bb;
        this.sname = sname;
        this.filePath = filePath;
        this.nseq = nseq;
    }

    // ── Entry points ────────────────────────────────────────────────────────

    /** Writes the SVG figure for a single sequence. */
    void writeSvg(String reportfile, int k, int n, int len, int dw, int dh, int[] seqslen)
            throws IOException {
        int b = Math.min(bb.size(), MAX_CLUSTERS);
        int z = calculateClusterStep(b);
        float dotSize = calculateDotSize(b);
        int width = calculateWidth(k, len, dw, MAX_IMAGE_DIMENSION, MIN_IMAGE_WIDTH);
        int height = calculateHeight(b, z, dh, MAX_IMAGE_DIMENSION, MIN_IMAGE_HEIGHT, STEP_PADDING);

        SaveSVG(reportfile, k, n, len, b, z, width, height, dotSize, seqslen);
    }

    /**
     * Writes the SVG figure for a combined run, addressed in global long
     * coordinates so that a concatenation above ~2.1 Gb is still drawable.
     */
    void writeSvgLong(String reportfile, int k, long len, int dw, int dh,
            long[] seqslen, ArrayList<long[]> bbL) throws IOException {
        int b = Math.min(bbL.size(), MAX_CLUSTERS);
        int z = calculateClusterStep(b);
        float dotSize = calculateDotSize(b);
        // width is computed inline in long to avoid the int overflow that
        // calculateWidth(int l, ...) would hit for a >2.1 Gb concatenation.
        double rawWidth = (dw > 0) ? dw : k * Math.sqrt((double) len);
        int width = (int) Math.max(Math.min(rawWidth, MAX_IMAGE_DIMENSION), MIN_IMAGE_WIDTH);
        int height = calculateHeight(b, z, dh, MAX_IMAGE_DIMENSION, MIN_IMAGE_HEIGHT, STEP_PADDING);

        SaveSVGLong(reportfile, k, len, b, z, width, height, dotSize, seqslen, bbL);
    }

    /**
     * Writes the PNG figure. A failure to write the image is reported but not
     * propagated: the annotation and mask outputs are already on disk by then,
     * and losing the raster copy should not fail the run.
     */
    void writePng(String reportfile, int k, int n, int len, int dw, int dh, int[] seqslen) {
        int b = Math.min(bb.size(), MAX_CLUSTERS);
        int z = calculateClusterStep(b);
        float dotSize = calculateDotSize(b);
        int width = calculateWidth(k, len, dw, MAX_IMAGE_DIMENSION, MIN_IMAGE_WIDTH);
        int height = calculateHeight(b, z, dh, MAX_IMAGE_DIMENSION, MIN_IMAGE_HEIGHT, STEP_PADDING);

        try {
            SaveImage(reportfile, k, n, len, b, z, width, height, dotSize, seqslen);
        } catch (IOException e) {
            // Non-fatal: the GFF/mask outputs are already written; report the
            // image failure clearly on stderr rather than as a stdout success line.
            System.err.println("ERROR: failed to save image " + reportfile + " — " + e.getMessage());
        }
    }

    // ── Shared helpers ──────────────────────────────────────────────────────

    private static String esc(String s) {
        if (s == null) {
            return "";
        }
        StringBuilder out = new StringBuilder((int) (s.length() * 1.1));
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            switch (c) {
                case '&' ->
                    out.append("&amp;");
                case '<' ->
                    out.append("&lt;");
                case '>' ->
                    out.append("&gt;");
                case '"' ->
                    out.append("&quot;");
                case '\'' ->
                    out.append("&apos;");
                default ->
                    out.append(c);
            }
        }
        return out.toString();
    }

// Format integers with thousands separators like your PNG labels
    private static String formatThousands(int v) {
        return String.format("%,d", v);
    }
    private static String formatThousandsLong(long v) {
        return String.format("%,d", v);
    }
    private int calculateClusterStep(int b) {
        if (b > 500) {
            return 10;
        }
        if (b > 400) {
            return 11;
        }
        if (b > 300) {
            return 12;
        }
        if (b > 200) {
            return 13;
        }
        if (b > 100) {
            return 14;
        }
        return 16;
    }

    private int calculateWidth(int k, int l, int dw, double maxImageDimension, double minImageWidth) {
        double width = k * Math.sqrt(l);
        if (dw > 0) {
            width = dw;
        }
        return (int) Math.max(Math.min(width, maxImageDimension), minImageWidth);
    }

    private int calculateHeight(int b, int z, int dh, int maxImageDimension, int minImageHeight, int stepPadding) {
        int height = b * z + stepPadding;
        if (dh > 0) {
            height = dh;
        }
        return Math.max(Math.min(height, maxImageDimension), minImageHeight);
    }

    private float calculateDotSize(int b) {
        float dotSize = 20 - (b / 100.0f);
        return Math.max(12.0f, dotSize);
    }

    // ── SVG ─────────────────────────────────────────────────────────────────

    private void SaveSVG(String svgfile, int k, int n, int l, int b, int z, int width, int height, float dotSize, int[] seqslen) throws IOException {
        final double nucleotidesPerPixel = (double) width / l;
        if (svgfile.length() == 0) {
            svgfile = filePath + "_" + (n + 1) + ".svg";
            if (nseq == 1) {
                svgfile = filePath + ".svg";
            }
        } else {
            svgfile = svgfile + ".svg";
        }

        System.out.println("Saving SVG " + (width + 100) + "x" + (height + 200) + " : " + svgfile);
        StringBuilder sb = new StringBuilder(1 << 20); // pre-allocate
        // SVG header
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        sb.append("<svg xmlns=\"http://www.w3.org/2000/svg\" ")
                .append("xmlns:xlink=\"http://www.w3.org/1999/xlink\" ")
                .append("width=\"").append(width + 100).append("\" ")
                .append("height=\"").append(height + 200).append("\" ")
                .append("viewBox=\"0 0 ").append(width + 100).append(" ").append(height + 200).append("\">\n");

        // Background
        sb.append("  <rect x=\"0\" y=\"0\" width=\"").append(width + 100).append("\" height=\"").append(height + 200).append("\" fill=\"#FFFFFF\"/>\n");

        // Styles (adjustable)
        sb.append("  <style><![CDATA[\n")
                .append("    .axis { stroke:#000; stroke-width:1; }\n")
                .append("    .tick { stroke:#000; stroke-width:1; }\n")
                .append("    .labelBig { font-family:monospace; font-size:18px; font-weight:bold; fill:#000; }\n")
                .append("    .labelSmall { font-family:monospace; font-size:8px; fill:#000; }\n")
                .append("    .labelMed { font-family:monospace; font-size:16px; fill:#000; }\n")
                .append("    .brown { stroke:#663300; }\n") // Brown
                .append("    .blue { stroke:#0000FF; }\n")
                .append("    .red { stroke:#FF0000; }\n")
                .append("    .darkgreen { stroke:#006600; }\n")
                .append("  ]]></style>\n");

        // Top ruler line
        sb.append("  <line class=\"axis\" x1=\"50\" y1=\"55\" x2=\"").append(50 + width).append("\" y2=\"55\"/>\n");

        // Ruler ticks + numeric labels (emulating drawLinesAndLabels)
        int f = k + 5;
        int w = width / f;
        int d = l / f;

        for (int i = 0; i <= f; i++) {
            int xTick = 50 + i * w;
            sb.append("  <line class=\"tick\" x1=\"").append(xTick).append("\" y1=\"45\" x2=\"").append(xTick).append("\" y2=\"55\"/>\n");
            int v = 1 + i * d;
            if (v > l) {
                v = l;
            }
            // Map cumulative coordinate into per-sequence coordinate if multiple seqs
            if (seqslen != null && seqslen.length > 0) {
                for (int j = 1; j < seqslen.length; j++) {
                    if (v >= seqslen[j - 1] && v <= seqslen[j]) {
                        v = 1 + v - seqslen[j - 1];
                        break;
                    }
                }
            }
            sb.append("  <text class=\"labelMed\" x=\"").append(40 + i * w).append("\" y=\"44\">").append(formatThousands(v)).append("</text>\n");
        }

        // Sequence boundary ticks and names (top-left area)
        if (seqslen != null && seqslen.length > 0) {
            int x1 = 0;
            for (int i = 0; i < seqslen.length; i++) {
                int tx = (int) (x1 * nucleotidesPerPixel);
                int xTick = 50 + tx;
                sb.append("  <line class=\"tick\" x1=\"").append(xTick).append("\" y1=\"1\" x2=\"").append(xTick).append("\" y2=\"55\"/>\n");
                sb.append("  <text class=\"labelBig\" x=\"").append(xTick + 15).append("\" y=\"18\">").append(esc(sname[i])).append("</text>\n");
                sb.append("  <text class=\"labelBig\" x=\"").append(xTick + 5).append("\" y=\"50\">1</text>\n");
                x1 = seqslen[i];
            }
        } else {
            sb.append("  <line class=\"tick\" x1=\"50\" y1=\"1\" x2=\"50\" y2=\"20\"/>\n");
            sb.append("  <text class=\"labelBig\" x=\"65\" y=\"18\">").append(esc(sname[n])).append("</text>\n");
        }

        // Gray/Brown baseline segments at y=60 for each cluster interval (emulating first pass in drawClusters)
        // And colored cluster segments at per-cluster y (see clusterRowY).
        for (int i = 0; i < b; i++) {
            int[] z7 = bb.get(i);

            // baseline (brown)
            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) Math.round(z7[j] * nucleotidesPerPixel);
                int x2;
                if (z7[j + 1] > 0) {
                    x2 = 50 + (int) Math.round((z7[j] + z7[j + 1]) * nucleotidesPerPixel);
                } else {
                    x2 = 50 + (int) Math.round((z7[j] - z7[j + 1]) * nucleotidesPerPixel);
                }
                sb.append("  <line class=\"brown\" x1=\"").append(x1).append("\" y1=\"60\" x2=\"").append(x2).append("\" y2=\"60\" stroke-width=\"").append(dotSize).append("\"/>\n");
            }

            int y = clusterRowY(i, z);
            // colored spans
            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) Math.round(z7[j] * nucleotidesPerPixel);
                int x2;
                String cssClass;
                if (z7[j + 1] > 0) {
                    x2 = 50 + (int) Math.round((z7[j] + z7[j + 1]) * nucleotidesPerPixel);
                    cssClass = (i == 0) ? "darkgreen" : "blue";
                } else {
                    x2 = 50 + (int) Math.round((z7[j] - z7[j + 1]) * nucleotidesPerPixel);
                    cssClass = "red";
                }
                sb.append("  <line class=\"").append(cssClass).append("\" x1=\"").append(x1).append("\" y1=\"").append(y).append("\" x2=\"").append(x2).append("\" y2=\"").append(y).append("\" stroke-width=\"").append(dotSize).append("\"/>\n");
                if (i > 1) {
                    sb.append("  <text class=\"labelSmall\" x=\"").append(x2 + 1).append("\" y=\"").append(y).append("\">").append(i).append("</text>\n");
                }
            }
        }
        sb.append("</svg>\n");
        try (BufferedWriter w1 = new BufferedWriter(new FileWriter(svgfile))) {
            w1.write(sb.toString());
        }
    }
    private void SaveSVGLong(String svgfile, int k, long l, int b, int z, int width, int height, float dotSize, long[] seqslen, ArrayList<long[]> bbL) throws IOException {
        final double nucleotidesPerPixel = (double) width / l;
        if (svgfile.length() == 0) {
            svgfile = filePath + "_1.svg";
            if (nseq == 1) {
                svgfile = filePath + ".svg";
            }
        } else {
            svgfile = svgfile + ".svg";
        }

        System.out.println("Saving SVG " + (width + 100) + "x" + (height + 200) + " : " + svgfile);
        StringBuilder sb = new StringBuilder(1 << 20); // pre-allocate
        // SVG header
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        sb.append("<svg xmlns=\"http://www.w3.org/2000/svg\" ")
                .append("xmlns:xlink=\"http://www.w3.org/1999/xlink\" ")
                .append("width=\"").append(width + 100).append("\" ")
                .append("height=\"").append(height + 200).append("\" ")
                .append("viewBox=\"0 0 ").append(width + 100).append(" ").append(height + 200).append("\">\n");

        // Background
        sb.append("  <rect x=\"0\" y=\"0\" width=\"").append(width + 100).append("\" height=\"").append(height + 200).append("\" fill=\"#FFFFFF\"/>\n");

        // Styles (adjustable)
        sb.append("  <style><![CDATA[\n")
                .append("    .axis { stroke:#000; stroke-width:1; }\n")
                .append("    .tick { stroke:#000; stroke-width:1; }\n")
                .append("    .labelBig { font-family:monospace; font-size:18px; font-weight:bold; fill:#000; }\n")
                .append("    .labelSmall { font-family:monospace; font-size:8px; fill:#000; }\n")
                .append("    .labelMed { font-family:monospace; font-size:16px; fill:#000; }\n")
                .append("    .brown { stroke:#663300; }\n") // Brown
                .append("    .blue { stroke:#0000FF; }\n")
                .append("    .red { stroke:#FF0000; }\n")
                .append("    .darkgreen { stroke:#006600; }\n")
                .append("  ]]></style>\n");

        // Top ruler line
        sb.append("  <line class=\"axis\" x1=\"50\" y1=\"55\" x2=\"").append(50 + width).append("\" y2=\"55\"/>\n");

        // Ruler ticks + numeric labels (emulating drawLinesAndLabels)
        int f = k + 5;
        int w = width / f;
        long d = l / f;

        for (int i = 0; i <= f; i++) {
            int xTick = 50 + i * w;
            sb.append("  <line class=\"tick\" x1=\"").append(xTick).append("\" y1=\"45\" x2=\"").append(xTick).append("\" y2=\"55\"/>\n");
            long v = 1 + (long) i * d;
            if (v > l) {
                v = l;
            }
            // Map cumulative coordinate into per-sequence coordinate if multiple seqs
            if (seqslen != null && seqslen.length > 0) {
                for (int j = 1; j < seqslen.length; j++) {
                    if (v >= seqslen[j - 1] && v <= seqslen[j]) {
                        v = 1 + v - seqslen[j - 1];
                        break;
                    }
                }
            }
            sb.append("  <text class=\"labelMed\" x=\"").append(40 + i * w).append("\" y=\"44\">").append(formatThousandsLong(v)).append("</text>\n");
        }

        // Sequence boundary ticks and names (top-left area)
        if (seqslen != null && seqslen.length > 0) {
            long x1 = 0;
            for (int i = 0; i < seqslen.length; i++) {
                int tx = (int) (x1 * nucleotidesPerPixel);
                int xTick = 50 + tx;
                sb.append("  <line class=\"tick\" x1=\"").append(xTick).append("\" y1=\"1\" x2=\"").append(xTick).append("\" y2=\"55\"/>\n");
                sb.append("  <text class=\"labelBig\" x=\"").append(xTick + 15).append("\" y=\"18\">").append(esc(sname[i])).append("</text>\n");
                sb.append("  <text class=\"labelBig\" x=\"").append(xTick + 5).append("\" y=\"50\">1</text>\n");
                x1 = seqslen[i];
            }
        } else {
            sb.append("  <line class=\"tick\" x1=\"50\" y1=\"1\" x2=\"50\" y2=\"20\"/>\n");
            sb.append("  <text class=\"labelBig\" x=\"65\" y=\"18\">").append(esc(sname[0])).append("</text>\n");
        }

        // Baseline (brown) + colored cluster segments, same layout rules as SaveSVG.
        for (int i = 0; i < b; i++) {
            long[] z7 = bbL.get(i);

            // baseline (brown)
            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) Math.round(z7[j] * nucleotidesPerPixel);
                int x2;
                if (z7[j + 1] > 0) {
                    x2 = 50 + (int) Math.round((z7[j] + z7[j + 1]) * nucleotidesPerPixel);
                } else {
                    x2 = 50 + (int) Math.round((z7[j] - z7[j + 1]) * nucleotidesPerPixel);
                }
                sb.append("  <line class=\"brown\" x1=\"").append(x1).append("\" y1=\"60\" x2=\"").append(x2).append("\" y2=\"60\" stroke-width=\"").append(dotSize).append("\"/>\n");
            }

            int y = clusterRowY(i, z);
            // colored spans
            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) Math.round(z7[j] * nucleotidesPerPixel);
                int x2;
                String cssClass;
                if (z7[j + 1] > 0) {
                    x2 = 50 + (int) Math.round((z7[j] + z7[j + 1]) * nucleotidesPerPixel);
                    cssClass = (i == 0) ? "darkgreen" : "blue";
                } else {
                    x2 = 50 + (int) Math.round((z7[j] - z7[j + 1]) * nucleotidesPerPixel);
                    cssClass = "red";
                }
                sb.append("  <line class=\"").append(cssClass).append("\" x1=\"").append(x1).append("\" y1=\"").append(y).append("\" x2=\"").append(x2).append("\" y2=\"").append(y).append("\" stroke-width=\"").append(dotSize).append("\"/>\n");
                if (i > 1) {
                    sb.append("  <text class=\"labelSmall\" x=\"").append(x2 + 1).append("\" y=\"").append(y).append("\">").append(i).append("</text>\n");
                }
            }
        }
        sb.append("</svg>\n");
        try (BufferedWriter w1 = new BufferedWriter(new FileWriter(svgfile))) {
            w1.write(sb.toString());
        }
    }

    // ── PNG ─────────────────────────────────────────────────────────────────

    private void SaveImage(String pngfile, int k, int n, int l, int b, int z, int width, int height, float dotSize, int[] seqslen) throws IOException {
        final int DPI = 1200;
        final double inchToMeter = 0.0254;
        double nucleotidesPerPixel = (double) width / l;

        if (pngfile.length() == 0) {
            pngfile = filePath + "_" + (n + 1) + ".png";
            if (nseq == 1) {
                pngfile = filePath + ".png";
            }
        } else {
            pngfile = pngfile + ".png";
        }

        System.out.println("Saving picture " + (width + 100) + "x" + (height + 200) + " : " + pngfile);
        BufferedImage image = new BufferedImage(width + 100, height + 200, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();
        g2d.setStroke(new BasicStroke(dotSize));
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, width + 100, height + 200);
        g2d.setColor(Color.BLACK);
        g2d.setFont(new Font("Monospaced", Font.BOLD, 25));

        drawLinesAndLabels(g2d, k, l, width, seqslen);
        if (seqslen.length > 0) {
            int x1 = 0;
            for (int i = 0; i < seqslen.length; i++) {
                x1 = (int) (x1 * nucleotidesPerPixel);
                g2d.drawLine(x1 + 50, 1, x1 + 50, 55);
                g2d.drawString(sname[i], x1 + 65, 18);
                g2d.drawString("1", x1 + 55, 50);
                x1 = seqslen[i];
            }
        } else {
            g2d.drawLine(50, 1, 50, 20);
            g2d.drawString(sname[n], 65, 18);
        }

        drawClusters(g2d, b, z, nucleotidesPerPixel);
        g2d.dispose();

        // PNG writer
        Iterator<ImageWriter> writers = ImageIO.getImageWritersByFormatName("png");
        if (!writers.hasNext()) {
            throw new IllegalStateException("No PNG writer found");
        }

        ImageWriter writer = writers.next();
        File outputFile = new File(pngfile);
        try (ImageOutputStream ios = ImageIO.createImageOutputStream(outputFile)) {
            writer.setOutput(ios);
            ImageWriteParam param = writer.getDefaultWriteParam();

            IIOMetadata metadata = writer.getDefaultImageMetadata(ImageTypeSpecifier.createFromBufferedImageType(BufferedImage.TYPE_INT_RGB), param);
            if (metadata.isReadOnly() || !metadata.isStandardMetadataFormatSupported()) {
                System.err.println("Warning: can't write metadata for DPI");
            } else {
                double pixelsPerMeter = DPI / inchToMeter;
                IIOMetadataNode pHYs_node = new IIOMetadataNode("pHYs");
                pHYs_node.setAttribute("pixelsPerUnitXAxis", Integer.toString((int) pixelsPerMeter));
                pHYs_node.setAttribute("pixelsPerUnitYAxis", Integer.toString((int) pixelsPerMeter));
                pHYs_node.setAttribute("unitSpecifier", "meter");
                IIOMetadataNode root = new IIOMetadataNode("javax_imageio_png_1.0");
                root.appendChild(pHYs_node);

                metadata.mergeTree("javax_imageio_png_1.0", root);
            }
            writer.write(metadata, new IIOImage(image, null, metadata), param);
        }

        writer.dispose();
    }

    private void drawLinesAndLabels(Graphics2D g2d, int k, int l, int width, int[] seqslen) {
        g2d.drawLine(50, 55, width + 50, 55); // top line (x1, y, x2, y)
        int f = k + 5;
        int w = width / f;
        int d = l / f;
        for (int i = 0; i <= f; i++) {
            g2d.drawLine(i * w + 50, 45, i * w + 50, 55);
            int v = 1 + i * d;
            if (v > l) {
                v = l;
            }
            for (int j = 1; j < seqslen.length; j++) {
                if (v >= seqslen[j - 1] && v <= seqslen[j]) {
                    v = 1 + v - seqslen[j - 1];
                    break;
                }
            }
            g2d.drawString(String.format("%,d", v), 63 + i * w, 44);
        }
    }

    /**
     * Vertical position of cluster row {@code i} (row spacing {@code z}), shared
     * by the SVG and PNG renderers so the two output formats lay their rows out
     * identically. Rows 0 and 1 sit at fixed heights near the top; rows ≥2 step
     * down by {@code z}. Previously the PNG path used a different formula
     * (110 / 120 + i*z) than the SVG path (130 / 180 + i*z), so the same run
     * produced visually different figures in the two formats.
     */
    private static int clusterRowY(int i, int z) {
        return (i == 0) ? 90 : (i == 1) ? 130 : 180 + i * z;
    }

    private void drawClusters(Graphics2D g2d, int b, int z, double w1) {
        // Color DarkRed = new Color(153, 0, 0); //https://teaching.csse.uwa.edu.au/units/CITS1001/colorinfo.html
        Color DarkGreen = new Color(0, 102, 0);
        Color Brown = new Color(102, 51, 0);
        g2d.setFont(new Font("Monospaced", Font.PLAIN, 13));
        for (int i = 0; i < b; i++) {
            int[] z7 = bb.get(i);

            // Gray lines at height 22
            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = (z7[j + 1] > 0) ? 50 + (int) ((z7[j] + z7[j + 1]) * w1) : 50 + (int) ((z7[j] - z7[j + 1]) * w1);
                g2d.setColor(Brown);
                g2d.drawLine(x1, 60, x2, 60); // draw dark gray line (x1, y, x2, y)
            }

            int y = clusterRowY(i, z);

            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 50;
                if (z7[j + 1] > 0) {
                    x2 = x2 + (int) ((z7[j] + z7[j + 1]) * w1);
                    g2d.setColor(i == 0 ? DarkGreen : Color.BLUE);
                } else {
                    x2 = x2 + (int) ((z7[j] - z7[j + 1]) * w1);
                    g2d.setColor(Color.RED);
                }
                g2d.drawLine(x1, y, x2, y); // draw blue line
                if (i > 1) {
                    g2d.drawString(String.valueOf(i), x2 + 10, y);
                }
            }

        }
    }
}
