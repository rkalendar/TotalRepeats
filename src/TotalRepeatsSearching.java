import java.awt.BasicStroke;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import javax.imageio.*;
import javax.imageio.metadata.*;
import javax.imageio.stream.ImageOutputStream;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Iterator;

public final class TotalRepeatsSearching {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = new String[sname.length];
        for (int i = 0; i < sname.length; i++) {
            this.sname[i] = sname[i].split("\\s+", 2)[0];
        }
        nseq = seq.length;
    }

    public void SetRefSequences(String[] refseq, String[] refsname) {
        this.refseq = refseq;
        this.refsname = refsname;
    }

    public void SetReportFile(String a) {
        ReportFilePath = a;
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
        if (this.kmerln > 21) {
            this.kmerln = 21;
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

    public void SetSSRdetection(boolean i) {
        SSRdetection = i;
    }

    public void SetMaskGenerate(boolean i) {
        MaskOnly = i;
    }

    public void SetShowSeq(boolean i) {
        SeqShow = i;
    }

    public void SetImage(int w, int h) {
        if (w > 0) {
            iwidth = w;
        }
        if (h > 0) {
            iheight = h;
        }
    }

    private int[] CombineArrays(int[] u, int[] u2) {
        int newLength = u.length + u2.length;
        int[] combinedArray = new int[newLength];
        System.arraycopy(u, 0, combinedArray, 0, u.length);
        System.arraycopy(u2, 0, combinedArray, u.length, u2.length);
        return combinedArray;
    }

    // Extract the blocks (pairs: start, length) of a combined/global array that
    // fall inside one sequence's range [seqStart, seqEnd) and remap them to that
    // sequence's local coordinates. Negative lengths (reverse strand) are kept.
    private int[] SliceBlocksLocal(int[] blocks, int seqStart, int seqEnd) {
        int seqLen = seqEnd - seqStart;
        ArrayList<Integer> out = new ArrayList<>();
        for (int j = 0; j + 1 < blocks.length; j += 2) {
            int start = blocks[j];
            int len = blocks[j + 1];
            int absLen = Math.abs(len);
            if (start >= seqStart && start < seqEnd) {
                int localStart = start - seqStart;
                int localLen = absLen;
                if (localStart + localLen > seqLen) {
                    localLen = seqLen - localStart;   // clip to the sequence end
                }
                if (localLen > 0) {
                    out.add(localStart);
                    out.add(len < 0 ? -localLen : localLen);
                }
            }
        }
        int[] res = new int[out.size()];
        for (int t = 0; t < res.length; t++) {
            res[t] = out.get(t);
        }
        return res;
    }

    public void RunCombine(int k, boolean fst) throws IOException {
        startTime = System.nanoTime();
        int[] seqslen = new int[nseq];
        int[] u2 = new int[0];
        int[] ssr2 = new int[0];

        // Per-sequence statistics, remembered for the individual report headers
        // (the per-file output is built afterwards as a slice of the combined run).
        double[] repStat = new double[nseq];
        double[] ssrStat = new double[nseq];
        double[] gapLenStat = new double[nseq];
        double[] gapPctStat = new double[nseq];

        String[] seqs = seq;     // keep the individual sequences for the per-file reports

        bb = new ArrayList<>();

        int sz = 0;
        for (int i = 0; i < nseq; i++) {
            repeatslen = 0;
            gapslen = 0;
            int l = seq[i].length();

            System.out.println("\n" + sname[i]);
            System.out.println("Target sequence length = " + l + " nt");

            LowComplexitySequence m1 = new LowComplexitySequence();
            m1.FindAllSSRs(seq[i], telomers, SSRdetection);
            byte[] ssrmsk = m1.MapBytes();
            int[] ssr = m1.IntBlocks();
            int[] copy = Arrays.copyOf(ssr, ssr.length);

            if (i > 0) {
                for (int j = 0; j < ssr.length; j += 2) {
                    ssr[j] = ssr[j] + sz;
                }
            }
            ssr2 = CombineArrays(ssr2, ssr);

            ssrglobal = m1.GetTotalRepeats();
            MaskingSequence ms = new MaskingSequence();

            int[] u = ms.mask(seq[i], ssrmsk, kmerln, minlenseq);
            repeatslen = ms.repeatLength();
            gapslen = ms.gapsLength();
            repeatslen = (repeatslen * 100) / (l - gapslen);
            ssrglobal = (ssrglobal * 100) / (l - gapslen);
            gaps = (gapslen * 100) / l;
            System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
            System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
            System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)");
            maskduration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Masking time taken: " + maskduration + " seconds\n");

            // remember this sequence's own statistics for the individual report
            repStat[i] = repeatslen;
            ssrStat[i] = ssrglobal;
            gapLenStat[i] = gapslen;
            gapPctStat[i] = gaps;

            filePath = filesPath[i];
            SavingMask3("", i, u, copy);

            if (i > 0) {
                for (int j = 0; j < u.length; j += 2) {
                    u[j] = u[j] + sz;
                }
            }
            u2 = CombineArrays(u2, u);
            sz = sz + l;
            seqslen[i] = sz;
        }

        bb.add(ssr2);

        String s = String.join("", seq);
        int l = s.length();
        seq = new String[]{s};
        s = "";
        System.out.println("\nClustering started...");
        ClusteringMasking(seq[0], u2, fst);

        if (bb != null) {
            SavingGFF(ReportFilePath, 0, l, seqslen);
            SavingPicture(ReportFilePath, k, 0, l, iwidth, iheight, seqslen);
            SavingSVG(ReportFilePath, k, 0, l, iwidth, iheight, seqslen);//(int k, int n, int len, int dw, int dh, int[] seqslen)

            // --- Individual (per-file) reports and pictures ---
            // Built as exact slices of the COMBINED clustering so that each
            // sequence's individual report/picture matches its region in the
            // combined one: the same families keep the same cluster index (hence
            // the same colour, row and ClusterID) and reference labels; only the
            // blocks are restricted to this sequence and remapped to its own local
            // coordinates. refclust stays the combined one (ordering is preserved).
            ArrayList<int[]> bbCombined = bb;
            seq = seqs;                          // restore the individual sequences (for SeqShow)
            int start = 0;
            for (int i = 0; i < nseq; i++) {
                int end = seqslen[i];
                int li = seqs[i].length();

                ArrayList<int[]> bbLocal = new ArrayList<>(bbCombined.size());
                for (int[] z7 : bbCombined) {
                    bbLocal.add(SliceBlocksLocal(z7, start, end));   // same order/count as combined
                }
                bb = bbLocal;

                // restore this sequence's own statistics for the report header
                repeatslen = repStat[i];
                ssrglobal = ssrStat[i];
                gapslen = gapLenStat[i];
                gaps = gapPctStat[i];

                filePath = filesPath[i];
                SavingGFF(filesPath[i], i, li, new int[0]);
                SavingPicture(filesPath[i], k, i, li, iwidth, iheight, new int[0]);
                SavingSVG(filesPath[i], k, i, li, iwidth, iheight, new int[0]);

                start = end;
            }
            bb = bbCombined;
        }
    }


    public void RunCombineMask(int k, boolean fst) throws IOException {
        startTime = System.nanoTime();
        int[] seqslen = new int[nseq];
        int[] u2 = new int[0];
        int[] ssr2 = new int[0];

        bb = new ArrayList<>();

        int sz = 0;
        for (int i = 0; i < nseq; i++) {
            repeatslen = 0;
            gapslen = 0;
            int l = seq[i].length();

            startTime = System.nanoTime();
            System.out.println("Target sequence length = " + l + " nt");
            System.out.println("\n" + sname[i]);
            System.out.println("Target sequence length = " + l + " nt");

            LowComplexitySequence m1 = new LowComplexitySequence();
            m1.FindAllSSRs(seq[i].toLowerCase(), telomers, SSRdetection);
            byte[] ssrmsk = m1.MapBytes();
            int[] ssr = m1.IntBlocks();
            ssrglobal = m1.GetTotalRepeats();
            bb.add(ssr);

            MaskResult fc = new MaskResult();
            int[] u = fc.ReadMask(seq[i], gap, minlenseq, ssrmsk);
            repeatslen = fc.getRepeatsLen();
            gapslen = fc.getGaps();

            if (i > 0) {
                for (int j = 0; j < ssr.length; j += 2) {
                    ssr[j] = ssr[j] + sz;
                }
            }
            ssr2 = CombineArrays(ssr2, ssr);

            repeatslen = (repeatslen * 100) / (l - gapslen);
            ssrglobal = (ssrglobal * 100) / (l - gapslen);
            gaps = (gapslen * 100) / l;
            System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
            System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
            System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)");

            maskduration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Masking time taken: " + maskduration + " seconds\n");

            if (i > 0) {
                for (int j = 0; j < u.length; j += 2) {
                    u[j] = u[j] + sz;
                }
            }
            u2 = CombineArrays(u2, u);

            sz = sz + l;
            seqslen[i] = sz;
        }

        bb.add(ssr2);

        String s = String.join("", seq);
        int l = s.length();
        seq = new String[]{s};
        s = "";
        System.out.println("\nClustering started...");
        ClusteringMasking(seq[0], u2, fst);

        if (bb != null) {
            SavingGFF(ReportFilePath, 0, l, seqslen);
            SavingPicture(ReportFilePath, k, 0, l, iwidth, iheight, seqslen);
        }
    }

    public void RunUniquesMaskSaving(int gap, int minLenSeq) throws IOException {
        for (int i = 0; i < nseq; i++) {

            // Compute runs of UPPERCASE (unique) once
            MaskResult mr = new MaskResult();
            int[] runs = mr.ReadUpperMask(seq[i], gap, minLenSeq); // [start0, len0, start1, len1, ...]

            // Pick output path
            final String maskedFile = (nseq == 1) ? (filePath + ".report") : (filePath + "_" + (i + 1) + ".report");
            System.out.println("Saving report file: " + maskedFile);
            Path out = Path.of(maskedFile);

            // Ensure parent directory exists (if any)
            Path parent = out.getParent();
            if (parent != null) {
                Files.createDirectories(parent);
            }

            // Buffered ASCII writer for speed and clarity
            try (BufferedWriter w = Files.newBufferedWriter(out, StandardCharsets.US_ASCII)) {
                w.write(filePath);
                w.write("\n");
                for (int j = 0; j + 1 < runs.length; j += 2) {
                    int start = runs[j];        // 0-based, inclusive
                    int len = runs[j + 1];

                    // Keep the threshold logic consistent with ReadUpperMask (>=)
                    if (len > minLenSeq && start > -1 && start + len <= seq[i].length()) {
                        // 1-based, inclusive coordinates for the header
                        int start1 = start + 1;
                        int end1 = start + len; // inclusive end, correct (not start + len - 1 + 1)
                        // Header
                        // Example: >12345-67890 
                        w.write(">");
                        w.write(Integer.toString(start1));
                        w.write("-");
                        w.write(Integer.toString(end1));
                        w.newLine();
                        w.write(seq[i], start, len);
                        w.newLine();
                    }
                }
            }
        }
    }

    public void RunHomologyMasking(int k) throws IOException {
        startTime = System.nanoTime();
        byte[][] rssr = new byte[nseq][];
        int[][] ssr = new int[nseq][];
        int[] ssrlen = new int[nseq];

        for (int i = 0; i < nseq; i++) {
            LowComplexitySequence m1 = new LowComplexitySequence();
            m1.FindAllSSRs(seq[i], telomers, SSRdetection);
            ssr[i] = m1.IntBlocks();
            ssrlen[i] = m1.GetTotalRepeats();
            rssr[i] = m1.MapBytes();
        }

        MaskingSequences ms = new MaskingSequences();
        ArrayList<int[]> rs = ms.mask(seq, rssr, kmerln, minlenseq, true);

        long[] rplen = ms.repeatLength();
        long[] gplen = ms.gapsLength();

        for (int i = 0; i < nseq; i++) {
            int[] u = rs.get(i);
            int l = seq[i].length();

            repeatslen = (rplen[i] * 100) / (l - gplen[i]);
            double gps = (gplen[i] * 100) / l;
            double ssrln = (ssrlen[i] * 100) / (l - gplen[i]);

            System.out.println(sname[i]);
            System.out.println("Target sequence length = " + l + " nt");
            System.out.println("Sequence coverage by total repeats=" + String.format("%.2f", repeatslen) + "%");
            System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrln) + "%");
            System.out.println("Sequence gap (bp)=" + (int) gplen[i] + " (" + String.format("%.4f", gps) + "%)\n");

            filePath = filesPath[i];
            SavingMask(i, u, ssr[i]);
        }

        maskduration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Masking time taken: " + maskduration + " seconds\n");

    }

    private void SavingMask2(int n, byte[] c, int x1, int x2) throws IOException {
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

    public void RunCombine2(int k, boolean fst) throws IOException {
        startTime = System.nanoTime();
        int[] seqslen = new int[nseq];
        int sz = 0;
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            sz = sz + l;
            seqslen[i] = sz; // end coordinate 
            System.out.println("\n" + sname[i]);
            System.out.println("Target sequence length = " + l + " nt");
        }

        String[] seqs = seq;            // keep individual sequences for per-file output
        String s = String.join("", seq);
        int l = s.length();
        seq = new String[]{s};
        byte[] c = s.toUpperCase().getBytes();
        s = "";

        repeatslen = 0;
        gapslen = 0;
        bb = new ArrayList<>();

        LowComplexitySequence m1 = new LowComplexitySequence();
        m1.FindAllSSRs(seq[0], telomers, SSRdetection);
        byte[] ssrmsk = m1.MapBytes();
        int[] ssr = m1.IntBlocks();
        ssrglobal = m1.GetTotalRepeats();

        for (int i = 0; i < nseq - 1; i++) {
            for (int j = seqslen[i] - gap - 1; j < seqslen[i]; j++) {
                ssrmsk[j] = 99;
            }
        }

        bb.add(ssr);
        MaskingSequence ms = new MaskingSequence();

        int[] u = ms.mask(seq[0], ssrmsk, kmerln, minlenseq);
        repeatslen += ms.repeatLength();
        gapslen = ms.gapsLength();
        repeatslen = (repeatslen * 100) / (l - gapslen);
        ssrglobal = (ssrglobal * 100) / (l - gapslen);
        gaps = (gapslen * 100) / l;

        System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
        System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
        System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)\n");
        for (int j = 0; j < u.length; j += 2) {
            int from = Math.max(0, u[j]);
            int to = Math.min(l, u[j] + u[j + 1]);
            for (int i = from; i < to; i++) {
                byte b = c[i];
                if (b >= 'A' && b < 'Z') {
                    c[i] = (byte) (b + 32);
                }
            }
        }
        for (int j = 0; j < ssr.length; j += 2) {
            int from = Math.max(0, ssr[j]);
            int to = Math.min(l, ssr[j] + ssr[j + 1]);
            for (int i = from; i < to; i++) {
                byte b = c[i];
                if (b >= 'A' && b < 'Z') {
                    c[i] = (byte) (b + 32);
                }
            }
        }

        int x1 = 0;
        for (int i = 0; i < nseq; i++) {
            filePath = filesPath[i];
            SavingMask2(i, c, x1, seqslen[i]);
            x1 = seqslen[i];
        }
        if (u.length > 1) {
            System.out.println("Clustering started...");
            ClusteringMasking(seq[0], u, fst);
        }
        if (bb != null) {
            SavingGFF(ReportFilePath, 0, l, seqslen);
            SavingPicture(ReportFilePath, k, 0, l, iwidth, iheight, seqslen);
            SavingSVG(ReportFilePath, k, 0, l, iwidth, iheight, seqslen);//(int k, int n, int len, int dw, int dh, int[] seqslen)                    
        }

        // --- Individual (per-file) report and picture ---
        // The combined output above is now fully written. Partition the global
        // repeat blocks (u) and STR blocks (ssr) by sequence range, remap them to
        // local coordinates, cluster each sequence on its own and save a separate
        // GFF report, PNG and SVG per input file. seq is restored to the original
        // per-sequence array so SavingGFF/SavingPicture/SavingSVG see seq[i]/sname[i].
        seq = seqs;
        int start = 0;
        for (int i = 0; i < nseq; i++) {
            int end = seqslen[i];               // global end (exclusive)
            int li = seqs[i].length();

            int[] ssrLocal = SliceBlocksLocal(ssr, start, end);
            int[] uLocal = SliceBlocksLocal(u, start, end);

            // Per-file statistics (gaps counted as N/n bases in this sequence)
            long repBp = 0;
            for (int j = 1; j < uLocal.length; j += 2) {
                repBp += Math.abs(uLocal[j]);
            }
            long ssrBp = 0;
            for (int j = 1; j < ssrLocal.length; j += 2) {
                ssrBp += Math.abs(ssrLocal[j]);
            }
            long gapBp = 0;
            for (int p = 0; p < li; p++) {
                char ch = seqs[i].charAt(p);
                if (ch == 'N' || ch == 'n') {
                    gapBp++;
                }
            }
            long eff = (li - gapBp > 0) ? (li - gapBp) : li;
            repeatslen = (repBp * 100.0) / eff;
            ssrglobal = (ssrBp * 100.0) / eff;
            gapslen = gapBp;
            gaps = (gapBp * 100.0) / li;

            System.out.println("\n" + sname[i]);
            System.out.println("Target sequence length = " + li + " nt");
            System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
            System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
            System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)");

            bb = new ArrayList<>();
            bb.add(ssrLocal);                   // STR blocks (own coordinates)
            if (uLocal.length > 1) {
                System.out.println("Clustering started (" + sname[i] + ")...");
                ClusteringMasking(seqs[i], uLocal, fst);
            }
            filePath = filesPath[i];
            SavingGFF(filesPath[i], i, li, new int[0]);
            SavingPicture(filesPath[i], k, i, li, iwidth, iheight, new int[0]);
            SavingSVG(filesPath[i], k, i, li, iwidth, iheight, new int[0]);

            start = end;
        }
    }

    public void RunThroughMask(int k, boolean fst) throws IOException {
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            repeatslen = 0;
            ssrglobal = 0;
            gapslen = 0;
            bb = new ArrayList<>();

            if (l > minlenseq) {
                startTime = System.nanoTime();
                System.out.println("Target sequence length = " + l + " nt");

                LowComplexitySequence m1 = new LowComplexitySequence();
                m1.FindAllSSRs(seq[i].toLowerCase(), telomers, SSRdetection);
                byte[] ssrmsk = m1.MapBytes();
                int[] ssr = m1.IntBlocks();
                ssrglobal = m1.GetTotalRepeats();
                bb.add(ssr);

                MaskResult fc = new MaskResult();
                int[] u = fc.ReadMask(seq[i], gap, minlenseq, ssrmsk);
                repeatslen = fc.getRepeatsLen();
                gapslen = fc.getGaps();

                repeatslen = (repeatslen * 100) / (l - gapslen);
                ssrglobal = (ssrglobal * 100) / (l - gapslen);
                gaps = (gapslen * 100) / l;

                System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
                System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
                System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)\n");

                maskduration = (System.nanoTime() - startTime) / 1000000000;
                System.out.println("Time taken for mask analysis: " + maskduration + " seconds\n");

                if (u.length > 1) {
                    System.out.println("Clustering started...");
                    ClusteringMasking(seq[i], u, fst);
                }
                if (bb != null) {
                    SavingGFF("", i, l, new int[0]);
                    SavingPicture("", k, i, l, iwidth, iheight, new int[0]);
                    SavingSVG("", k, i, l, iwidth, iheight, new int[0]);//(int k, int n, int len, int dw, int dh, int[] seqslen)                            
                }
            }
        }
    }

    public void RunThroughGFF(String inputGFFfile, int k) throws IOException {
        seq = new String[1];
        sname = new String[1];
        nseq = 1;
        int[] x = new int[1];
        startTime = System.nanoTime();

        OpenGFFfile fc = new OpenGFFfile(inputGFFfile);
        bb = fc.getData();
        x[0] = fc.getSeqLen();
        sname[0] = fc.getName();
        System.out.println("Target sequence length = " + x[0] + " nt");

        if (bb != null) {
            SavingPicture("", k, 0, x[0], iwidth, iheight, x);
            SavingSVG("", k, 0, x[0], iwidth, iheight, x);//(int k, int n, int len, int dw, int dh, int[] seqslen)               
        }

    }

    public void Run(int k, boolean fst) throws IOException {
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            repeatslen = 0;
            ssrglobal = 0;
            gapslen = 0;
            bb = new ArrayList<>();

            if (l > minlenseq) {
                startTime = System.nanoTime();
                System.out.println("Target sequence length = " + l + " nt");

                LowComplexitySequence2 m1 = new LowComplexitySequence2();
                m1.FindAllSSRs(seq[i], telomers, SSRdetection);
                byte[] ssrmsk = m1.MapBytes();
                int[] ssr = m1.IntBlocks();

                ssrglobal = m1.GetTotalRepeats();

                bb.add(ssr);

                MaskingSequence ms = new MaskingSequence();

                int[] u = ms.mask(seq[i], ssrmsk, kmerln, minlenseq);
                repeatslen += ms.repeatLength();
                gapslen = ms.gapsLength();
                repeatslen = (repeatslen * 100) / (l - gapslen);
                ssrglobal = (ssrglobal * 100) / (l - gapslen);
                gaps = (gapslen * 100) / l;

                System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
                System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
                System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)\n");

                maskduration = (System.nanoTime() - startTime) / 1000000000;
                System.out.println("Masking time taken: " + maskduration + " seconds\n");

                //   if (u.length > 1) {
                if (MaskOnly) {
                    SavingMask(i, u, ssr);
                } else {
                    SavingMask(i, u, ssr);
                    System.out.println("Clustering started...");
                    ClusteringMasking(seq[i], u, fst);
                }
                // }

                if (!MaskOnly && bb != null) {
                    SavingGFF("", 0, l, new int[0]);
                    SavingPicture("", k, i, l, iwidth, iheight, new int[0]);
                    SavingSVG("", k, i, l, iwidth, iheight, new int[0]);//(int k, int n, int len, int dw, int dh, int[] seqslen)
                }
            }
        }
    }

    public void RunAlignmentMask(int k, boolean fst) throws IOException {
        for (int i = 0; i < nseq; i++) {
            int l = seq[i].length();
            repeatslen = 0;
            ssrglobal = 0;
            gapslen = 0;
            bb = new ArrayList<>();

            if (l > minlenseq) {
                startTime = System.nanoTime();
                System.out.println("Target sequence length = " + l + " nt");

                LowComplexitySequence m1 = new LowComplexitySequence();
                m1.FindAllSSRs(seq[i], telomers, SSRdetection);
                int[] ssr = m1.IntBlocks();

                ssrglobal = m1.GetTotalRepeats();

                bb.add(ssr);

                MaskingPairwiseAlignmentSequence ms = new MaskingPairwiseAlignmentSequence();
                int[] u = ms.mask(seq[i], kmerln, minlenseq);
                SaveMask(i, ms.getByteMask());

                repeatslen += ms.repeatLength();
                gapslen = l - ms.noGapsLength();
                repeatslen = (repeatslen * 100) / (l - gapslen);
                ssrglobal = (ssrglobal * 100) / (l - gapslen);
                gaps = (gapslen * 100) / l;

                System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
                System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
                System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)\n");

                maskduration = (System.nanoTime() - startTime) / 1000000000;
                System.out.println("Masking time taken: " + maskduration + " seconds\n");

                if (u.length > 1) {
                    if (!MaskOnly) {
                        System.out.println("Clustering started...");
                        ClusteringMasking(seq[i], u, fst);
                    }
                }

                if (!MaskOnly && bb != null) {
                    SavingGFF("", 0, l, new int[0]);
                    SavingPicture("", k, i, l, iwidth, iheight, new int[0]);
                    SavingSVG("", k, i, l, iwidth, iheight, new int[0]);//(int k, int n, int len, int dw, int dh, int[] seqslen)
                }
            }
        }
    }

    private void SaveMask(int n, byte[] m) throws IOException {
        int l = seq[n].length();
        if (gapslen == 0) {
            gapslen = l;
        }
        double z = gapslen;
        double v = gapslen;

        long duration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Time taken: " + duration + " seconds\n");

        String maskedfile = (nseq == 1) ? filePath + ".msk" : filePath + "_" + (n + 1) + ".msk";

        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            System.out.println("Saving masked file: " + maskedfile);

            byte[] c = seq[n].getBytes();
            for (int i = 0; i < l; i++) {
                if (m[i] == 0) {
                    z--;
                    c[i] = (byte) (c[i] - 32);
                }
            }
            z = (z * 100 / v);
            System.out.println("Sequence coverage by repeats = " + String.format("%.2f", z) + "%");
            fileWriter.write(">" + sname[n] + " TotalRepeats: Sequence coverage by repeats = " + String.format("%.2f", z) + "%\n");
            String seqStr = new String(c);
            for (int i = 0; i < seqStr.length(); i += 70) {
                int end = Math.min(i + 70, seqStr.length());
                fileWriter.write(seqStr.substring(i, end));
                fileWriter.write("\n");
            }

        }
    }

    private int ClusteringMasking(String seq, int[] z2, boolean fst) {
        int[][] d;   // d[j][0] = x1; d[j][1] = length;
        int[] q;     // cluster ID for each block
        int ncl;

        SequencesClustering sc = new SequencesClustering(seq, refseq, z2, fst);
        d = sc.getSequenceOffsets();
        q = sc.getClusterIds();
        refclust = sc.getReferenceIds(); // ID+1 for each cluster
        ncl = sc.getClusterCount();

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

    private void SavingMask(int n, int[] m, int[] ssr) throws IOException {
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

    private void SavingMask3(String maskedfile, int n, int[] m, int[] ssr) throws IOException {
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

    private void SavingSVG(String reportfile, int k, int n, int len, int dw, int dh, int[] seqslen) throws IOException {
        int maxClusters = 500;
        int maxImageDimension = 120000;
        int minImageWidth = 4000;
        int minImageHeight = 100;
        int stepPadding = 20;
        int b = Math.min(bb.size(), maxClusters);
        int z = calculateClusterStep(b);
        float dotSize = calculateDotSize(b);
        int width = calculateWidth(k, len, dw, maxImageDimension, minImageWidth);
        int height = calculateHeight(b, z, dh, maxImageDimension, minImageHeight, stepPadding);

        SaveSVG(reportfile, k, n, len, b, z, width, height, dotSize, seqslen);
    }

// Escapes XML-sensitive characters for text nodes/attributes
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
        // And colored cluster segments at per-cluster y.
        // y position rule matches your code: y = (i > 10) ? 190 + i*z : 88 + i*20;
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

            //  int y = (i > 10) ? 190 + (i * z) : 88 + (i * 20);
            int y = (i == 0) ? 90 : (i == 1) ? 130 : 180 + i * z;
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

    private void SavingGFF(String reportfile, int n, int l, int[] h) throws IOException {
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
                            s0 = dna.ComplementDNA2(s0);
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

    private void SavingPicture(String reportfile, int k, int n, int len, int dw, int dh, int[] seqslen) throws IOException {
        int maxClusters = 500;
        int maxImageDimension = 120000;
        int minImageWidth = 4000;
        int minImageHeight = 100;
        int stepPadding = 20;

        // Adjust number of clusters `b`
        int b = Math.min(bb.size(), maxClusters); // Maximum of 1000 clusters

        // Adjust `z` (step between clusters) based on `b`
        int z = calculateClusterStep(b);
        // Calculate dot size
        float dotSize = calculateDotSize(b);
        //float dotSize =10;//'z;   

        // Calculate width and height
        int width = calculateWidth(k, len, dw, maxImageDimension, minImageWidth);
        int height = calculateHeight(b, z, dh, maxImageDimension, minImageHeight, stepPadding);

        try {
            SaveImage(reportfile, k, n, len, b, z, width, height, dotSize, seqslen);
        } catch (IOException e) {
            System.out.println("Error saving the picture: " + e.getMessage());
        }
    }

    private int calculateClusterStep(int b) {
        if (b > 500) {
            return 4;
        }
        if (b > 400) {
            return 5;
        }
        if (b > 300) {
            return 6;
        }
        if (b > 200) {
            return 10;
        }
        if (b > 100) {
            return 12;
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

            int y = (i == 0) ? 90 : (i == 1) ? 110 : 120 + i * z;
            //int y = 90 + i * z;

            for (int j = 0; j < z7.length - 1; j += 2) {
                int x1 = 50 + (int) (z7[j] * w1);
                int x2 = 50;
                if (z7[j + 1] > 0) {
                    x2 = x2 + (int) ((z7[j] + z7[j + 1]) * w1);
                    g2d.setColor(
                            i == 0 ? DarkGreen
                                    : i > 1 ? Color.BLUE
                                            : Color.DARK_GRAY
                    );
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
    private long maskduration = 0;
    private double gaps = 0;
    private double gapslen = 0;
    private double repeatslen = 0;
    private double ssrglobal = 0;
    private long startTime;
    private int nseq = 0;
    private int iwidth = 0;
    private int iheight = 0;
    private int minlenseq = 90;      // Minimal repeat block size
    private int kmerln = 19;         // kmer=12-21
    private int flanks = 20;
    private int gap = 21;            // gap between repeat blocks, gap=kmer
    private final int telomers = 14; // Kmax=11 ->SSR  Kmax=14 -> telomers //=17
    private boolean SeqShow;
    private boolean SSRdetection = true;
    private boolean MaskOnly;
    private String filePath;
    private String ReportFilePath;
    private int[] refclust;
    private String[] filesPath;
    private String[] seq;
    private String[] sname;
    private String[] refseq;
    private String[] refsname;
    private ArrayList<int[]> bb;
}
