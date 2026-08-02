import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;

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

    public void SetPangenome(boolean b) {
        this.pangenome = b;
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

    /** Which measure clusters the blocks: SequencesClustering.MODE_PROFILE / MODE_CONTAIN. */
    public void SetClusterMode(int i) {
        clusterMode = i;
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

    // ── long-coordinate helpers for the combined (>2.1 Gb) run ──────────────────
    /**
     * Copies a LOCAL block array (pairs: start, length; lengths may be negative
     * for the reverse strand) into a GLOBAL long block array, shifting every
     * start by {@code sz}. Lengths (and their sign) are preserved. The input
     * array is not modified, so the caller can keep using the local coordinates
     * afterwards.
     */
    private long[] shiftToGlobal(int[] blocks, long sz) {
        long[] g = new long[blocks.length];
        for (int j = 0; j + 1 < blocks.length; j += 2) {
            g[j] = (long) blocks[j] + sz;   // start → global
            g[j + 1] = blocks[j + 1];       // length (sign preserved)
        }
        return g;
    }

    /**
     * Concatenates two long block arrays.
     */
    private long[] concatLong(long[] a, long[] b) {
        long[] r = new long[a.length + b.length];
        System.arraycopy(a, 0, r, 0, a.length);
        System.arraycopy(b, 0, r, a.length, b.length);
        return r;
    }

    /**
     * Long-coordinate variant: extracts the blocks of
     * a GLOBAL long array that fall inside one sequence's range [seqStart,
     * seqEnd) and remaps them to that sequence's LOCAL int coordinates. A
     * single sequence is always < 2.1 Gb, so the local start/length fit in an
     * int. Negative lengths (reverse strand) are kept.
     */
    private int[] sliceBlocksLocalLong(long[] blocks, long seqStart, long seqEnd) {
        long seqLen = seqEnd - seqStart;
        ArrayList<Integer> out = new ArrayList<>();
        for (int j = 0; j + 1 < blocks.length; j += 2) {
            long start = blocks[j];
            long len = blocks[j + 1];
            long absLen = Math.abs(len);
            if (start >= seqStart && start < seqEnd) {
                long localStart = start - seqStart;
                long localLen = absLen;
                if (localStart + localLen > seqLen) {
                    localLen = seqLen - localStart;   // clip to the sequence end
                }
                if (localLen > 0) {
                    out.add((int) localStart);
                    out.add(len < 0 ? -(int) localLen : (int) localLen);
                }
            }
        }
        int[] res = new int[out.size()];
        for (int t = 0; t < res.length; t++) {
            res[t] = out.get(t);
        }
        return res;
    }

    /**
     * Where a combined run gets its repeat blocks. The rest of the pipeline —
     * global long coordinates, joint clustering, the combined report and the
     * per-file slices — is identical either way, which is why the two modes
     * share {@link #runCombined}.
     */
    private enum CombineSource {
        /** Detect repeats in the sequences themselves ({@code -collate}). */
        SEQUENCES,
        /** Read repeats from soft-masked sequences supplied as input ({@code -combinemask}). */
        MASKS
    }

    /** Genome-wide analysis with each sequence masked individually, then clustered jointly. */
    public void RunSynchronizingClassification(int k, boolean fst) throws IOException {
        runCombined(k, fst, CombineSource.SEQUENCES);
    }

    /** Genome-wide analysis taking previously masked sequences as the input. */
    public void RunCombineMask(int k, boolean fst) throws IOException {
        runCombined(k, fst, CombineSource.MASKS);
    }

    /**
     * The shared body of the two combined modes. They differ only in how each
     * sequence is masked — see {@link CombineSource} — and in whether a per-file
     * mask is written, which only makes sense when this run produced it.
     */
    private void runCombined(int k, boolean fst, CombineSource source) throws IOException {
        startTime = System.nanoTime();

        // Global coordinates are long: the concatenation of all sequences may exceed
        // the ~2.1 Gb int/String limit even though each individual chromosome fits in
        // a String. seqslen[i] is the cumulative global end of sequence i; u2/ssr2 are
        // the masked-repeat / STR blocks remapped to global coordinates.
        long[] seqslen = new long[nseq];
        long[] u2 = new long[0];
        long[] ssr2 = new long[0];

        // Per-sequence statistics, remembered for the individual report headers
        // (the per-file output is built afterwards as a slice of the combined run).
        double[] repStat = new double[nseq];
        double[] ssrStat = new double[nseq];
        double[] gapLenStat = new double[nseq];
        double[] gapPctStat = new double[nseq];

        String[] seqs = seq;     // keep the individual sequences for the per-file reports

        long sz = 0;
        for (int i = 0; i < nseq; i++) {
            repeatslen = 0;
            gapslen = 0;
            int l = seq[i].length();   // a single chromosome always fits in an int

            System.out.println("\n" + sname[i]);
            System.out.println("Target sequence length = " + l + " nt");

            LowComplexitySequence2 m1 = new LowComplexitySequence2();
            // A mask-derived input carries its repeats as lower case, which the STR
            // detector must not read as a case distinction.
            m1.FindAllSSRs(source == CombineSource.MASKS ? seq[i].toLowerCase() : seq[i],
                    telomers, SSRdetection);
            byte[] ssrmsk = m1.MapBytes();
            int[] ssr = m1.IntBlocks();
            long seqStart = System.nanoTime();   // this sequence's own masking clock
            int[] copy = Arrays.copyOf(ssr, ssr.length);   // LOCAL blocks for this file's .msk

            // Append this sequence's STR blocks to the combined run in GLOBAL long
            // coordinates (shift local starts by sz). The local ssr/copy stay untouched.
            ssr2 = concatLong(ssr2, shiftToGlobal(ssr, sz));

            ssrglobal = m1.GetTotalRepeats();

            int[] u;
            if (source == CombineSource.MASKS) {
                // The input is already masked: read the lower-case runs back out.
                MaskResult fc = new MaskResult();
                u = fc.ReadMask(seq[i], gap, minlenseq, ssrmsk);
                repeatslen = fc.getRepeatsLen();
                gapslen = fc.getGaps();
            } else {
                MaskingSequence ms = new MaskingSequence();
                u = ms.mask(seq[i], ssrmsk, kmerln, minlenseq);
                repeatslen = ms.repeatLength();
                gapslen = ms.gapsLength();
            }
            repeatslen = (repeatslen * 100) / (l - gapslen);
            ssrglobal = (ssrglobal * 100) / (l - gapslen);
            gaps = (gapslen * 100) / l;
            System.out.println("Sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
            System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
            System.out.println("Sequence gap (bp)=" + (int) gapslen + " (" + String.format("%.4f", gaps) + "%)");
            maskduration = (System.nanoTime() - seqStart) / 1000000000;
            System.out.println("Masking time taken: " + maskduration + " seconds\n");

            // remember this sequence's own statistics for the individual report
            repStat[i] = repeatslen;
            ssrStat[i] = ssrglobal;
            gapLenStat[i] = gapslen;
            gapPctStat[i] = gaps;

            // Per-file mask is written from LOCAL coordinates over seq[i] — but only
            // when this run did the masking. In -combinemask the masks are the input,
            // so writing them back out would merely copy them.
            if (source == CombineSource.SEQUENCES) {
                filePath = filesPath[i];
                SavingMask3("", i, u, copy);
            }

            // Append this sequence's masked-repeat blocks to the combined run (global).
            u2 = concatLong(u2, shiftToGlobal(u, sz));
            sz = sz + l;            // long accumulation (no int overflow at >2.1 Gb)
            seqslen[i] = sz;
        }

        // Virtual concatenation: a SeqStore replaces String.join("", seq) and never
        // builds a single >2.1 Gb String, so the "Requested string length exceeds VM
        // limit" OutOfMemoryError can no longer occur here.
        SeqStore store = new SeqStore(seqs);
        long l = store.length();

        // Combined cluster table in GLOBAL long coordinates. Index 0 is the STR row
        // (matching the old bb[0]); ClusteringMaskingCombined then appends UCRP + the
        // families (so the row/colour/ClusterID layout is identical to before).
        ArrayList<long[]> bbL = new ArrayList<>();
        bbL.add(ssr2);

        System.out.println("\nClustering started...");
        ClusteringMaskingCombined(store, u2, fst, bbL);

        if (bbL != null) {
            // Combined report + picture use long coordinates and read the sequence
            // through the SeqStore (so SeqShow works across the whole concatenation).
            SavingGFFLong(ReportFilePath, l, seqslen, bbL, store);
            SavingSVGLong(ReportFilePath, k, l, iwidth, iheight, seqslen, bbL);
            SavingPangenomeCombined(ReportFilePath, seqslen, bbL, l);   // pangenome: core / accessory / unique families

            // --- Individual (per-file) reports and pictures ---
            // Built as exact slices of the COMBINED clustering so that each sequence's
            // individual report/picture matches its region in the combined one: the
            // same families keep the same cluster index (hence the same colour, row and
            // ClusterID) and reference labels; only the blocks are restricted to this
            // sequence and remapped to its own LOCAL int coordinates (always < 2.1 Gb),
            // which lets the unchanged int-based savers be reused as-is. refclust stays
            // the combined one (ordering is preserved).
            seq = seqs;                          // individual sequences (for SeqShow)
            long start = 0;
            for (int i = 0; i < nseq; i++) {
                long end = seqslen[i];
                int li = seqs[i].length();

                ArrayList<int[]> bbLocal = new ArrayList<>(bbL.size());
                for (long[] z7 : bbL) {
                    bbLocal.add(sliceBlocksLocalLong(z7, start, end)); // same order/count as combined
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
        long[] seqslen = new long[nseq];   // global end offset of each sequence (for slicing u2)

        long sz = 0;
        for (int i = 0; i < nseq; i++) {
            LowComplexitySequence2 m1 = new LowComplexitySequence2();
            m1.FindAllSSRs(seq[i], telomers, SSRdetection);
            ssr[i] = m1.IntBlocks();
            ssrlen[i] = m1.GetTotalRepeats();
            rssr[i] = m1.MapBytes();
            sz += seq[i].length();
            seqslen[i] = sz;
        }

        // Global cross-sequence repeat masking with the single masker: ONE shared
        // k-mer map over ALL sequences (the same cross-sequence detection the former
        // MaskingSequences gave), returning GLOBAL long blocks. rssr is the per-
        // sequence STR mask. No >2.1 Gb String/array is ever built.
        MaskingSequence ms = new MaskingSequence();
        long[] u2 = ms.maskCombined(seq, rssr, kmerln, minlenseq);
        long repBpAll = ms.repeatLength();   // combined totals (used for the summary line)
        long gapBpAll = ms.gapsLength();

        long start = 0;
        for (int i = 0; i < nseq; i++) {
            long end = seqslen[i];
            int l = seq[i].length();

            // per-sequence masked-repeat blocks, remapped to LOCAL int coordinates.
            // Identical to the old rs.get(i), since maskCombined emits the same blocks.
            int[] u = sliceBlocksLocalLong(u2, start, end);

            // per-sequence repeat / gap recomputed from the slice — equal to the old
            // rplen[i] / gplen[i] (gaps = N/n bases, as MaskingSequences counted them).
            long repBp = 0;
            for (int j = 1; j < u.length; j += 2) {
                repBp += Math.abs(u[j]);
            }
            long gapBp = 0;
            for (int p = 0; p < l; p++) {
                char ch = seq[i].charAt(p);
                if (ch == 'N' || ch == 'n') {
                    gapBp++;
                }
            }

            // same (integer-division) arithmetic and output format as before
            repeatslen = (repBp * 100) / (l - gapBp);
            double gps = (gapBp * 100) / l;
            double ssrln = (ssrlen[i] * 100) / (l - gapBp);

            System.out.println(sname[i]);
            System.out.println("Target sequence length = " + l + " nt");
            System.out.println("Sequence coverage by total repeats=" + String.format("%.2f", repeatslen) + "%");
            System.out.println("Short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrln) + "%");
            System.out.println("Sequence gap (bp)=" + (int) gapBp + " (" + String.format("%.4f", gps) + "%)\n");

            filePath = filesPath[i];
            SavingMask(i, u, ssr[i]);

            start = end;
        }

        // Optional combined summary (now available from the single masker).
        // Remove these two lines to match the original output exactly.
        long effAll = (sz - gapBpAll > 0) ? (sz - gapBpAll) : sz;
        System.out.println("Combined coverage by total repeats=" + String.format("%.2f", (repBpAll * 100.0) / effAll) + "%");

        maskduration = (System.nanoTime() - startTime) / 1000000000;
        System.out.println("Masking time taken: " + maskduration + " seconds\n");
    }


    // ============================================================================
//  RunCombining — GLOBAL (cross-sequence) repeat masking via maskCombined.
// ----------------------------------------------------------------------------
//  Repeat masking is now ONE pass over all sequences with a shared k-mer map
//  (MaskingSequence.maskCombined), so repeats shared between different sequences
//  are detected. STR/SSR detection stays per sequence (intra-sequence, as in
//  RunCombine). No sequence is ever concatenated, so the total may exceed 2.1 Gb.
//
//  MaskingSequences class and no separate ByteStore overload needed — maskCombined
//  already byte-codes each chunk internally via normalise()/Tables.dx2.
// ============================================================================
    public void RunCombining(int k, boolean fst) throws IOException {
        startTime = System.nanoTime();

        long[] seqslen = new long[nseq];
        long[] ssr2 = new long[0];          // STR blocks, GLOBAL long coords
        String[] seqs = seq;                // keep individual sequences for per-file output

        // ── 1) STR detection per sequence (intra-sequence; same as RunCombine) ──
        // ssrmsk[i] is the per-base STR mask fed to maskCombined (so STR positions are
        // excluded from repeat coverage, exactly like the single-sequence path).
        byte[][] ssrmsk = new byte[nseq][];
        long sz = 0, ssrBpAll = 0;
        for (int i = 0; i < nseq; i++) {
            int l = seqs[i].length();       // a single chromosome always fits in an int
            System.out.println("\n" + sname[i]);
            System.out.println("Target sequence length = " + l + " nt");

            LowComplexitySequence2 m1 = new LowComplexitySequence2();
            m1.FindAllSSRs(seqs[i], telomers, SSRdetection);
            ssrmsk[i] = m1.MapBytes();
            ssr2 = concatLong(ssr2, shiftToGlobal(m1.IntBlocks(), sz));
            ssrBpAll += m1.GetTotalRepeats();

            sz += l;
            seqslen[i] = sz;
        }

        // ── 2) ONE global repeat masking over ALL sequences ──
        // Shared k-mer map across every sequence -> cross-sequence repeats detected;
        // no window crosses a junction -> no repeat spans a boundary. Equivalent to
        // "merge into one sequence, block junctions, mask once", but without building
        // a >2.1 Gb String/array. Returns blocks already in GLOBAL long coordinates.
        MaskingSequence ms = new MaskingSequence();
        long[] u2 = ms.maskCombined(seqs, ssrmsk, kmerln, minlenseq);
        long repBpAll = ms.repeatLength();
        long gapBpAll = ms.gapsLength();

        /*      
  // ── 2) ONE global repeat masking: shared k-mer map across ALL sequences ──
    // 'true' = half-kmer sensitivity (matches the single-sequence MaskingSequence).
    MaskingSequences ms = new MaskingSequences();
    ArrayList<int[]> rs = ms.mask(seqs, ssrmsk, kmerln, minlenseq, true);
    long[] rplen = ms.repeatLength();   // per-sequence masked-repeat bp
    long[] gplen = ms.gapsLength();     // per-sequence gap (N) bp
    
    // ── 3) promote per-sequence LOCAL blocks to GLOBAL long coords + aggregate ──
    long[] u2 = new long[0];
    long off = 0, repBpAll = 0, gapBpAll = 0;
    for (int i = 0; i < nseq; i++) {
        u2 = concatLong(u2, shiftToGlobal(rs.get(i), off));
        off += seqs[i].length();
        repBpAll += rplen[i];
        gapBpAll += gplen[i];
    }
         */
        // ── combined header ──
        printCombinedHeader(repBpAll, ssrBpAll, gapBpAll, sz);

        // ── 3) per-file .msk from sliced LOCAL blocks ──
        long start = 0;
        for (int i = 0; i < nseq; i++) {
            long end = seqslen[i];
            int[] uLoc = sliceBlocksLocalLong(u2, start, end);
            int[] ssrLoc = sliceBlocksLocalLong(ssr2, start, end);
            lowercaseAndSaveMsk(i, seqs[i], uLoc, ssrLoc);
            start = end;
        }

        // ── 4) clustering over the virtual concatenation + reports (long pipeline) ──
        SeqStore store = new SeqStore(seqs);
        long l = store.length();
        seq = seqs;

        ArrayList<long[]> bbL = new ArrayList<>();
        bbL.add(ssr2);                                  // index 0 = STR row (canonical layout)
        if (u2.length > 1) {
            System.out.println("Clustering started...");
            ClusteringMaskingCombined(store, u2, fst, bbL);
        }

        if (bbL != null) {
            SavingGFFLong(ReportFilePath, l, seqslen, bbL, store);
            SavingSVGLong(ReportFilePath, k, l, iwidth, iheight, seqslen, bbL);
            SavingPangenomeCombined(ReportFilePath, seqslen, bbL, l);
            savePerFileReports(k, seqs, seqslen, ssr2, u2, bbL);
        }
    }

// ─── shared private helpers (keep RunCombining short) ─────────────────────────
    private void printCombinedHeader(long repBp, long ssrBp, long gapBp, long lenAll) {
        long eff = (lenAll - gapBp > 0) ? (lenAll - gapBp) : lenAll;
        repeatslen = (repBp * 100.0) / eff;
        ssrglobal = (ssrBp * 100.0) / eff;
        gapslen = gapBp;
        gaps = (gapBp * 100.0) / lenAll;
        System.out.println("\nCombined: sequence coverage by repeats=" + String.format("%.2f", repeatslen) + "%");
        System.out.println("Combined: short tandem repeat (STR) sequence coverage=" + String.format("%.2f", ssrglobal) + "%");
        System.out.println("Combined: sequence gap (bp)=" + gapslen + " (" + String.format("%.4f", gaps) + "%)");
    }

// Lowercase the masked-repeat + STR positions over an UPPERCASE copy of seqs[i]
// and write the per-file .msk, same byte format as the previous RunCombining.


// Per-file reports as exact slices of the COMBINED clustering — same families keep
// the same cluster index/colour/ID; identical to the current RunCombining tail.
    private void savePerFileReports(int k, String[] seqs, long[] seqslen,
            long[] ssr2, long[] u2, ArrayList<long[]> bbL) throws IOException {
        seq = seqs;
        long start = 0;
        for (int i = 0; i < nseq; i++) {
            long end = seqslen[i];
            int li = seqs[i].length();

            int[] ssrLocal = sliceBlocksLocalLong(ssr2, start, end);
            int[] uLocal = sliceBlocksLocalLong(u2, start, end);

            long repBp = 0;
            for (int j = 1; j < uLocal.length; j += 2) {
                repBp += Math.abs(uLocal[j]);
            }
            long ssrBp = 0;
            for (int j = 1; j < ssrLocal.length; j += 2) {
                ssrBp += Math.abs(ssrLocal[j]);
            }
            long gapBp = 0;
            for (int pp = 0; pp < li; pp++) {
                char ch = seqs[i].charAt(pp);
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

            ArrayList<int[]> bbLocal = new ArrayList<>(bbL.size());
            for (long[] z7 : bbL) {
                bbLocal.add(sliceBlocksLocalLong(z7, start, end));
            }
            bb = bbLocal;

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

                LowComplexitySequence2 m1 = new LowComplexitySequence2();
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

                LowComplexitySequence2 m1 = new LowComplexitySequence2();
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

    // Single-sequence clustering. A single sequence always fits in one String
    // (and therefore one int-addressable segment), so it is wrapped in a one-segment
    // SeqStore and the int block offsets are promoted to the global long coordinates
    // the new SequencesClustering expects. The long results are safely narrowed back
    // to int for bb (every coordinate here is < 2.1 Gb). For the >2.1 Gb combined
    // case use ClusteringMaskingCombined instead.
    private int ClusteringMasking(String seq, int[] z2, boolean fst) {
        long[][] d;  // d[j][0] = global start; d[j][1] = length
        int[] q;     // cluster ID for each block
        int ncl;

        SeqStore store = new SeqStore(new String[]{seq});
        long[] offsets = new long[z2.length];
        for (int t = 0; t < z2.length; t++) {
            offsets[t] = z2[t];
        }

        SequencesClustering sc = new SequencesClustering(store, refseq, offsets, fst, clusterMode);
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
                    z.add((int) d[i][0]);
                    z.add((int) d[i][1]);
                }
                if (-q[i] == j) {
                    z.add((int) d[i][0]);
                    z.add(-(int) d[i][1]);
                }
            }

            bb.add(z.stream().mapToInt(Integer::intValue).toArray());
        }
        return bb.size();
    }

    // Combined clustering for the virtually-concatenated SeqStore: identical logic
    // to ClusteringMasking, but keeps GLOBAL long coordinates (the combined space can
    // exceed the ~2.1 Gb int limit). Cluster rows are appended to {@code bbL} in the
    // same order as the single-sequence path: index 0 is reserved by the caller for
    // the STR row, this method appends the UCRP row (cluster id 2) and then every
    // family (cluster ids 3..ncl), keeping empty rows as placeholders so that the
    // bb index aligns with the cluster id and with refclust exactly as before.
    private void ClusteringMaskingCombined(SeqStore store, long[] offsets, boolean fst,
            ArrayList<long[]> bbL) {
        SequencesClustering sc = new SequencesClustering(store, refseq, offsets, fst, clusterMode);
        long[][] d = sc.getSequenceOffsets();
        int[] q = sc.getClusterIds();
        refclust = sc.getReferenceIds(); // ID+1 for each cluster
        int ncl = sc.getClusterCount();

        if (ncl < 1) {
            return;
        }

        for (int i = 0; i < q.length; i++) {
            if (q[i] == 0) {
                q[i] = 2;
            }
        }

        for (int j = 2; j <= ncl; j++) {
            ArrayList<Long> z = new ArrayList<>();
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
            long[] row = new long[z.size()];
            for (int t = 0; t < row.length; t++) {
                row[t] = z.get(t);
            }
            bbL.add(row);
        }
    }



    // ===================================================================

    /**
     * The renderer for the current state of this run. A fresh instance is built
     * per figure because {@code bb} and {@code filePath} are reassigned as the
     * run walks through the input files.
     */
    private RepeatFigure figure() {
        return new RepeatFigure(bb, sname, filePath, nseq);
    }

    /**
     * The mask writer for the current file path. Built per file, like the other
     * helpers, so it sees {@code filePath} as it stands at that moment.
     */
    private MaskWriter masks() {
        return new MaskWriter(filePath, seq, sname, nseq, repeatslen);
    }

    private void SavingMask(int n, int[] m, int[] ssr) throws IOException {
        masks().writeMask(n, m, ssr);
    }

    private void SavingMask3(String maskedfile, int n, int[] m, int[] ssr) throws IOException {
        masks().writeMaskTo(maskedfile, n, m, ssr);
    }

    private void SavingMask2(int n, byte[] c, int x1, int x2) throws IOException {
        masks().writeBytes(n, c, x1, x2);
    }

    /**
     * Soft-masks {@code s} over the repeat and STR blocks and writes it as this
     * file's mask. The current file path is repointed first — deliberately, and
     * left pointing there afterwards, because the combined pangenome report
     * falls back to it when no explicit report base is given.
     */
    private void lowercaseAndSaveMsk(int i, String s, int[] u, int[] ssr) throws IOException {
        int l = s.length();
        byte[] ci = s.toUpperCase().getBytes();
        MaskWriter.softMask(ci, u, l);
        MaskWriter.softMask(ci, ssr, l);
        filePath = filesPath[i];
        SavingMask2(i, ci, 0, l);
    }

    /**
     * The annotation writer for the current state of this run. Like
     * {@link #figure()} it is built per report, so it captures the coverage
     * statistics as they stand for the sequence being written.
     */
    private AnnotationWriter annotations() {
        return new AnnotationWriter(
                new AnnotationWriter.Inputs(seq, sname, filesPath, filePath, nseq, bb, refclust, refsname),
                new AnnotationWriter.Stats(kmerln, minlenseq, flanks, gap, SeqShow,
                        repeatslen, ssrglobal, gapslen, gaps, maskduration, startTime, pangenome));
    }

    private void SavingGFF(String reportfile, int n, int l, int[] h) throws IOException {
        annotations().writeTable(reportfile, n, l, h);
    }

    private void SavingGFFLong(String reportfile, long l, long[] h, ArrayList<long[]> bbL, SeqStore store) throws IOException {
        annotations().writeTableLong(reportfile, l, h, bbL, store);
    }

    private void SavingPangenomeCombined(String reportBase, long[] seqslen, ArrayList<long[]> bbL, long l) throws IOException {
        annotations().writePangenome(reportBase, seqslen, bbL, l);
    }

    private void SavingSVG(String reportfile, int k, int n, int len, int dw, int dh, int[] seqslen) throws IOException {
        figure().writeSvg(reportfile, k, n, len, dw, dh, seqslen);
    }

    private void SavingSVGLong(String reportfile, int k, long len, int dw, int dh, long[] seqslen, ArrayList<long[]> bbL) throws IOException {
        figure().writeSvgLong(reportfile, k, len, dw, dh, seqslen, bbL);
    }

    private void SavingPicture(String reportfile, int k, int n, int len, int dw, int dh, int[] seqslen) {
        figure().writePng(reportfile, k, n, len, dw, dh, seqslen);
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
    private int clusterMode = SequencesClustering.MODE_PROFILE;   // 4-mer profile (-vector) by default; -contain selects containment
    private String filePath;
    private String ReportFilePath;
    private int[] refclust;
    private String[] filesPath;
    private String[] seq;
    private String[] sname;
    private String[] refseq;
    private String[] refsname;
    private ArrayList<int[]> bb;
    private boolean pangenome = true;   // generate pangenome (core/accessory/unique) report in combined runs
}
