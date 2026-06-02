import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Pangenomic analysis of the repeat clusters produced by the combined
 * (multi-file) TotalRepeats run.
 *
 * The combined clustering groups homologous repeat instances detected across
 * the concatenated set of input sequences into families (CRP clusters). By
 * mapping every instance back to the input sequence it came from (via the
 * cumulative sequence-boundary array {@code seqslen}), each family gets a
 * presence/absence profile over the analyzed genomes. From that profile the
 * families are split into:
 *
 * <ul>
 *   <li><b>core</b> &mdash; present in all genomes (the shared content);</li>
 *   <li><b>accessory</b> (shell) &mdash; present in some but not all genomes;</li>
 *   <li><b>unique</b> (cloud / singleton) &mdash; present in a single genome
 *       (the differing content).</li>
 * </ul>
 *
 * The class is self-contained: it only needs the cluster blocks, the
 * sequence boundaries and the genome names, so it can be reused by any of the
 * combined run modes.
 */
public final class PangenomeAnalysis {

    /** One repeat family (CRP cluster) with its per-genome presence. */
    private static final class Family {
        final int clusterId;        // matches the ClusterID column of the GFF report
        final String refLabel;      // reference annotation, or "-"
        final int[] count;          // instances per genome
        final long[] bp;            // covered bp per genome
        int frequency;              // number of genomes in which it is present
        int totalMembers;
        long totalBp;
        String category;            // core | accessory | unique

        Family(int clusterId, String refLabel, int nGenomes) {
            this.clusterId = clusterId;
            this.refLabel = refLabel;
            this.count = new int[nGenomes];
            this.bp = new long[nGenomes];
        }
    }

    private final int nGenomes;
    private final String[] names;
    private final long[] seqslen;           // cumulative end coordinates per genome (global, long)

    private final List<Family> families = new ArrayList<>();

    // Aggregate, non-family content per genome.
    private final int[] ucrpCount;          // genome-specific unclustered repeats
    private final long[] ucrpBp;
    private final int[] strCount;           // short tandem repeats
    private final long[] strBp;

    private int coreCount, accessoryCount, uniqueCount;
    private long coreBp, accessoryBp, uniqueBp;

    public PangenomeAnalysis(ArrayList<long[]> bb, long[] seqslen, String[] names,
                             int[] refclust, String[] refsname) {
        this.nGenomes = names.length;
        this.names = names;
        this.seqslen = seqslen;
        this.ucrpCount = new int[nGenomes];
        this.ucrpBp = new long[nGenomes];
        this.strCount = new int[nGenomes];
        this.strBp = new long[nGenomes];

        // bb index 0 = STR, index 1 = UCRP, index >= 2 = CRP families.
        for (int i = 0; i < bb.size(); i++) {
            long[] z = bb.get(i);
            if (z == null) {
                continue;
            }
            int clusterId = i + 1;       // same numbering as the GFF ClusterID column

            if (i == 0) {                // STR aggregate
                for (int j = 0; j + 1 < z.length; j += 2) {
                    int g = genomeOf(z[j]);
                    strCount[g]++;
                    strBp[g] += Math.abs(z[j + 1]);
                }
                continue;
            }
            if (i == 1) {                // UCRP: genome-specific unclustered repeats
                for (int j = 0; j + 1 < z.length; j += 2) {
                    int g = genomeOf(z[j]);
                    ucrpCount[g]++;
                    ucrpBp[g] += Math.abs(z[j + 1]);
                }
                continue;
            }

            // CRP homology family
            String refLabel = "-";
            if (refclust != null && clusterId < refclust.length && refclust[clusterId] > 0
                    && refsname != null && refclust[clusterId] - 1 < refsname.length) {
                refLabel = refclust[clusterId] + ":" + refsname[refclust[clusterId] - 1];
            }
            Family f = new Family(clusterId, refLabel, nGenomes);
            for (int j = 0; j + 1 < z.length; j += 2) {
                int g = genomeOf(z[j]);
                f.count[g]++;
                f.bp[g] += Math.abs(z[j + 1]);
            }
            for (int g = 0; g < nGenomes; g++) {
                if (f.count[g] > 0) {
                    f.frequency++;
                }
                f.totalMembers += f.count[g];
                f.totalBp += f.bp[g];
            }
            if (f.frequency == 0) {
                continue;               // defensive: empty family
            }
            if (f.frequency == nGenomes) {
                f.category = "core";
                coreCount++;
                coreBp += f.totalBp;
            } else if (f.frequency == 1) {
                f.category = "unique";
                uniqueCount++;
                uniqueBp += f.totalBp;
            } else {
                f.category = "accessory";
                accessoryCount++;
                accessoryBp += f.totalBp;
            }
            families.add(f);
        }
    }

    /** Maps a global start coordinate to its genome index (same rule as SavingGFF). */
    private int genomeOf(long start) {
        for (int w = 0; w < seqslen.length; w++) {
            if (seqslen[w] > start) {
                return w;
            }
        }
        return nGenomes - 1;            // fallback (coordinate at the very end)
    }

    private int totalFamilies() {
        return families.size();
    }

    private static String pct(long part, long whole) {
        return (whole > 0) ? String.format("%.2f", 100.0 * part / whole) + "%" : "0.00%";
    }

    /** Writes {@code base + "_pangenome.txt"} and {@code base + "_pangenome.tsv"}. */
    public void write(String base) throws IOException {
        writeSummary(base + "_pangenome.txt");
        writeMatrix(base + "_pangenome.tsv");
    }

    private void writeSummary(String file) throws IOException {
        int total = totalFamilies();
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome report: " + file);

            w.write("#TotalRepeats pangenome analysis\n");
            w.write("#Genomes (N=" + nGenomes + "):\n");
            for (int g = 0; g < nGenomes; g++) {
                w.write("#  " + (g + 1) + "\t" + names[g] + "\n");
            }
            w.write("#\n");
            w.write("#Definitions:\n");
            w.write("#  Family (CRP) = cluster of homologous repeat instances across the analyzed sequences.\n");
            w.write("#  Present in a genome = the family has >= 1 instance located in that genome.\n");
            w.write("#  core = present in all " + nGenomes + " genomes; accessory = present in 2.."
                    + (nGenomes - 1) + "; unique = present in exactly 1.\n#\n");

            w.write("== Summary ==\n");
            w.write("Repeat families (CRP clusters): " + total + "\n");
            w.write("  core      (in " + nGenomes + " genomes): " + coreCount
                    + "\t(" + pct(coreCount, total) + ")\tbp=" + coreBp + "\n");
            w.write("  accessory (in 2.." + (nGenomes - 1) + "):       " + accessoryCount
                    + "\t(" + pct(accessoryCount, total) + ")\tbp=" + accessoryBp + "\n");
            w.write("  unique    (in 1 genome):       " + uniqueCount
                    + "\t(" + pct(uniqueCount, total) + ")\tbp=" + uniqueBp + "\n");

            // Soft-core (>= 95% of genomes) - informational, useful when N is large.
            int softThresh = (int) Math.ceil(0.95 * nGenomes);
            int softcore = 0;
            for (Family f : families) {
                if (f.frequency >= softThresh) {
                    softcore++;
                }
            }
            w.write("  soft-core (>= 95% = >= " + softThresh + " genomes): " + softcore + "\n\n");

            // Family frequency spectrum (the classic pangenome U-shape).
            int[] specCount = new int[nGenomes + 1];
            long[] specBp = new long[nGenomes + 1];
            for (Family f : families) {
                specCount[f.frequency]++;
                specBp[f.frequency] += f.totalBp;
            }
            w.write("== Family frequency spectrum ==\n");
            w.write("Present_in_genomes\tFamilies\tBp\n");
            for (int k = 1; k <= nGenomes; k++) {
                w.write(k + "\t" + specCount[k] + "\t" + specBp[k] + "\n");
            }
            w.write("\n");

            // Per-genome summary.
            w.write("== Per-genome summary ==\n");
            w.write("Genome\tFamilies_present\tCore\tAccessory\tUnique_to_this\tCRP_bp\tUCRP_bp\tSTR_bp\n");
            for (int g = 0; g < nGenomes; g++) {
                int present = 0, core = 0, acc = 0, uniq = 0;
                long crpBp = 0;
                for (Family f : families) {
                    if (f.count[g] > 0) {
                        present++;
                        crpBp += f.bp[g];
                        switch (f.category) {
                            case "core" -> core++;
                            case "accessory" -> acc++;
                            case "unique" -> uniq++;
                        }
                    }
                }
                w.write(names[g] + "\t" + present + "\t" + core + "\t" + acc + "\t" + uniq
                        + "\t" + crpBp + "\t" + ucrpBp[g] + "\t" + strBp[g] + "\n");
            }
            w.write("\n");

            // Pairwise shared families (count and Jaccard similarity).
            w.write("== Pairwise shared families (shared / Jaccard) ==\n");
            w.write("genome_i\tgenome_j\tshared\tunion\tjaccard\n");
            for (int a = 0; a < nGenomes; a++) {
                for (int bIdx = a + 1; bIdx < nGenomes; bIdx++) {
                    int shared = 0, union = 0;
                    for (Family f : families) {
                        boolean pa = f.count[a] > 0;
                        boolean pb = f.count[bIdx] > 0;
                        if (pa && pb) {
                            shared++;
                        }
                        if (pa || pb) {
                            union++;
                        }
                    }
                    String jac = (union > 0) ? String.format("%.4f", (double) shared / union) : "0.0000";
                    w.write(names[a] + "\t" + names[bIdx] + "\t" + shared + "\t" + union + "\t" + jac + "\n");
                }
            }
            w.write("\n");

            // Core families.
            w.write("== Core families (common to all genomes) ==\n");
            w.write("ClusterID\tRefLabel\tMembers\tBp\n");
            for (Family f : families) {
                if ("core".equals(f.category)) {
                    w.write(f.clusterId + "\t" + f.refLabel + "\t" + f.totalMembers + "\t" + f.totalBp + "\n");
                }
            }
            w.write("\n");

            // Genome-specific (unique) families - the differing content.
            w.write("== Genome-specific (unique) families ==\n");
            for (int g = 0; g < nGenomes; g++) {
                w.write("[" + names[g] + "]\n");
                w.write("ClusterID\tRefLabel\tMembers\tBp\n");
                for (Family f : families) {
                    if ("unique".equals(f.category) && f.count[g] > 0) {
                        w.write(f.clusterId + "\t" + f.refLabel + "\t" + f.count[g] + "\t" + f.bp[g] + "\n");
                    }
                }
                w.write("\n");
            }
        }
    }

    /** Machine-readable presence/absence matrix (one row per CRP family). */
    private void writeMatrix(String file) throws IOException {
        try (BufferedWriter w = new BufferedWriter(new FileWriter(file))) {
            System.out.println("Saving pangenome matrix: " + file);
            StringBuilder header = new StringBuilder("ClusterID\tType\tCategory\tFrequency\tRefLabel");
            for (int g = 0; g < nGenomes; g++) {
                header.append("\t").append(names[g]);
            }
            header.append("\tTotalMembers\tTotalBp\n");
            w.write(header.toString());

            for (Family f : families) {
                StringBuilder row = new StringBuilder();
                row.append(f.clusterId).append("\tCRP\t").append(f.category)
                        .append("\t").append(f.frequency).append("\t").append(f.refLabel);
                for (int g = 0; g < nGenomes; g++) {
                    row.append("\t").append(f.count[g]);
                }
                row.append("\t").append(f.totalMembers).append("\t").append(f.totalBp).append("\n");
                w.write(row.toString());
            }
        }
    }
}
