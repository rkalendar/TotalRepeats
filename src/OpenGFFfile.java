import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class OpenGFFfile {

    private final ArrayList<int[]> bb;
    private int seqlen;
    private String name;

    public int getSeqLen() {
        return seqlen;
    }

    public String getName() {
        return name;
    }

    public ArrayList<int[]> getData() {
        return bb;
    }

    public OpenGFFfile(String gfffile) throws IOException {
        // Bucket pairs (start, stop/±stop) by ClusterID, preserving first-seen order
        Map<Integer, List<Integer>> buckets = new LinkedHashMap<>();
        try (BufferedReader br = Files.newBufferedReader(Paths.get(gfffile))) {
            String line;
            name = "";
            while ((line = br.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }

                if (line.startsWith("Sequence length (bp)=")) {
                    String[] d = line.split("=");
                    seqlen = Integer.parseInt(d[1].trim());
                    continue;
                }

//Repeats search for:
//Seqid	Repeat	ClusterID	Start	Stop	Length	Strand	Phase
//NC_017854.1 Pyricularia oryzae 70-15 chromosome 7, whole genome shotgun sequence	STR	1	1	170	170	+	
                if (line.startsWith("Seqid")) {
                    continue;
                }

                String[] d = line.split("\t");
                if (d.length < 7) {
                    continue;  
                }
                // Columns (0-based):
                // 0: Seqid, 1: Repeat, 2: ClusterID, 3: Start, 4: Stop, 5: Length, 6: Strand, 7: Phase (optional)
                boolean plusStrand = "+".equals(d[6].trim());
                int clusterId = Integer.parseInt(d[2].trim());
                int x = Integer.parseInt(d[3].trim()) - 1;     // 0-based start originally
                int ln = Integer.parseInt(d[5].trim());        // length              
                int signedStop = plusStrand ? ln : -ln;
                if ("".equals(name)) {
                    name = d[0];
                }

                List<Integer> bucket = buckets.computeIfAbsent(clusterId, k -> new ArrayList<>());
                bucket.add(x);
                bucket.add(signedStop);
            }
        }

        // Convert buckets to int[] arrays
        bb = new ArrayList<>(buckets.size());
        for (List<Integer> vals : buckets.values()) {
            int[] arr = vals.stream().mapToInt(Integer::intValue).toArray();
            bb.add(arr);
        }
    }
}
