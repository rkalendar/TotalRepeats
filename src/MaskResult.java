import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MaskResult {

    private int gapslen;
    private int repeatslen;

    public MaskResult() {
    }

    public int getGaps() {
        return gapslen;
    }

    public int getRepeatsLen() {
        return repeatslen;
    }

    public int[] ReadMask(String seq, int gap, int minLenSeq, byte[] ssrmsk) {

        byte[] b = seq.getBytes();
        int n = b.length;

        List<Integer> runs = new ArrayList<>();
        int seqIdx = -1;          // increases only when we see letters
        boolean inLower = false;
        int runStart = -1;

        for (int i = 0; i < n; i++) {
            byte c = b[i];
            boolean isUpper = (c > 64 && c < 90); //https://www.rapidtables.com/code/text/ascii-table.html
            boolean isLower = (c > 96 && c < 122);
            if (c == 78 || c == 110) {
                gapslen++;
            }
            if (isUpper || isLower) {
                seqIdx++;
                if (isLower) {
                    repeatslen++;
                    if (ssrmsk[seqIdx] == 0) {
                        if (!inLower) {
                            inLower = true;
                            runStart = seqIdx;
                        }
                    }

                } else { // uppercase
                    if (inLower) {
                        if (runs.size() > 1) {
                            int e = runs.get(runs.size() - 1);
                            if (e + gap > seqIdx - 1) {
                                runs.set(runs.size() - 1, seqIdx - 1);
                            } else {
                                runs.add(runStart);
                                runs.add(seqIdx - 1); // end of lowercase run is previous letter
                            }
                        } else {
                            runs.add(runStart);
                            runs.add(seqIdx - 1); // end of lowercase run is previous letter
                        }
                        inLower = false;
                    }
                }
            }
        }
        if (inLower) { // tail run
            runs.add(runStart);
            runs.add(seqIdx);
        }

        if (runs.isEmpty()) {
            return new int[0];
        }

        // Merge runs if the gap (number of uppercase letters between them) is <= gap
        // i.e., if next.start - curr.end - 1 <= gap
        List<Integer> merged = new ArrayList<>(runs.size());
        int curS = runs.get(0);
        int curE = runs.get(1);
        for (int p = 2; p < runs.size(); p += 2) {
            int s = runs.get(p);
            int e = runs.get(p + 1);
            if (curE + gap > s) {
                curE = e; // extend
            } else {
                merged.add(curS);
                merged.add(curE);
                curS = s;
                curE = e;
            }
        }
        merged.add(curS);
        merged.add(curE);

        // Filter by min length and emit result
        int[] out = new int[merged.size()];
        int w = 0;
        for (int p = 0; p < merged.size(); p += 2) {
            int s = merged.get(p);
            int e = merged.get(p + 1);
            int len = 1 + e - s;
            if (len >= minLenSeq) {
                out[w++] = s;
                out[w++] = len; // length
            }
        }

        return Arrays.copyOf(out, w);
    }

    public int[] ReadUpperMask(String seq, int gap, int minLenSeq) {
        byte[] b = seq.getBytes();
        int n = b.length;

        List<Integer> runs = new ArrayList<>();
        int seqIdx = -1;          // increases only when we see letters
        boolean inUpper = false;
        int runStart = -1;

        for (int i = 0; i < n; i++) {
            byte c = b[i];
            boolean isUpper = (c >= 65 && c <= 90);   // A–Z
            boolean isLower = (c >= 97 && c <= 122);  // a–z

            if (isUpper || isLower) {
                seqIdx++;
                if (isUpper) { // now uppercase = repeats           
                    if (!inUpper) {
                        inUpper = true;
                        runStart = seqIdx;
                    }

                } else { // lowercase
                    if (inUpper) {
                        if (runs.size() > 1) {
                            int e = runs.get(runs.size() - 1);
                            if (e + gap > seqIdx - 1) {
                                runs.set(runs.size() - 1, seqIdx - 1);
                            } else {
                                runs.add(runStart);
                                runs.add(seqIdx - 1);
                            }
                        } else {
                            runs.add(runStart);
                            runs.add(seqIdx - 1);
                        }
                        inUpper = false;
                    }
                }
            }
        }
        if (inUpper) { // tail run
            runs.add(runStart);
            runs.add(seqIdx);
        }

        if (runs.isEmpty()) {
            return new int[0];
        }

        // Merge runs if the gap (number of lowercase letters between them) is <= gap
        List<Integer> merged = new ArrayList<>(runs.size());
        int curS = runs.get(0);
        int curE = runs.get(1);
        for (int p = 2; p < runs.size(); p += 2) {
            int s = runs.get(p);
            int e = runs.get(p + 1);
            if (curE + gap > s) {
                curE = e;
            } else {
                merged.add(curS);
                merged.add(curE);
                curS = s;
                curE = e;
            }
        }
        merged.add(curS);
        merged.add(curE);

        // Filter by min length and emit result
        int[] out = new int[merged.size()];
        int w = 0;
        for (int p = 0; p < merged.size(); p += 2) {
            int s = merged.get(p);
            int e = merged.get(p + 1);
            int len = 1 + e - s;
            if (len >= minLenSeq) {
                out[w++] = s;
                out[w++] = len;
            }
        }
        return Arrays.copyOf(out, w);
    }

}
