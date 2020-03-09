package utility;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/*
 *  FastaParser.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 *
 * Benjamin Kaestle, 3.2020
 */

public class FastaParser {
    private static Map<String, String[]> reduce = new HashMap<String, String[]>() {{
        put("L", new String[]{"L", "L", "L", "L", "L"});
        put("V", new String[]{"L", "L", "L", "L", "L"});
        put("I", new String[]{"L", "L", "L", "L", "L"});
        put("M", new String[]{"L", "L", "L", "L", "L"});
        put("C", new String[]{"C", "C", "L", "L", "L"});
        put("A", new String[]{"A", "A", "A", "A", "L"});
        put("G", new String[]{"G", "G", "A", "A", "L"});
        put("S", new String[]{"S", "S", "S", "A", "L"});
        put("T", new String[]{"T", "S", "S", "A", "L"});
        put("P", new String[]{"P", "P", "P", "A", "L"});
        put("F", new String[]{"F", "F", "F", "F", "L"});
        put("Y", new String[]{"F", "F", "F", "F", "L"});
        put("W", new String[]{"W", "F", "F", "F", "L"});
        put("E", new String[]{"E", "E", "E", "E", "E"});
        put("D", new String[]{"D", "E", "E", "E", "E"});
        put("N", new String[]{"N", "E", "E", "E", "E"});
        put("Q", new String[]{"Q", "E", "E", "E", "E"});
        put("K", new String[]{"K", "K", "K", "E", "E"});
        put("R", new String[]{"K", "K", "K", "E", "E"});
        put("H", new String[]{"H", "H", "H", "E", "E"});
    }};

    public static String[] getExample() {
        String path = "main/resources/test.fasta";
        File file = new File(path);
        try {
            return parseFasta(file, "DNA", -1);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static String[] getBigExample() {
        String path = "main/resources/sequence.fasta";
        File file = new File(path);
        try {
            return parseFasta(file, "DNA", -1);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     *
     * @param inputFile
     * @param alphabet
     * @param redu - -1:no reduction 0: reduction to 15; 1: reduction to 10; 2: reduction to 8; 3: reduction to 4; 4: reduction to 2
     * @return
     * @throws IOException
     */
    public static String[] parseFasta(File inputFile, String alphabet, int redu) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String st;
        String header = "";
        StringBuilder fasta = new StringBuilder();
        StringBuilder reverse = new StringBuilder();
        int sequences = 1;
        boolean nonAlphabetical = false;
        while ((st = br.readLine()) != null) {
            if (!st.startsWith(">")) {
                if (alphabet.equals("AA") && redu != -1) {
                    String y = "";
                    for (Character c : st.toCharArray()) {
                        y += reduce.get(c.toString())[redu];
                    }
                    st = y;
                }
                st = st.trim();
                fasta.append(st);
                if (!consitsOfAlphabet(st.trim(), alphabet)) {
                    nonAlphabetical = true;
                }
                if (alphabet.equals("DNA")) {
                    for (char c : st.toCharArray()) {
                        switch (c) {
                            case 'A':
                                reverse.append("T");
                                break;
                            case 'T':
                                reverse.append("A");

                                break;
                            case 'G':
                                reverse.append("C");

                                break;
                            case 'C':
                                reverse.append("G");

                                break;
                            default:
                                reverse.append("N");
                                break;
                        }
                    }
                }
            } else if (header == "") {
                header = st;
            } else {
                sequences++;
            }
        }
        if (sequences != 1) {
            System.out.println("WARNING: input file \"" + inputFile.getPath() + "\" contains " + sequences + " sequences");
        }
        if (nonAlphabetical) {
            System.out.println("WARNING: input file \"" + inputFile.getPath() + "\" contains letters other than the alphabet of " + alphabet);
        }
        return new String[]{header, fasta.toString(), reverse.reverse().toString()};
    }

    /**
     * checks if a Sequence consist of a given alphabet or not (DNA includes 'N' as well)
     *
     * @param s        - sequence that gets checked
     * @param alphabet - one of DNA ode AA (null returns true)
     * @return - boolean if sequence consists of given alphabet
     */
    public static boolean consitsOfAlphabet(String s, String alphabet) {
        if (alphabet == null) return true;
        if (s.equals("")) return true;
        switch (alphabet) {
            case "DNA":
                return ((!s.equals(""))
                        && (s != null)
                        && (s.matches("[ATGCN]+")));
            case "AA":
                return ((!s.equals(""))
                        && (s != null)
                        && (s.matches("[GALMFWKQESPVICYHRNDT]+")));
            default:
                return false;
        }
    }
}

