package utility;

import java.io.*;

public class FastaParser {
    public static String[] getExample (){
        String path = "main/resources/test.fasta";
        File file = new File(path);
        try {
            return parseFasta(file);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
    public static String[] getBigExample (){
        String path = "main/resources/sequence.fasta";
        File file = new File(path);
        try {
            return parseFasta(file);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }
    public static String[] parseFasta(File inputFile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String st;
        String header = "";
        StringBuilder fasta=new StringBuilder();
        while ((st = br.readLine()) != null) {
            if (!st.startsWith(">"))
                fasta.append(st.trim());
            else
                header = st;
        }
        return new String[]{header, fasta.toString()};
    }
}

