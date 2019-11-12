package phylotree;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class ParseRefSeq {

    public static ArrayList<CompleteGenome> getAllCompleteGenomes(String filepath) {
        BufferedReader reader;
        String[] split;
        ArrayList<CompleteGenome> result = new ArrayList<>();
        try {
            reader = new BufferedReader(new FileReader(filepath));
            String line = reader.readLine();
            while (line != null) {
                if (!line.startsWith("#")) {
                    split = line.split("\t");
                    if(split[11].trim().equals("Complete Genome")){
                        CompleteGenome x = new CompleteGenome(split[19]+"/"+split[19].split("/")[split[19].split("/").length-1]+"_genomic.fna.gz", Integer.parseInt(split[14].split("/")[0]), split[0], Integer.parseInt(split[5]));
                        result.add(x);
                    }

                }
                line = reader.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return result;
    }
}
