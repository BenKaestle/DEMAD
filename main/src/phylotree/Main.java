package phylotree;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class Main {

    public static void main (String[] args){
        System.out.println(ParseRefSeq.getAllCompleteGenomes(args[1]).size());
        downloadFromYear(Integer.parseInt(args[0]), ParseRefSeq.getAllCompleteGenomes(args[1]));
    }

    public static void downloadFromYear(int year, ArrayList<CompleteGenome> genomes){
        Process p;
        for (CompleteGenome genome : genomes){
            if (genome.year == year){
                try {
                    p = Runtime.getRuntime().exec("wget "+genome.downloadLink);
                    p.waitFor();
                    p = Runtime.getRuntime().exec("gunzip "+genome.downloadLink);
                    p.waitFor();
                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

    }
}
