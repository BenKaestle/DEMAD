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
        Process p, q;
        int count=0;
        for (CompleteGenome genome : genomes) {
            if (genome.year == year) {
                count++;
            }
        }
        int counter=0;
        for (CompleteGenome genome : genomes){
            if (genome.year == year){
                try {
                    counter++;
                    System.out.println(counter + " / "+ count);
                    p = Runtime.getRuntime().exec("wget "+genome.downloadLink);
                    p.waitFor();
                    q = Runtime.getRuntime().exec("gunzip "+genome.downloadLink);
                    q.waitFor();
                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

    }
}
