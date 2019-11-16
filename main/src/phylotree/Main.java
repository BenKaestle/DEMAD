package phylotree;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class Main {

    public static void main (String[] args){
        System.out.println(ParseRefSeq.getAllCompleteGenomes(args[1]).size());
        downloadFromYear(Integer.parseInt(args[0]), ParseRefSeq.getAllCompleteGenomes(args[1]));
    }

    private static void applyMash(String filepath, String mashfilepath){
        ArrayList<String> files = getAllFilesFromDir(filepath);
        String command = mashfilepath +" sketch ";
        for (String s : files){
            command.concat(s+" ");
        }
        command.concat("-o "+filepath+"/sketchout -p 100");
    }

    public static ArrayList<String> getAllFilesFromDir (String filepath){
        final File folder = new File(filepath);

        ArrayList<String> result = new ArrayList<>();

        search(".*\\.fna", folder, result);
        Collections.sort(result);
        return result;
    }

    public static void search(final String pattern, final File folder, ArrayList<String> result) {
        for (final File f : folder.listFiles()) {

            if (f.isDirectory()) {
                search(pattern, f, result);
            }

            if (f.isFile()) {
                if (f.getName().matches(pattern)) {
                    result.add(f.getAbsolutePath());
                }
            }

        }
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
