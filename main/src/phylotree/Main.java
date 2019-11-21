package phylotree;

import mash.MashDistance;
import sun.tools.jar.CommandLine;
import utility.WriteReadObject;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

public class Main {

    public static void main(String[] args) throws IOException {
//        System.out.println(ParseRefSeq.getAllCompleteGenomes(args[1]).size());
//        downloadFromYear(Integer.parseInt(args[0]), ParseRefSeq.getAllCompleteGenomes(args[1]));

//        applyMash(args[0]);
        combineLists(args);
    }

    private static void combineLists(String[] args) throws FileNotFoundException {
        FileWriter fileWriter = null;
        ArrayList<BufferedReader> listOfReaders = new ArrayList<>();
        String line;
        String result = "";
        String subresult = "";
        boolean next = true;
        String[] split;
        for (int i = 1; i < args.length; i++) {
            listOfReaders.add(new BufferedReader(new FileReader(args[i])));
        }
        try {
            fileWriter = new FileWriter(args[0] + ".txt");
            while (next) {
                for (BufferedReader r : listOfReaders) {
                    line = r.readLine();
                    if (line != null) {
                        split = line.split("\t");
                        if (split.length == 1) subresult = line;
                        else if (split.length == 5) {
                            subresult = split[2] + "\t" + split[3] + "\t" + split[4].split("/")[0];
                        } else {
                            System.out.println("wtf");
                        }
                        result += subresult + "\t";
                    }
                }
                if (result == "")
                    next = false;
                else
                    fileWriter.write(result + "\n");
                result = "";
            }

            fileWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    private static void constructPhyloList(String[] args) {
        ArrayList<String> result = new ArrayList<>();
        ArrayList<CompleteGenome> x = getFromYear(2019, ParseRefSeq.getAllCompleteGenomes(args[0]));
        for (CompleteGenome c : x) {
            for (CompleteGenome v : x) {
                result.add(v.taxId + "\t" + v.taxId2 + "\t" + c.taxId + "\t" + c.taxId2);
            }
        }
        WriteReadObject.writeTxt("out_test", result);

    }

    private static void applyMash(String filepath) throws IOException {
        Process p = null;
        try {
            p = Runtime.getRuntime().exec("/home/kaestle/Mash/mash-Linux64-v2.2/mash dist " + filepath + " " + filepath + " -p 100\n");
            BufferedReader buf = new BufferedReader(new InputStreamReader(
                    p.getInputStream()));
            String line = "";
            FileWriter fileWriter = null;
            fileWriter = new FileWriter(filepath + "_out.txt");
            int i = 0;
            while ((line = buf.readLine()) != null) {
                i++;
                if (i % 1000 == 0) System.out.println(i);
                fileWriter.write(line + "\n");
            }
            fileWriter.close();


        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public static ArrayList<String> getAllFilesFromDir(String filepath) {
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

    public static void downloadFromYear(int year, ArrayList<CompleteGenome> genomes) {
        Process p, q;
        int count = 0;
        for (CompleteGenome genome : genomes) {
            if (genome.year == year) {
                count++;
            }
        }
        int counter = 0;
        for (CompleteGenome genome : genomes) {
            if (genome.year == year) {
                try {
                    counter++;
                    System.out.println(counter + " / " + count);
                    p = Runtime.getRuntime().exec("wget " + genome.downloadLink);
                    p.waitFor();
                    q = Runtime.getRuntime().exec("gunzip " + genome.downloadLink);
                    q.waitFor();
                } catch (IOException | InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public static ArrayList<CompleteGenome> getFromYear(int year, ArrayList<CompleteGenome> genomes) {
        ArrayList<CompleteGenome> result = new ArrayList<>();
        for (CompleteGenome genome : genomes) {
            if (genome.year == year) {
                result.add(genome);
            }
        }
        Collections.sort(result, Comparator.comparing(a -> a.name));
        return result;
    }
}
