package test;

import utility.BloomFilter;
import utility.FastaParser;
import utility.Sketch;
import utility.WriteReadObject;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Testing {
    public static void main(String[] args) throws IOException {
//        Sketch sketch = new Sketch(new long[3], "header",3);
//        Sketch ssketch = new Sketch(new long[4], "header2",4);
//        Sketch[] sketches = new Sketch[]{sketch,ssketch};
//        WriteReadObject.writeObjectToFile(sketches,"main/resources/test");
//        Sketch[] s = (Sketch[])WriteReadObject.readObjectFromFile("main/resources/test");
//        System.out.println(s[0].getGenome_size());
//        System.out.println(s[1].getGenome_size());
        String x="";
        printTable(WriteReadObject.readTSVtable(x));
        WriteReadObject.writePhylipFile(WriteReadObject.readTSVtable(x),"original_mash_result");
    }

    //01.fna 02.fna 03.fna 04.fna 05.fna 06.fna 07.fna 08.fna 09.fna 10.fna 11.fna 12.fna 13.fna 14.fna 15.fna 16.fna 17.fna 18.fna 19.fna 20.fna 21.fna 22.fna 23.fna 24.fna 25.fna 26.fna 27.fna 28.fna 29.fna 30.fna 31.fna 32.fna 33.fna 34.fna 35.fna 36.fna 37.fna 38.fna 39.fna 40.fna 41.fna 42.fna 43.fna 44.fna 45.fna 46.fna 47.fna 48.fna 49.fna 50.fna

    private static void printTable(String[][] strings) {
        for (int i=0;i<strings.length;i++){
            for (int j=0;j<strings[0].length;j++){
                System.out.print(strings[i][j]+"\t");
            }
            System.out.println();
        }
    }

    public static void usingFileWriter() throws IOException
    {
        String fileContent = "4\n" +
                "LS1         0  0.083  0.25  0.458\n" +
                "LS2         0.083  0  0.167  0.392\n" +
                "LS3         0.25  0.167  0  0.392\n" +
                "LS4         0.458  0.392  0.392  0";

        FileWriter fileWriter = new FileWriter("main/test.phylip");
        fileWriter.write(fileContent);
        fileWriter.close();
    }


    public static int optimalK(long genome_size, float prob){
        return (int) Math.ceil(Math.log(genome_size*(1-prob)/prob)/Math.log(4));
    }

    public static long[][] getBinomialDistribution(int min, int max, long total) {
        int n = max - min;
        long[][] ret = new long[2][n + 1];
        int mean = (n + 1) / 2;
        float p = 1;
        if (n > 0) {
            p = (float) mean / (float) n;
        }

        long count = 0;
        for (int i = 0; i <= n; i++) {
            double p_i = combination(n, i) * Math.pow(p, i)
                    * Math.pow((1 - p), (n - i));
            long count_i = (long) (total * p_i);
            ret[0][i] = i + min;
            ret[1][i] = count_i;
            count += count_i;
        }

        while (count < total) {
            int i = 1;
            ret[1][i]++;
            count++;
        }

        return ret;
    }
    // calculate the combination
    // the value would be very large, so store it in the type of double
    public static double combination(int n, int k) {
        double ret = 1;
        while (k > 0) {
            ret = ret * ((double) n / (double) k);
            k--;
            n--;
        }
        return ret;
    }



    public static void test(){
        String str = "abcdbdbcbcbbcadbcda";
        long y = 0;
        for (int i =0; i<21; i++){

            long x = (long)4*(long)Math.pow(5,49-i)%(long)Math.pow(2,32);
            y+=x;
            System.out.println(y);
        }
    }

    public static String[] splitBy(int cores, String sequence, int kmersize){
        String[] split = new String[cores];
        float splitHere= sequence.length()/(float)cores;
        int addDown = (kmersize-1)/2;
        int addUp = (int) Math.ceil((kmersize-1)/2);

        for (int i=0;i<cores;i++){
            int start=(int)(splitHere*i)-addDown;
            int end =(int)(splitHere*(i+1))+addUp;
            if (start<0)
                start=0;
            if (end>sequence.length())
                end = sequence.length();
            split[i] = sequence.substring(start,end);
        }
        return split;
    }


}
