package test;

import utility.BloomFilter;
import utility.FastaParser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class Testing {
    public static void main(String[] args) throws IOException {
        System.out.println(optimalK(3000000000l,0.01f));
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
