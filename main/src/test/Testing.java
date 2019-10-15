package test;

import utility.FastaParser;

import java.io.File;
import java.io.IOException;

public class Testing {
    public static void main(String[] args) throws IOException {
        String path = "main/resources/test.fasta";
        File file = new File(path);
        System.out.println(FastaParser.parseFasta(file));
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
}
