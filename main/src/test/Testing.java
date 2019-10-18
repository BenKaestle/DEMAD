package test;

import utility.FastaParser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class Testing {
    public static void main(String[] args) throws IOException {
        
        ArrayList<String> x = new ArrayList<>();
        x.add("b");
        x.add("c");
        x.add("a");
        System.out.println(x.remove(0));
        System.out.println(x.remove(0));
        System.out.println(x.remove(0));
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
