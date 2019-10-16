package mash;

import hash_functions.Rabin_Karp;
import hash_functions.Test_Hash;
import utility.FastaParser;

import java.util.Arrays;
import java.util.concurrent.*;

class Task implements Callable<long[]>
{
    private String name;
    private CountDownLatch latch;
    private int sketchSize;
    private int kmerSize;
    private String sequence;
    private long hashModulo;

    public Task(String name, CountDownLatch latch, int sketchSize, int kmerSize, String sequence, long hashModulo) {
        this.name = name;
        this.latch = latch;
        this.sketchSize = sketchSize;
        this.kmerSize = kmerSize;
        this.sequence = sequence;
        this.hashModulo = hashModulo;
    }

    public String getName() {
        return name;
    }

    //todo is it faster to just reverse the whole sequence and take the smaller one?
    public static String canonical_kmer (String kmer){
        String reverse = "";
        for (int i=kmer.length()-1; i>=0;i--){
            switch (kmer.charAt(i)){
                case 'A':
                    reverse+="T";
                    break;
                case 'T':
                    reverse+="A";
                    break;
                case 'G':
                    reverse+="C";
                    break;
                case 'C':
                    reverse+="G";
                    break;
            }
        }
        return kmer.compareTo(reverse) > 0 ? reverse : kmer;
    }

    @Override
    public long[] call() {
        long[] sketch = new long[sketchSize];
        Arrays.fill(sketch, Integer.MAX_VALUE);
        for (int i = 0; i <= sequence.length() - kmerSize; i++) {
            Test_Hash test_hash = new Test_Hash();
            long hash = test_hash.hash(canonical_kmer(sequence.substring(i, i + kmerSize)),hashModulo);
            if (hash < sketch[0]) {
                boolean found = false;
                for (int j = 1; j < sketchSize; j++) {
                    if (hash < sketch[j]) {
                        sketch[j - 1] = sketch[j];
                    } else {
                        sketch[j - 1] = hash;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    sketch[sketchSize - 1] = hash;
                }
            }
        }
        latch.countDown();
        return sketch;
//            Long duration = (long) (Math.random() * 50);
//            System.out.println("Doing a task during : " + name);
//            TimeUnit.SECONDS.sleep(duration);


    }
}



public class Mash {
    public static void main(String[] args){
//        Rabin_Karp rabin_karp = new Rabin_Karp(11,0);
//        System.out.println(rabin_karp.hash("ATCAT", (long) Math.pow(2,32)));
//        System.out.println(rabin_karp.hash("ATGTGCGATTCGATTCGATTA", (long) Math.pow(2,32)));
//        System.out.println(rabin_karp.hash("GGGGGGGGGGGGGGGGGGGGG", (long) Math.pow(2,32)));
//        System.out.println( (long) Math.pow(2,64));

        String[] x = new String[]{FastaParser.getExample()};
        mash_sketch(x,100,15,4, (long)Math.pow(2,32));
    }
    public static int[] mash_sketch (String[] sequences, int sketchSize, int kmerSize, int cores, long hashModulo){
        String[][] split = new String[sequences.length][cores];
        int i=0;
        for (String sequence : sequences){
            split[i] = splitBy(cores, sequence, kmerSize);
            i++;
        }

        runInParallel(cores,sketchSize,kmerSize,split, hashModulo);


        return null;
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

    public static void runInParallel(int cores, int sketchSize, int kmerSize, String[][] sequences, long hashModulo){
        int threads=sequences.length*cores;
        ExecutorService pool = Executors.newFixedThreadPool(cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<long[]>[][] results = new Future[sequences.length][cores];
        for (int i = 0; i < threads; i++)
        {
            Task task = new Task("Task " + i, latch,sketchSize, kmerSize, sequences[i/cores][i%cores], hashModulo);
            System.out.println("A new task has been added : " + task.getName());
            results[i/cores][i%cores] = pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        pool.shutdown();
        long[][] combined = new long[cores][sketchSize];
        long[][] sketches = new long[sequences.length][sketchSize];
        int seq_num=0;
        for (Future<long[]>[] x : results){
            int task_num=0;
            for (Future<long[]> xx : x) {
                try {
                    combined[task_num] = xx.get();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    e.printStackTrace();
                }
                task_num++;
            }
            sketches[seq_num] = combine(combined);
            seq_num++;
        }


    }

    private static long[] combine(long[][] combined) {
        int cores = combined.length;
        int sketch_size = combined[0].length;
        int[] markers = new int[cores];
        Arrays.fill(markers,sketch_size-1);
        long[] result = new long[sketch_size];
        int result_marker=0;
        while (result_marker<sketch_size){
            long min = Long.MAX_VALUE;
            int[] min_args = new int[cores];
            for(int i=0; i<cores;i++){
                if(min > combined[i][markers[i]]){
                    min_args = new int[cores];
                    min_args[i] = 1;
                    min = combined[i][markers[i]];
                }
                if(min == combined[i][markers[i]]) {
                    min_args[i] = 1;
                }
            }
            for(int i=0; i<cores;i++){
                if(min_args[i] == 1){
                    markers[i]--;
                }
            }
            result[result_marker] = min;
            result_marker++;
        }
        return result;
    }


}
