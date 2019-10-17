package mash;

import hash_functions.Murmur3;
import hash_functions.Rabin_Karp;
import hash_functions.Test_Hash;
import utility.FastaParser;

import java.util.Arrays;
import java.util.concurrent.*;

class MergeTask implements Callable<long[]>
{
    private long[][] combined;
    private CountDownLatch latch;
    private String name;

    public MergeTask(String name,long[][] combined, CountDownLatch latch) {
        this.combined = combined;
        this.latch = latch;
        this.name = name;
    }

    public String getName() {
        return name;
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

    @Override
    public long[] call() throws Exception {
        latch.countDown();
        return combine(combined);
    }
}
class KmerTask implements Callable<long[]>
{
    private String name;
    private CountDownLatch latch;
    private int sketchSize;
    private int kmerSize;
    private String sequence;
    private long hashModulo;
    private String hashType;

    public KmerTask(String name, CountDownLatch latch, int sketchSize, int kmerSize, String sequence, long hashModulo, String hashType) {
        this.name = name;
        this.latch = latch;
        this.sketchSize = sketchSize;
        this.kmerSize = kmerSize;
        this.sequence = sequence;
        this.hashModulo = hashModulo;
        this.hashType = hashType;
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
        int count = 0;
        long[] sketch = new long[sketchSize];
        Arrays.fill(sketch, Integer.MAX_VALUE);
        for (int i = 0; i <= sequence.length() - kmerSize; i++) {

            long hash=0;
            switch (hashType){
                case "Test":
                    Test_Hash test_hash = new Test_Hash();
                    hash = test_hash.hash(canonical_kmer(sequence.substring(i, i + kmerSize)),hashModulo);
                    break;
                case "Rabin_Karp":
                    Rabin_Karp rabin_karp = new Rabin_Karp(11,5);
                    hash = rabin_karp.hash(canonical_kmer(sequence.substring(i, i + kmerSize)),hashModulo);
                    break;
                case "Murmur3":
                    Murmur3 murmur3 = new Murmur3();
                    hash = murmur3.hash(canonical_kmer(sequence.substring(i, i + kmerSize)),hashModulo);
                    break;


            };
            count++;
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
        System.out.println("kmers hashed: " + count);
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

        String[] x = new String[]{FastaParser.getBigExample()};
        System.out.println("eingelesen");
        mash_sketch(x,1000,21,4, (long)Math.pow(2,32), "Rabin_Karp");
    }
    public static int[] mash_sketch(String[] sequences, int sketchSize, int kmerSize, int cores, long hashModulo, String hashType){
        String[][] split = new String[sequences.length][cores];
        int i=0;
        for (String sequence : sequences){
            split[i] = splitBy(cores, sequence, kmerSize);
            i++;
        }
        System.out.println("splitted");

        runInParallel(cores,sketchSize,kmerSize,split, hashModulo, hashType);


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

    public static void runInParallel(int cores, int sketchSize, int kmerSize, String[][] sequences, long hashModulo, String hashType){
        int threads=sequences.length*cores;
        ExecutorService pool = Executors.newFixedThreadPool(cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<long[]>[][] results = new Future[sequences.length][cores];
        for (int i = 0; i < threads; i++)
        {
            KmerTask task = new KmerTask("Task " + i, latch,sketchSize, kmerSize, sequences[i/cores][i%cores], hashModulo, hashType);
            System.out.println("A new task has been added : " + task.getName());
            results[i/cores][i%cores] = pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        long[][] combined = new long[cores][sketchSize];
        Future<long[]>[] sketches = new Future[sequences.length];
        threads = sequences.length;
        latch = new CountDownLatch(threads);
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
            MergeTask task = new MergeTask("Task " + seq_num,combined, latch);
            System.out.println("A new task has been added : " + task.getName());
            sketches[seq_num] = pool.submit(task);
            seq_num++;
        }

        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        pool.shutdown();
        for (Future<long[]> x :sketches){
            long[] y=null;
            try {
                y = x.get();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
            for (long z : y){
                System.out.println(z);
            }
        }
    }


}
