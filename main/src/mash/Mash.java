package mash;

import hash_functions.Murmur3;
import hash_functions.Rabin_Karp;
import hash_functions.Test_Hash;
import utility.FastaParser;
import utility.Sketch;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
final class KmerTask implements Callable<ArrayList<Sketch>>
{
    private String name;
    private CountDownLatch latch;
    private int sketchSize;
    private int kmerSize;
    private long hashModulo;
    private String hashType;
    private ArrayList<String> listOfSequences;
    private Sketch sketch;
    private ArrayList<Sketch> sketches;

    public KmerTask(String name, CountDownLatch latch, int sketchSize, int kmerSize, long hashModulo, String hashType, ArrayList<String> listOfSequences) {
        this.name = name;
        this.latch = latch;
        this.sketchSize = sketchSize;
        this.kmerSize = kmerSize;
        this.hashModulo = hashModulo;
        this.hashType = hashType;
        this.listOfSequences = listOfSequences;
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
    public synchronized ArrayList<Sketch> call() {
        sketches = new ArrayList<>();
        while (!listOfSequences.isEmpty()) {
            String path = listOfSequences.remove(0);

            String[] sequenceInfo = new String[2];
            try {
                sequenceInfo = FastaParser.parseFasta(new File(path));
            } catch (IOException e) {
                e.printStackTrace();
            }
            int count = 0;
            String sequence = sequenceInfo[1];
            long[] sketchHashes = new long[sketchSize];
            sketch = new Sketch(sketchHashes, sequenceInfo[0]);
            Arrays.fill(sketchHashes, Integer.MAX_VALUE);
            long hash = 0;
            for (int i = 0; i <= sequence.length() - kmerSize; i++) {

                switch (hashType) {
                    case "Test":
                        Test_Hash test_hash = new Test_Hash();
                        hash = test_hash.hash(canonical_kmer(sequence.substring(i, i + kmerSize)), hashModulo);
                        break;
                    case "Rabin_Karp":
                        if (hash == 0)
                            hash = Rabin_Karp.hash(canonical_kmer(sequence.substring(i, i + kmerSize)), hashModulo, 11, 5);
                        else
                            hash = Rabin_Karp.hashNext(sequence.charAt(i - 1), sequence.charAt(i + kmerSize), hash, hashModulo, 11, 5, kmerSize);
                        break;
                    case "Murmur3":
                        hash = Murmur3.murmurhash3_x86_32(canonical_kmer(sequence.substring(i, i + kmerSize)), 0, kmerSize, 1);
                        break;


                }
                ;
                count++;
                if (hash < sketchHashes[0]) {
                    boolean found = false;
                    for (int j = 1; j < sketchSize; j++) {
                        if (hash < sketchHashes[j]) {
                            sketchHashes[j - 1] = sketchHashes[j];
                        } else {
                            sketchHashes[j - 1] = hash;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        sketchHashes[sketchSize - 1] = hash;
                    }
                }
            }
            System.out.println("kmers hashed: " + count + " by thread " + this.name);
            latch.countDown();
            sketches.add(sketch);
        }
        return sketches;
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

        ArrayList<String> x = new ArrayList<>();
        for (int i =0;i<20;i++)
            x.add("main/resources/sequence.fasta");
        long startTime = System.currentTimeMillis();
        mash_sketch(x,1000,21,4, (long)Math.pow(2,32), "Murmur3");
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        System.out.println(duration);

    }
    public static Sketch[] mash_sketch(ArrayList<String> sequences, int sketchSize, int kmerSize, int cores, long hashModulo, String hashType){


        runInParallel(cores,sketchSize,kmerSize,sequences,hashModulo,hashType);


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

    public static void runInParallel(int cores, int sketchSize, int kmerSize, ArrayList<String> sequences, long hashModulo, String hashType){
        int threads=sequences.size();
        ExecutorService pool = Executors.newFixedThreadPool(cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<Sketch>>[] results = new Future[cores];
        for (int i = 0; i < cores; i++)
        {
            KmerTask task = new KmerTask("Task " + i, latch,sketchSize, kmerSize, hashModulo, hashType, sequences);
            System.out.println("A new task has been added : " + task.getName());
            results[i] = pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        pool.shutdown();
        Sketch[] sketches = new Sketch[threads];
        int pointer = 0;
        for(Future<ArrayList<Sketch>> sketch : results){
            try {
                for(Sketch s : sketch.get()){
                    sketches[pointer] = s;
                    pointer++;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }



        for (Sketch s : sketches){
            for (long l : s.getHashes()){
                System.out.println(l);
            }
        }
//        long[][] combined = new long[cores][sketchSize];
//        Future<long[]>[] sketches = new Future[sequences.length];
//        threads = sequences.length;
//        latch = new CountDownLatch(threads);
//        int seq_num=0;
//        for (Future<long[]>[] x : results){
//            int task_num=0;
//            for (Future<long[]> xx : x) {
//                try {
//                    combined[task_num] = xx.get();
//                } catch (InterruptedException e) {
//                    e.printStackTrace();
//                } catch (ExecutionException e) {
//                    e.printStackTrace();
//                }
//                task_num++;
//            }
//            MergeTask task = new MergeTask("Task " + seq_num,combined, latch);
//            System.out.println("A new task has been added : " + task.getName());
//            sketches[seq_num] = pool.submit(task);
//            seq_num++;
//        }
//
//        try {
//            latch.await();
//        } catch (InterruptedException e) {
//            e.printStackTrace();
//        }

//        pool.shutdown();
//        for (Future<long[]> x :sketches){
//            long[] y=null;
//            try {
//                y = x.get();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            } catch (ExecutionException e) {
//                e.printStackTrace();
//            }
//            for (long z : y){
//                System.out.println(z);
//            }
//        }
    }


}
