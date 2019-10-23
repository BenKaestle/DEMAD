package mash;

import hash_functions.Murmur3;
import hash_functions.Rabin_Karp;
import hash_functions.Test_Hash;
import utility.*;

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
final class DistTask implements Callable<ArrayList<Distance>>
{

    private String name;
    private Sketch sketch_1;
    private Sketch sketch_2;
    private Sketch[] sketches;
    private long current_hash_1;
    private long current_hash_2;
    private int sketch_length;
    private int same_hash_counter;
    private ArrayList<Distance> distances;
    private ArrayList<Sketch[]> pairOfSketches;
    private float jaccard_index, p_value, mash_distance;
    private int kmer_length;
    private CountDownLatch latch;
    private float r, pkx, pky;
    private final int ALPHABET_SIZE = 4;
    private int genome_size_1, genome_size_2;

    public DistTask(String name, ArrayList<Sketch[]> pairOfSketches, CountDownLatch latch, int kmer_length) {
        this.name = name;
        this.pairOfSketches = pairOfSketches;
        this.latch = latch;
        this.kmer_length = kmer_length;
    }

    public String getName() {
        return name;
    }

    @Override
    public synchronized ArrayList<Distance> call() throws Exception {
        distances = new ArrayList<>();
        int pointer_1, pointer_2;
        while (!pairOfSketches.isEmpty()) {
            pointer_1=0;
            pointer_2=0;
            sketches = pairOfSketches.remove(0);
            sketch_1=sketches[0];
            sketch_2=sketches[1];
            sketch_length = sketch_1.getHashes().length;
            same_hash_counter=0;
            for (int i =0; i<sketch_length ;i++){
                current_hash_1 = sketch_1.getHashes()[pointer_1];
                current_hash_2 = sketch_2.getHashes()[pointer_2];
                if (current_hash_1>current_hash_2){
                    pointer_1++;
                }
                if (current_hash_1<current_hash_2){
                    pointer_2++;
                }
                if (current_hash_1==current_hash_2){
                    same_hash_counter++;
                    pointer_1++;
                    pointer_2++;
                }
            }
            System.out.println(same_hash_counter);
            jaccard_index = (float) same_hash_counter/sketch_length; //todo jaccard estimate = jaccard index???
            mash_distance = -1f/kmer_length * (float) Math.log((2*jaccard_index)/(1+jaccard_index));
            genome_size_1 = sketch_1.getGenome_size();
            genome_size_2 = sketch_2.getGenome_size();
            pkx = 1- (float)Math.pow(1-Math.pow(ALPHABET_SIZE,-kmer_length),genome_size_1);
            pky = 1- (float)Math.pow(1-Math.pow(ALPHABET_SIZE,-kmer_length),genome_size_2);
            r = pkx*pky/(pkx+pky-pkx*pky);
            p_value=1;
            System.out.println(r);
            for (int i =0; i<same_hash_counter;i++){
                System.out.println(p_value);
                p_value -= (float)(binom(sketch_length,i)*Math.pow(r,i)*Math.pow((1-r),(sketch_length-i)));
            }
            distances.add(new Distance(sketch_1.getHeader(), sketch_2.getHeader(), jaccard_index, p_value, mash_distance));
            latch.countDown();
        }
        return distances;
    }
    public static double binom(int n, int k) {
        double res = 1;
        while (k > 0) {
            res = res * ((double) n / (double) k);
            k--;
            n--;
        }
        return res;
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
    private BloomFilter bloomFilter;
    private SketchSet sketchSet;
    private String kmer;

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
                default:
                    reverse+="N";
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
            sketch = new Sketch(sketchHashes, sequenceInfo[0], sequence.length());
            Arrays.fill(sketchHashes, Integer.MAX_VALUE);
            long hash = 0;
            bloomFilter = new BloomFilter(16,4);
            sketchSet = new SketchSet(5);
            for (int i = 0; i <= sequence.length() - kmerSize; i++) {

//                switch (hashType) {
//                    case "Test":
//                        Test_Hash test_hash = new Test_Hash();
//                        hash = test_hash.hash(canonical_kmer(sequence.substring(i, i + kmerSize)), hashModulo);
//                        break;
//                    case "Rabin_Karp":
//                        if (hash == 0)
//                            hash = Rabin_Karp.hash(canonical_kmer(sequence.substring(i, i + kmerSize)), hashModulo, 11, 5);
//                        else
//                            hash = Rabin_Karp.hashNext(sequence.charAt(i - 1), sequence.charAt(i + kmerSize), hash, hashModulo, 11, 5, kmerSize);
//                        break;
//                    case "Murmur3":
//                        hash = Murmur3.murmurhash3_x86_32(canonical_kmer(sequence.substring(i, i + kmerSize)), 0, kmerSize, 1);
//                        break;
//
//
//                }
                kmer = canonical_kmer(sequence.substring(i, i + kmerSize));
                hash = Murmur3.murmurhash3_x86_32(kmer, 0, kmerSize, 1);
                count++;

                if (hash < sketchHashes[0]) {
                    if(sketchSet.contains(hash)) {
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
                    else{
                        sketchSet.add(hash);
                    }
                }
            }
            System.out.println("kmers hashed: " + count + " by thread " + this.name);
            latch.countDown();
            sketches.add(sketch);
        }
        return sketches;
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
//        for (int i =0;i<2;i++)
        x.add("main/resources/sequence2.fasta");
        x.add("main/resources/sequence.fasta");

        int sketchSize = 1000;
        int kmerSize = 10;
        int cores = 4;


        long startTime = System.currentTimeMillis();
        mash_dist(mash_sketch(x,sketchSize,kmerSize,cores, (long)Math.pow(2,32), "Murmur3"),cores, kmerSize)[0].print();
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        System.out.println(duration);

    }
    private static Sketch[] mash_sketch(ArrayList<String> sequences, int sketchSize, int kmerSize, int cores, long hashModulo, String hashType){


        return mash_sketch_in_parallel(cores,sketchSize,kmerSize,sequences,hashModulo,hashType);
    }

    private static Sketch[] mash_sketch_in_parallel(int cores, int sketchSize, int kmerSize, ArrayList<String> sequences, long hashModulo, String hashType){
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
        return sketches;
    }

    private static Distance[] mash_dist (Sketch[] sketches, int cores, int kmer_length){
        int length = sketches.length;
        ArrayList<Sketch[]> pairOfSketches = new ArrayList<>();
        for (int i =0; i<length;i++){
            for (int j =i+1; j<length;j++){
                pairOfSketches.add(new Sketch[]{sketches[i],sketches[j]});
            }

        }
        return mash_dist_in_parallel(pairOfSketches, cores, kmer_length);


    }

    private static Distance[] mash_dist_in_parallel(ArrayList<Sketch[]> pairOfSketches, int cores, int kmer_length) {
        int threads=pairOfSketches.size();
        ExecutorService pool = Executors.newFixedThreadPool(cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<Distance>>[] results = new Future[cores];
        for (int i = 0; i < cores; i++)
        {
            DistTask task = new DistTask("Task " + i, pairOfSketches, latch, kmer_length);
            System.out.println("A new task has been added : " + task.getName());
            results[i] = pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        pool.shutdown();
        Distance[] distances = new Distance[threads];
        int pointer = 0;
        for(Future<ArrayList<Distance>> distance : results){
            try {
                for(Distance d : distance.get()){
                    distances[pointer] = d;
                    pointer++;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return distances;
    }

    /**
     * for a given genome size and the desired probability prob of observing a random k-mer
     * @param genome_size
     * @param prob - probability of observing a random k-mer
     * @return
     */
    public static int optimalK(long genome_size, float prob){
        return (int) Math.ceil(Math.log(genome_size*(1-prob)/prob)/Math.log(4));
    }




}
