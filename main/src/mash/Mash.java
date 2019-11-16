package mash;

import hash_functions.Murmur3;
import utility.*;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.*;
import org.apache.commons.cli.*;


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
    private float jaccard_index, mash_distance;
    private double p_value;
    private InputParameters parameters;
    private CountDownLatch latch;
    private float r, pkx, pky;
    private final int ALPHABET_SIZE = 4;
    private int genome_size_1, genome_size_2;

    public DistTask(String name, CountDownLatch latch, InputParameters parameters) {
        this.name = name;
        this.latch = latch;
        this.parameters = parameters;
    }

    public String getName() {
        return name;
    }

    @Override
    public ArrayList<Distance> call() throws Exception {
        distances = new ArrayList<>();
        int pointer_1, pointer_2;
        while (!parameters.pairsOfSketches.isEmpty()) {
            sketches = parameters.pairsOfSketches.get();
            if (sketches==null) break;
            sketch_1=sketches[0];
            sketch_2=sketches[1];
            pointer_1=sketch_1.getHashes().length-1;
            pointer_2=sketch_2.getHashes().length-1;
            sketch_length = Math.min(sketch_1.getHashes().length, sketch_2.getHashes().length);
            same_hash_counter=0;
            for (int i =0; i<sketch_length ;i++){
                current_hash_1 = sketch_1.getHashes()[pointer_1];
                current_hash_2 = sketch_2.getHashes()[pointer_2];
                if (current_hash_1>current_hash_2){
                    pointer_2--;
                }
                if (current_hash_1<current_hash_2){
                    pointer_1--;
                }
                if (current_hash_1==current_hash_2){
                    same_hash_counter++;
                    pointer_1--;
                    pointer_2--;
                }
            }
            jaccard_index = (float) same_hash_counter/sketch_length; //todo jaccard estimate = jaccard index???
//            if (jaccard_index==0f)
//                mash_distance=1;
//            else
//                mash_distance = -1f/parameters.kmerSize * (float) Math.log((2*jaccard_index)/(1+jaccard_index));
            mash_distance = (jaccard_index==0f)? 1:-1f/parameters.kmerSize * (float) Math.log((2*jaccard_index)/(1+jaccard_index));
            genome_size_1 = sketch_1.getGenome_size();
            genome_size_2 = sketch_2.getGenome_size();
            pkx = 1- (float)Math.pow(1-Math.pow(ALPHABET_SIZE,-parameters.kmerSize),genome_size_1);
            pky = 1- (float)Math.pow(1-Math.pow(ALPHABET_SIZE,-parameters.kmerSize),genome_size_2);
            r = pkx*pky/(pkx+pky-pkx*pky);
            p_value=1;
            for (int i =0; i<same_hash_counter;i++){
                p_value -= (binom(sketch_length,i)*Math.pow(r,i)*Math.pow((1-r),(sketch_length-i)));
            }
            distances.add(new Distance(sketch_1.getHeader(), sketch_2.getHeader(), jaccard_index, p_value, mash_distance, sketch_1.getFilename(), sketch_2.getFilename(), same_hash_counter));
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
    private Sketch sketch;
    private ArrayList<Sketch> sketches;
    private BloomFilter bloomFilter;
    private SketchSet sketchSet;
    private String kmer;
    private String reverseKmer;
    private InputParameters parameters;
    private int sequenceLength;

    public KmerTask(String name, CountDownLatch latch, InputParameters parameters) {
        this.name = name;
        this.latch = latch;
        this.parameters = parameters;
    }

    public String getName() {
        return name;
    }

    //todo is it faster to just reverse the whole sequence and take the smaller one?
    public static String canonical_kmer (String kmer){
        StringBuilder reverse = new StringBuilder();
        for (int i=kmer.length()-1; i>=0;i--){
            switch (kmer.charAt(i)){
                case 'A':
                    reverse.append("T");
                    break;
                case 'T':
                    reverse.append("A");

                    break;
                case 'G':
                    reverse.append("C");

                    break;
                case 'C':
                    reverse.append("G");

                    break;
                default:
                    reverse.append("N");
                    break;
            }
        }
        return kmer.compareTo(reverse.toString()) > 0 ? reverse.toString() : kmer;
    }

    /**
     * for a given genome size and the desired probability prob of observing a random k-mer calculate the optimal k value
     * @param genome_size
     * @param prob - probability of observing a random k-mer
     * @return
     */
    public static int optimalK(long genome_size, float prob){
        return (int) Math.ceil(Math.log(genome_size*(1-prob)/prob)/Math.log(4));
    }

    @Override
    public ArrayList<Sketch> call() throws Exception {
        sketches = new ArrayList<>();
        String[] sequenceInfo = new String[2];
        Murmur3.LongPair longPair = new Murmur3.LongPair();
        while (!parameters.sequences.isEmpty()) {
            long[] sketchHashes = new long[parameters.sketchSize];
            String path = parameters.sequences.get();
            if (path==null) break;
            try {
                sequenceInfo = FastaParser.parseFasta(new File(path));
            } catch (IOException e) {
                e.printStackTrace();
            }
            int count = 0;
            String sequence = sequenceInfo[1];
            String reverseSequence = sequenceInfo[2];
            if (optimalK(sequence.length(),parameters.lowKThreshold)>parameters.kmerSize){
                System.out.println("WARNING: For the k-mer size used ("+parameters.kmerSize+"), the random match probability is above the specified warning threshold ("+parameters.lowKThreshold+
                        ") for the sequence \""+path+"\" of size "+sequence.length()+". Distances to this sequence may be underestimated as a result. To meet the threshold of "+parameters.lowKThreshold+
                        ", a k-mer size of at least "+optimalK(sequence.length(),parameters.lowKThreshold)+" is required. See: -k, -w.");
            }
            Arrays.fill(sketchHashes, Long.MAX_VALUE);
            sketch = new Sketch(sketchHashes, sequenceInfo[0],path, parameters.kmerSize, parameters.hashFunction, sequence.length(), parameters.seed);
            long hash = 0;
            bloomFilter = new BloomFilter(parameters.bloomFilterSize,parameters.bloomFilterHashes);
            sketchSet = new SketchSet(5);
            sequenceLength = sequence.length();
            if (parameters.hashFunction==1) { //murmur3 x64_128
                for (int i = 0; i <= sequence.length() - parameters.kmerSize; i++) {
                    kmer = sequence.substring(i,i+parameters.kmerSize);
                    reverseKmer = reverseSequence.substring(sequenceLength-i-parameters.kmerSize, sequenceLength-i);
                    kmer = kmer.compareTo(reverseKmer) > 0 ? reverseKmer : kmer;
//                    kmer = canonical_kmer(sequence.substring(i, i + parameters.kmerSize));
                    try {
                        Murmur3.murmurhash3_x64_128(kmer.getBytes("UTF-8"), 0, parameters.kmerSize, parameters.seed, longPair);
                        hash = longPair.val1;
                    } catch (UnsupportedEncodingException e) {
                        e.printStackTrace();
                    }
                    count++;

                    if (hash < sketchHashes[0]) {
                        if ((parameters.bloomFilter && bloomFilter.contains(kmer))||!parameters.bloomFilter) {
                            boolean contains = false;
                            for (int j = 0; j < parameters.sketchSize; j++) {
                                if (sketchHashes[j] == hash) {
                                    contains = true;
                                    break;
                                }
                            }
                            if (!contains) {
                                boolean last = false;
                                for (int j = 1; j < parameters.sketchSize; j++) {
                                    if (hash < sketchHashes[j]) {
                                        sketchHashes[j - 1] = sketchHashes[j];
                                    } else if (hash == sketchHashes[j]) {
                                        System.out.println("duplicate");
                                        System.out.println(hash);
                                    } else {
                                        sketchHashes[j - 1] = hash;
                                        last = true;
                                        break;
                                    }
                                }
                                if (!last) {
                                    sketchHashes[parameters.sketchSize - 1] = hash;
                                }
                            }
                        } else {
                            bloomFilter.add(kmer);
                        }
                    }
                }
            }
            else if (parameters.hashFunction==0) { //murmur3 x86_32
                for (int i = 0; i <= sequence.length() - parameters.kmerSize; i++) {
                    kmer = canonical_kmer(sequence.substring(i, i + parameters.kmerSize));
                    hash = Murmur3.murmurhash3_x86_32(kmer, 0, parameters.kmerSize, parameters.seed);

                    count++;

                    if (hash < sketchHashes[0]) {
                        if (parameters.bloomFilter && bloomFilter.contains(kmer)||!parameters.bloomFilter) {
                            boolean contains = false;
                            for (int j = 0; j < parameters.sketchSize; j++) {
                                if (sketchHashes[j] == hash) {
                                    contains = true;
                                    break;
                                }
                            }
                            if (!contains) {
                                boolean last = false;
                                for (int j = 1; j < parameters.sketchSize; j++) {
                                    if (hash < sketchHashes[j]) {
                                        sketchHashes[j - 1] = sketchHashes[j];
                                    } else if (hash == sketchHashes[j]) {
                                        System.out.println("duplicate");
                                        System.out.println(hash);
                                    } else {
                                        sketchHashes[j - 1] = hash;
                                        last = true;
                                        break;
                                    }
                                }
                                if (!last) {
                                    sketchHashes[parameters.sketchSize - 1] = hash;
                                }
                            }
                        } else {
                            bloomFilter.add(kmer);
                        }
                    }
                }
            }
            else{
                throw new Exception("no such Hash-Function");
            }
            System.out.println("kmers hashed: " + count + " by thread " + this.name +"\t"+parameters.sequences.size()+" genomes left");
            latch.countDown();
            sketches.add(sketch);
        }
        return sketches;
    }
}



public class Mash {
    public static void main(String[] args){

        InputParameters parameters = new InputParameters();
        parameters.parseInput(args);
        long startTime = System.currentTimeMillis();
        if (parameters.type.equals("sketch")){
            WriteReadObject.writeObjectToFile(mash_sketch(parameters),parameters.outputFile.concat(".ske"));
        }
        else if (parameters.type.equals("dist")){
            Distance[] distances = null;
            if (parameters.sketches!=null){
                distances = mash_dist(parameters.sketches,parameters);
            }
            else{
                distances = mash_dist(mash_sketch(parameters),parameters);
            }
            if (parameters.tableOutput){
                printTable(tableOutput(distances,parameters.sequenceFiles));
            }else {
                Distance[] distances1 = new Distance[distances.length*2+parameters.sequenceFiles.size()];
                int i=0;
                for (i=0;i<distances.length;i++) {
                    distances1[2*i]=distances[i];
                    distances1[2*i+1]=new Distance(distances[i]);
                }
                i*=2;
                for (String s : parameters.sequenceFiles){
                    distances1[i] = new Distance(s,s,1,0,0,s,s,parameters.sketchSize);
                    i++;
                }
                Arrays.sort(distances1,Comparator.comparing(a -> a.getFilePath1()));
                Arrays.sort(distances1,Comparator.comparing(a -> a.getFilePath2()));
                System.out.println();
                WriteReadObject.writeTxtFile(distances, parameters.outputFile);
                for (Distance d : distances1) {

//                    System.out.print(d.getSameHashes()+", ");
                    d.printShort();
                }
                System.out.println();
            }
            if(parameters.outputPhylip){
                WriteReadObject.writePhylipFile(tableOutput(distances,parameters.sequenceFiles),parameters.outputFile);
            }
        }
        else if (parameters.type.equals("info")){
            for (Sketch sketch : parameters.sketches){
                System.out.println(sketch.toString()+"\n");
            }
        }



//        mash_dist(mash_sketch(parameters),parameters)[0].print();
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        System.out.println(duration);

    }

    private static void printTable(String[][] strings) {
        for (int i=0;i<strings.length;i++){
            for (int j=0;j<strings[0].length;j++){
                System.out.print(strings[i][j]+"\t");
            }
            System.out.println();
        }
    }

    private static String[][] tableOutput(Distance[] distances, ArrayList<String> sequenceFiles) {
        String[][] table = new String[sequenceFiles.size()+1][sequenceFiles.size()+1];
        table[0][0]="Mash-Distance";
        for (int i=1;i<sequenceFiles.size()+1;i++){
            table[i][0] = sequenceFiles.get(i-1);
            table[0][i] = sequenceFiles.get(i-1);
            table[i][i] = "0";
        }
        for (Distance distance : distances){
            int i,j;
            for (i=0;i<sequenceFiles.size();i++){
                if (distance.getFilePath1().equals(sequenceFiles.get(i))) break;
            }
            for (j=0;j<sequenceFiles.size();j++){
                if (distance.getFilePath2().equals(sequenceFiles.get(j))) break;
            }
            table[i+1][j+1] = String.valueOf(distance.getMash_distance());
            table[j+1][i+1] = String.valueOf(distance.getMash_distance());
        }
        return table;
    }

    private static Sketch[] mash_sketch(InputParameters parameters){


        return mash_sketch_in_parallel(parameters);
    }

    private static Sketch[] mash_sketch_in_parallel(InputParameters parameters){
        int threads=parameters.sequences.size();
        ExecutorService pool = Executors.newFixedThreadPool(parameters.cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<Sketch>>[] results = new Future[parameters.cores];
        for (int i = 0; i < parameters.cores; i++)
        {
            KmerTask task = new KmerTask("Task " + i, latch,parameters);
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

    private static Distance[] mash_dist (Sketch[] sketches, InputParameters parameters){
        int length = sketches.length;
        ArrayList<Sketch[]> pairOfSketches = new ArrayList<>();
        for (int i =0; i<length;i++){
            for (int j =i+1; j<length;j++){
                pairOfSketches.add(new Sketch[]{sketches[i],sketches[j]});
            }

        }
        parameters.pairsOfSketches = new SynchronizedList<>(pairOfSketches);
        return mash_dist_in_parallel(parameters);


    }

    private static Distance[] mash_dist_in_parallel(InputParameters parameters) {
        int threads=parameters.pairsOfSketches.size();
        ExecutorService pool = Executors.newFixedThreadPool(parameters.cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<Distance>>[] results = new Future[parameters.cores];
        for (int i = 0; i < parameters.cores; i++)
        {
            DistTask task = new DistTask("Task " + i, latch, parameters);
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






}
