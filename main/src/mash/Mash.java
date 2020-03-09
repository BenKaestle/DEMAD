package mash;

import hash_functions.*;
import utility.*;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import java.util.concurrent.*;


final class DistTask implements Callable<ArrayList<MashDistance>> {

    private String name;
    private MashSketch mashSketch_1;
    private long current_hash_1;
    private long current_hash_2;
    private int sketch_length;
    private int same_hash_counter;
    private ArrayList<MashDistance> mashDistances;
    private float jaccard_index, mash_distance;
    private double p_value;
    private InputParameters parameters;
    private CountDownLatch latch;
    private float r, pkx, pky;
    private int alphabetSize = 4;
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
    public ArrayList<MashDistance> call() {
        mashDistances = new ArrayList<>();
        int pointer_1, pointer_2;
        if (parameters.alphabet=="AA"){
            switch (parameters.reduceAlphabet){
                case -1:
                    alphabetSize =20;
                    break;
                case 0:
                    alphabetSize =15;
                    break;
                case 1:
                    alphabetSize =10;
                    break;
                case 2:
                    alphabetSize =8;
                    break;
                case 3:
                    alphabetSize =4;
                    break;
                case 4:
                    alphabetSize =2;
                    break;
                default:
                    alphabetSize = 20;
                    break;
            }
        }
        while (!parameters.mashSketchesSynch.isEmpty()) {
            mashSketch_1 = parameters.mashSketchesSynch.get();
            if (mashSketch_1 == null) break;
            for (MashSketch mashSketch_2 : parameters.mashSketches) {
                pointer_1 = mashSketch_1.getHashes().length - 1;
                pointer_2 = mashSketch_2.getHashes().length - 1;
                sketch_length = Math.min(mashSketch_1.getHashes().length, mashSketch_2.getHashes().length);
                same_hash_counter = 0;
                for (int i = 0; i < sketch_length; i++) {
                    current_hash_1 = mashSketch_1.getHashes()[pointer_1];
                    current_hash_2 = mashSketch_2.getHashes()[pointer_2];
                    if (current_hash_1 > current_hash_2) {
                        pointer_2--;
                    }
                    if (current_hash_1 < current_hash_2) {
                        pointer_1--;
                    }
                    if (current_hash_1 == current_hash_2) {
                        same_hash_counter++;
                        pointer_1--;
                        pointer_2--;
                    }
                }
                jaccard_index = (float) same_hash_counter / sketch_length;
                mash_distance = (jaccard_index == 0f) ? 1 : -1f / parameters.kmerSize * (float) Math.log((2 * jaccard_index) / (1 + jaccard_index));
                genome_size_1 = mashSketch_1.getGenome_size();
                genome_size_2 = mashSketch_2.getGenome_size();
                pkx = 1 - (float) Math.pow(1 - Math.pow(alphabetSize, -parameters.kmerSize), genome_size_1);
                pky = 1 - (float) Math.pow(1 - Math.pow(alphabetSize, -parameters.kmerSize), genome_size_2);
                r = pkx * pky / (pkx + pky - pkx * pky);
                p_value = 1;
                for (int i = 0; i < same_hash_counter; i++) {
                    p_value -= (binom(sketch_length, i) * Math.pow(r, i) * Math.pow((1 - r), (sketch_length - i)));
                }

                if (mash_distance<=0) mash_distance=0;
                if (p_value<0) p_value=0;
                mashDistances.add(new MashDistance(mashSketch_1.getHeader(), mashSketch_2.getHeader(), jaccard_index, p_value, mash_distance, mashSketch_1.getFilename(), mashSketch_2.getFilename(), same_hash_counter));
            }
            if(parameters.mashSketchesSynch.size()>0) {
                System.out.println("comparison finished by thread " + this.name + "\t" + (parameters.mashSketchesSynch.size() + parameters.cores) + " comparisons left");
            } else{
                System.out.println("comparison finished by thread " + this.name + "\tless than " + parameters.cores + " comparisons left");
            }
            latch.countDown();
        }
        return mashDistances;
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

final class KmerTask implements Callable<ArrayList<MashSketch>> {
    private String name;
    private CountDownLatch latch;
    private MashSketch mashSketch;
    private ArrayList<MashSketch> mashSketches;
    private BloomFilter bloomFilter;
    private String kmer;
    private String reverseKmer;
    private InputParameters parameters;
    private int sequenceLength;
    private HashFunction hashFunction;

    public KmerTask(String name, CountDownLatch latch, InputParameters parameters) {
        this.name = name;
        this.latch = latch;
        this.parameters = parameters;
    }

    public String getName() {
        return name;
    }
    /**
     * for a given genome size and the desired probability prob of observing a random k-mer calculate the optimal k value
     *
     * @param genome_size
     * @param prob        - probability of observing a random k-mer
     * @return
     */
    public static int optimalK(long genome_size, float prob) {
        return (int) Math.ceil(Math.log(genome_size * (1 - prob) / prob) / Math.log(4));
    }

    public static void addHash(long hash, long[] sketchHashes, InputParameters parameters, BloomFilter bloomFilter, String kmer) {
        if (hash < sketchHashes[0]) {
            if ((parameters.bloomFilter && bloomFilter.contains(kmer)) || !parameters.bloomFilter) {
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

    @Override
    public ArrayList<MashSketch> call() throws Exception {
        mashSketches = new ArrayList<>();
        String[] sequenceInfo = new String[2];
        Murmur3.LongPair longPair = new Murmur3.LongPair();
        while (!parameters.sequences.isEmpty()) {
            long[] sketchHashes = new long[parameters.sketchSize];
            String path = parameters.sequences.get();
            if (path == null) break;
            try {
                sequenceInfo = FastaParser.parseFasta(new File(path),parameters.alphabet,parameters.reduceAlphabet);
            } catch (IOException e) {
                e.printStackTrace();
            }
            int count = 0;
            String sequence = sequenceInfo[1];
            String reverseSequence = sequenceInfo[2];
            if (optimalK(sequence.length(), parameters.lowKThreshold) > parameters.kmerSize) {
                System.out.println("WARNING: For the k-mer size used (" + parameters.kmerSize + "), the random match probability is above the specified warning threshold (" + parameters.lowKThreshold +
                        ") for the sequence \"" + path + "\" of size " + sequence.length() + ". Distances to this sequence may be underestimated as a result. To meet the threshold of " + parameters.lowKThreshold +
                        ", a k-mer size of at least " + optimalK(sequence.length(), parameters.lowKThreshold) + " is required. See: -k, -w.");
            }
            Arrays.fill(sketchHashes, Long.MAX_VALUE);
            mashSketch = new MashSketch(sketchHashes, sequenceInfo[0], path, parameters.kmerSize, parameters.hashFunction, sequence.length(), parameters.seed);
            long hash = 0;
            if (parameters.bloomFilter)
                bloomFilter = new BloomFilter(parameters.bloomFilterSize, parameters.bloomFilterHashes);
            sequenceLength = sequence.length();

            if (parameters.hashFunction == 1) { //murmur3 x64_128
                hashFunction = new Murmur3_64(parameters.kmerSize, parameters.seed);
            } else if (parameters.hashFunction == 0) { //murmur3 x86_32
                hashFunction = new Murmur3_32(parameters.kmerSize, parameters.seed);
            } else if (parameters.hashFunction == 2) { //wang hash
                hashFunction = new Wang();
            } else if (parameters.hashFunction == 3) { //java hash_code
                hashFunction = new JavaHashCode();
            } else {
                throw new Exception("no such Hash-Function");
            }
            for (int i = 0; i <= sequence.length() - parameters.kmerSize; i++) {
                kmer = createKmer(i,sequence,reverseSequence);
                hash = hashFunction.hash(kmer);
                count++;
                addHash(hash, sketchHashes, parameters, bloomFilter, kmer);
            }




            if(parameters.sequences.size()>0) {
                System.out.println("kmers hashed: " + count + " by thread " + this.name + "\t" + (parameters.cores + parameters.sequences.size()) + " genomes left");
            } else{
                System.out.println("kmers hashed: " + count + " by thread " + this.name + "\t less than " + parameters.cores + " genomes left");
            }
            latch.countDown();
            mashSketches.add(mashSketch);
        }
        return mashSketches;
    }

    private String createKmer(int i, String sequence, String reverseSequence) {
        kmer = sequence.substring(i, i + parameters.kmerSize);
        if (reverseSequence.length()==0) return kmer;
        reverseKmer = reverseSequence.substring(sequenceLength - i - parameters.kmerSize, sequenceLength - i);
        kmer = kmer.compareTo(reverseKmer) > 0 ? reverseKmer : kmer;
        return kmer;
    }
}


public class Mash {
    public static void main(String[] args) {
        mash(args);
    }

    public static void mash(String[] args) {

        InputParameters parameters = new InputParameters();
        parameters.parseInput(args,true);
        long startTime = System.currentTimeMillis();
        if (parameters.type.equals("sketch")) {
            WriteReadObject.writeObjectToFile(mash_sketch(parameters), parameters.outputFile.concat(".msk"));
        } else if (parameters.type.equals("dist")) {
            MashDistance[] mashDistances;
            if (parameters.mashSketches != null) {
                mashDistances = mash_dist(parameters.mashSketches, parameters);
            } else {
                mashDistances = mash_dist(mash_sketch(parameters), parameters);
            }
            if (parameters.tableOutput) {
                printTable(tableOutput(mashDistances, parameters.sequenceFiles));
            } else {
                Arrays.sort(mashDistances, Comparator.comparing(a -> a.getFilePath1()));
                Arrays.sort(mashDistances, Comparator.comparing(a -> a.getFilePath2()));
                WriteReadObject.writeTxtFile(mashDistances, parameters.outputFile);
                for (MashDistance d : mashDistances) {
                    System.out.println(d.toString());
                }
                System.out.println();
            }
            if (parameters.outputPhylip) {
                WriteReadObject.writePhylipFile(tableOutput(mashDistances, parameters.sequenceFiles), parameters.outputFile);
            }
        } else if (parameters.type.equals("info")) {
            for (MashSketch mashSketch : parameters.mashSketches) {
                System.out.println(mashSketch.toString() + "\n");
            }
        }
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        System.out.println(duration);

    }

    private static void printTable(String[][] strings) {
        for (int i = 0; i < strings.length; i++) {
            for (int j = 0; j < strings[0].length; j++) {
                System.out.print(strings[i][j] + "\t");
            }
            System.out.println();
        }
    }

    private static String[][] tableOutput(MashDistance[] mashDistances, ArrayList<String> sequenceFiles) {
        String[][] table = new String[sequenceFiles.size() + 1][sequenceFiles.size() + 1];
        table[0][0] = "Mash-Distance";
        for (int i = 1; i < sequenceFiles.size() + 1; i++) {
            table[i][0] = sequenceFiles.get(i - 1);
            table[0][i] = sequenceFiles.get(i - 1);
            table[i][i] = "0";
        }
        for (MashDistance mashDistance : mashDistances) {
            int i, j;
            for (i = 0; i < sequenceFiles.size(); i++) {
                if (mashDistance.getFilePath1().equals(sequenceFiles.get(i))) break;
            }
            for (j = 0; j < sequenceFiles.size(); j++) {
                if (mashDistance.getFilePath2().equals(sequenceFiles.get(j))) break;
            }
            table[i + 1][j + 1] = String.format(Locale.ROOT,"%f", mashDistance.getMash_distance());
            table[j + 1][i + 1] = String.format(Locale.ROOT,"%f", mashDistance.getMash_distance());
        }
        return table;
    }

    private static MashSketch[] mash_sketch(InputParameters parameters) {


        return mash_sketch_in_parallel(parameters);
    }

    private static MashSketch[] mash_sketch_in_parallel(InputParameters parameters) {
        int threads = parameters.sequences.size();
        ExecutorService pool = Executors.newFixedThreadPool(parameters.cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<MashSketch>>[] results = new Future[parameters.cores];
        for (int i = 0; i < parameters.cores; i++) {
            KmerTask task = new KmerTask("Task " + i, latch, parameters);
            System.out.println("A new task has been added : " + task.getName());
            results[i] = pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        pool.shutdown();
        MashSketch[] mashSketches = new MashSketch[threads];
        int pointer = 0;
        for (Future<ArrayList<MashSketch>> sketch : results) {
            try {
                for (MashSketch s : sketch.get()) {
                    mashSketches[pointer] = s;
                    pointer++;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return mashSketches;
    }

    private static MashDistance[] mash_dist(MashSketch[] mashSketches, InputParameters parameters) {
        parameters.mashSketches = mashSketches;
        parameters.mashSketchesSynch = new SynchronizedList<MashSketch>(new ArrayList<MashSketch>(Arrays.asList(mashSketches)));
        return mash_dist_in_parallel(parameters);
    }

    private static MashDistance[] mash_dist_in_parallel(InputParameters parameters) {
        int threads = parameters.mashSketchesSynch.size();
        ExecutorService pool = Executors.newFixedThreadPool(parameters.cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<MashDistance>>[] results = new Future[parameters.cores];
        for (int i = 0; i < parameters.cores; i++) {
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
        MashDistance[] mashDistances = new MashDistance[(int) Math.pow(threads, 2)];
        int pointer = 0;
        for (Future<ArrayList<MashDistance>> distance : results) {
            try {
                for (MashDistance d : distance.get()) {
                    mashDistances[pointer] = d;
                    pointer++;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return mashDistances;
    }


}
