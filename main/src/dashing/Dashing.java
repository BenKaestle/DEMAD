package dashing;

import hash_functions.Murmur3;
import hash_functions.Wang;
import mash.MashSketch;
import utility.*;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Comparator;
import java.util.concurrent.*;

final class DistTask implements Callable<ArrayList<DashingDistance>>
{
    private String name;
    private DashingSketch dashingSketch_1;
    private float harmonicMean_1, harmonicMean_2, harmonicMean_c;
    private int[] combinedRegister;
    private int sketch_length;
    private int same_hash_counter;
    private ArrayList<DashingDistance> dashingDistances;
    private float jaccard_index, mash_distance;
    private double p_value;
    private InputParameters parameters;
    private CountDownLatch latch;
    private float r, pkx, pky;
    private final int ALPHABET_SIZE = 4;
    private int registerSize;
    final private float ALPHA_M = (float)(1/(2*Math.log(2)));

    public DistTask(String name, CountDownLatch latch, InputParameters parameters) {
        this.name = name;
        this.latch = latch;
        this.parameters = parameters;
        this.combinedRegister = new int[(int)Math.pow(2,parameters.prefixSize)];
    }

    public String getName() {
        return name;
    }

    public float harmonicMean(int[] register){
        registerSize=register.length;
        float denominator=0;
        for(int i=0;i<registerSize;i++){
            denominator+=Math.pow(2,-(register[i]+1));
        }
        return ALPHA_M*registerSize*registerSize/denominator;
    }

    @Override
    public ArrayList<DashingDistance> call(){
        dashingDistances = new ArrayList<>();
        while (!parameters.dashingSketchesSynch.isEmpty()) {
            dashingSketch_1 = parameters.dashingSketchesSynch.get();
            if (dashingSketch_1 ==null) break;
            for(DashingSketch dashingSketch_2:parameters.dashingSketches){
                if (dashingSketch_1.hasNoHarmonicMean()){
                    dashingSketch_1.setHarmonicMean(harmonicMean(dashingSketch_1.getRegister()));
                    dashingSketch_1.setNoHarmonicMean(false);
                }
                if (dashingSketch_2.hasNoHarmonicMean()){
                    dashingSketch_2.setHarmonicMean(harmonicMean(dashingSketch_2.getRegister()));
                    dashingSketch_2.setNoHarmonicMean(false);
                }
                harmonicMean_1=dashingSketch_1.getHarmonicMean();
                harmonicMean_2=dashingSketch_2.getHarmonicMean();
                for (int i =0;i<dashingSketch_1.getRegister().length;i++){
                    combinedRegister[i]=Math.max(dashingSketch_1.getRegister()[i],dashingSketch_2.getRegister()[i]);
                }
                harmonicMean_c = harmonicMean(combinedRegister);
                jaccard_index=(harmonicMean_1+harmonicMean_2-harmonicMean_c)/harmonicMean_c;
                if (jaccard_index<0)jaccard_index=0;
//                System.out.println(jaccard_index + " "+harmonicMean_1 + " "+harmonicMean_2 + " "+harmonicMean_c + " ");





                dashingDistances.add(new DashingDistance(dashingSketch_1.getHeader(), dashingSketch_2.getHeader(), jaccard_index, 0, 0, dashingSketch_1.getFilename(), dashingSketch_2.getFilename(), 0));
            }
            System.out.println("comparison finished by thread " + this.name + "\t" + (parameters.dashingSketchesSynch.size()+parameters.cores) + " comparisons left");
            latch.countDown();
        }
        return dashingDistances;
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



final class KmerTask implements Callable<ArrayList<DashingSketch>> {
    private String name;
    private CountDownLatch latch;
    private DashingSketch dashingSketch;
    private ArrayList<DashingSketch> dashingSketches;
    private BloomFilter bloomFilter;
    private SketchSet sketchSet;
    private String kmer;
    private String reverseKmer;
    private InputParameters parameters;
    private int sequenceLength;
    private String sequence;
    private String reverseSequence;
    private int registerSize;
    private int registerKey;


    public KmerTask(String name, CountDownLatch latch, InputParameters parameters) {
        this.name = name;
        this.latch = latch;
        this.parameters = parameters;
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

    public static int countZeros(long val, int prefixSize)
    {
        long y;
        int n = 64;
        y = val >> 32;
        if (y != 0) {
            n = n - 32;
            val = y;
        }
        y = val >> 16;
        if (y != 0) {
            n = n - 16;
            val = y;
        }
        y = val >> 8;
        if (y != 0) {
            n = n - 8;
            val = y;
        }
        y = val >> 4;
        if (y != 0) {
            n = n - 4;
            val = y;
        }
        y = val >> 2;
        if (y != 0) {
            n = n - 2;
            val = y;
        }
        y = val >> 1;
        if (y != 0)
            return Math.min(n - 2,64-prefixSize);
        return Math.min(n - (int)val, 64-prefixSize);
    }

    @Override
    public ArrayList<DashingSketch> call() throws Exception {
        dashingSketches = new ArrayList<>();
        String[] sequenceInfo = new String[2];
        Murmur3.LongPair longPair = new Murmur3.LongPair();
        while (!parameters.sequences.isEmpty()) {
            registerSize = (int) Math.pow(2, parameters.prefixSize);
            int[] register = new int[registerSize];
            String path = parameters.sequences.get();
            if (path == null) break;
            try {
                sequenceInfo = FastaParser.parseFasta(new File(path));
            } catch (IOException e) {
                e.printStackTrace();
            }
            int count = 0;
            sequence = sequenceInfo[1];
            reverseSequence = sequenceInfo[2];
            sequenceLength = sequence.length();
            if (optimalK(sequenceLength, parameters.lowKThreshold) > parameters.kmerSize) {
                System.out.println("WARNING: For the k-mer size used (" + parameters.kmerSize + "), the random match probability is above the specified warning threshold (" + parameters.lowKThreshold +
                        ") for the sequence \"" + path + "\" of size " + sequenceLength + ". Distances to this sequence may be underestimated as a result. To meet the threshold of " + parameters.lowKThreshold +
                        ", a k-mer size of at least " + optimalK(sequenceLength, parameters.lowKThreshold) + " is required. See: -k, -w.");
            }
            dashingSketch = new DashingSketch(register, sequenceInfo[0], path, parameters.kmerSize, parameters.hashFunction, sequenceLength, parameters.seed);
            long hash = 0;
            if (parameters.bloomFilter) bloomFilter = new BloomFilter(parameters.bloomFilterSize, parameters.bloomFilterHashes);
            sketchSet = new SketchSet(5);
            if (parameters.hashFunction == 1) { //murmur3 x64_128
                for (int i = 0; i <= sequenceLength - parameters.kmerSize; i++) {
                    kmer = sequence.substring(i, i + parameters.kmerSize);
                    reverseKmer = reverseSequence.substring(sequenceLength - i - parameters.kmerSize, sequenceLength - i);
                    kmer = kmer.compareTo(reverseKmer) > 0 ? reverseKmer : kmer;
                    try {
                        Murmur3.murmurhash3_x64_128(kmer.getBytes("UTF-8"), 0, parameters.kmerSize, parameters.seed, longPair);
                        hash = longPair.val1;
                    } catch (UnsupportedEncodingException e) {
                        e.printStackTrace();
                    }
                    count++;
                    if ((parameters.bloomFilter && bloomFilter.contains(kmer)) || !parameters.bloomFilter) {
                        registerKey=(int)(hash%(long)registerSize);
                        if (registerKey<0)registerKey=registerKey+registerSize;
                        register[registerKey]=Math.max(register[registerKey],countZeros(hash,parameters.prefixSize));
                    }
                    else{
                        bloomFilter.add(kmer);
                    }
                }
            }
            else if (parameters.hashFunction == 2) { //wang
                for (int i = 0; i <= sequenceLength - parameters.kmerSize; i++) {
                    kmer = sequence.substring(i, i + parameters.kmerSize);
                    reverseKmer = reverseSequence.substring(sequenceLength - i - parameters.kmerSize, sequenceLength - i);
                    kmer = kmer.compareTo(reverseKmer) > 0 ? reverseKmer : kmer;
                    hash = Wang.hash(kmer);
                    count++;
                    if ((parameters.bloomFilter && bloomFilter.contains(kmer)) || !parameters.bloomFilter) {
                        registerKey=(int)(hash%(long)registerSize);
                        if (registerKey<0)registerKey=registerKey+registerSize;
                        register[registerKey]=Math.max(register[registerKey],countZeros(hash,parameters.prefixSize));
                    }
                    else{
                        bloomFilter.add(kmer);
                    }
                }
            }
            System.out.println("kmers hashed: " + count + " by thread " + this.name +"\t"+(parameters.cores+parameters.sequences.size())+" genomes left");
            latch.countDown();
            dashingSketches.add(dashingSketch);
        }
        return dashingSketches;
    }

    public String getName() {
        return this.name;
    }
}

public class Dashing {
    public static void main(String[] args) {
        dashing(args);
    }

    public static void dashing(String[] args) {

        InputParameters parameters = new InputParameters();
        parameters.parseInputDashing(args);
        long startTime = System.currentTimeMillis();
        if (parameters.type.equals("sketch")) {
            WriteReadObject.writeObjectToFile(dashingSketch(parameters), parameters.outputFile.concat(".dsk"));
        } else if (parameters.type.equals("dist")) {
            DashingDistance[] dashingDistances = null;
            if (parameters.mashSketches != null) {
                dashingDistances = dashingDist(parameters.dashingSketches, parameters);
            } else {
                dashingDistances = dashingDist(dashingSketch(parameters), parameters);
            }
            if (parameters.tableOutput) {
                printTable(tableOutput(dashingDistances, parameters.sequenceFiles));
            } else {

//                Arrays.sort(dashingDistances, Comparator.comparing(a -> a.getFilePath1()));
//                Arrays.sort(dashingDistances, Comparator.comparing(a -> a.getFilePath2()));
//                System.out.println();
//                WriteReadObject.writeTxtFile(dashingDistances, parameters.outputFile);
                for (DashingDistance d : dashingDistances) {


//                    d.printShort();
                    System.out.println(d.toStringJaccard());
                }
                System.out.println();
            }
            if (parameters.outputPhylip) {
                WriteReadObject.writePhylipFile(tableOutput(dashingDistances, parameters.sequenceFiles), parameters.outputFile);
            }
        } else if (parameters.type.equals("info")) {
            for (DashingSketch dashingSketch : parameters.dashingSketches) {
                System.out.println(dashingSketch.toString() + "\n");
            }
        }


//        mash_dist(mash_sketch(parameters),parameters)[0].print();
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        System.out.println(duration);

    }

    private static DashingDistance[] dashingDist(DashingSketch[] dashingSketches, InputParameters parameters) {
        parameters.dashingSketches = dashingSketches;
        parameters.dashingSketchesSynch = new SynchronizedList<DashingSketch>(new ArrayList<DashingSketch>(Arrays.asList(dashingSketches)));
        int threads = dashingSketches.length;
        ExecutorService pool = Executors.newFixedThreadPool(parameters.cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<DashingDistance>>[] results = new Future[parameters.cores];
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
        DashingDistance[] dashingDistances = new DashingDistance[threads * threads];
        int pointer = 0;
        for (Future<ArrayList<DashingDistance>> distance : results) {
            try {
                for (DashingDistance d : distance.get()) {
                    dashingDistances[pointer] = d;
                    pointer++;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return dashingDistances;
    }

    private static DashingSketch[] dashingSketch(InputParameters parameters) {
        int threads = parameters.sequences.size();
        ExecutorService pool = Executors.newFixedThreadPool(parameters.cores);
        CountDownLatch latch = new CountDownLatch(threads);
        Future<ArrayList<DashingSketch>>[] results = new Future[parameters.cores];
        for (int i = 0; i < parameters.cores; i++) {
            dashing.KmerTask task = new dashing.KmerTask("Task " + i, latch, parameters);
            System.out.println("A new task has been added : " + task.getName());
            results[i] = pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        pool.shutdown();
        DashingSketch[] dashingSketches = new DashingSketch[threads];
        int pointer = 0;
        for (Future<ArrayList<DashingSketch>> sketch : results) {
            try {
                for (DashingSketch s : sketch.get()) {
                    dashingSketches[pointer] = s;
                    pointer++;
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
        return dashingSketches;
    }

    private static void printTable(String[][] tableOutput) {
        for (int i = 0; i < tableOutput.length; i++) {
            for (int j = 0; j < tableOutput[0].length; j++) {
                System.out.print(tableOutput[i][j] + "\t");
            }
            System.out.println();
        }
    }

    private static String[][] tableOutput(DashingDistance[] dashingDistances, ArrayList<String> sequenceFiles) {
        String[][] table = new String[sequenceFiles.size() + 1][sequenceFiles.size() + 1];
        table[0][0] = "Jaccard-Index";
        for (int i = 1; i < sequenceFiles.size() + 1; i++) {
            table[i][0] = sequenceFiles.get(i - 1);
            table[0][i] = sequenceFiles.get(i - 1);
            table[i][i] = "1";
        }
        for (DashingDistance dashingDistance : dashingDistances) {
            int i, j;
            for (i = 0; i < sequenceFiles.size(); i++) {
                if (dashingDistance.getFilePath1().equals(sequenceFiles.get(i))) break;
            }
            for (j = 0; j < sequenceFiles.size(); j++) {
                if (dashingDistance.getFilePath2().equals(sequenceFiles.get(j))) break;
            }
            table[i + 1][j + 1] = String.valueOf(dashingDistance.getJaccard_index());
            table[j + 1][i + 1] = String.valueOf(dashingDistance.getJaccard_index());
        }
        return table;
    }


}