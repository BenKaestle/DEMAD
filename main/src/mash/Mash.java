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

/*
 *  Mash.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/**
 *
 * Benjamin Kaestle, 3.2020
 */

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

    /**
     * callable method for mash distance tasks:
     * generates a list of MashDistance objects, as long as the Synchronized list containing MashSketch elements is
     * not empty. This way only a set number of Dist tasks are generated and reused
     * @return
     */
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

    /**
     * binomial function for n binom k
     * @param n
     * @param k
     * @return
     */
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

    /**
     * adds a hash value into the sketch of hashes if all conditions apply. additionally a bloom filter can be used
     * @param hash
     * @param sketchHashes
     * @param parameters
     * @param bloomFilter
     * @param kmer
     */
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
    /**
     * callable method for mash kmer tasks:
     * generates a list of MashSketch objects, as long as the Synchronized list containing input fasta filepaths is
     * not empty. This way only a set number of kmer tasks are generated and reused
     * @return
     */
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

    /**
     * for a given sliding window position i in a sequence with the reverse complement reverseSequence, the
     * lexicographically smaller one gets returned
     * @param i
     * @param sequence
     * @param reverseSequence
     * @return
     */
    private String createKmer(int i, String sequence, String reverseSequence) {
        kmer = sequence.substring(i, i + parameters.kmerSize);
        if (reverseSequence.length()==0) return kmer;
        reverseKmer = reverseSequence.substring(sequenceLength - i - parameters.kmerSize, sequenceLength - i);
        kmer = kmer.compareTo(reverseKmer) > 0 ? reverseKmer : kmer;
        return kmer;
    }
}


public class Mash {
    /**
     * main function of mash. parses inputs and performs all calculations and outputs
     * @param args
     */
    public static void mash(String[] args) {

        InputParameters parameters = new InputParameters();
        parameters.parseInput(args,true);
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
    }

    /**
     * prints for a 2d array strings the corresponding table into the command line
     * @param strings
     */
    private static void printTable(String[][] strings) {
        for (int i = 0; i < strings.length; i++) {
            for (int j = 0; j < strings[0].length; j++) {
                System.out.print(strings[i][j] + "\t");
            }
            System.out.println();
        }
    }

    /**
     * for all calculated MashDistances and the corresponding sequence file names a distance matrix gets generated and
     * returned
     * @param mashDistances
     * @param sequenceFiles
     * @return
     */
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

    /**
     * for the parsed input parameters several threads are created and organized
     * the resulting MashSketch elements are sorted into one list of Sketches and returned
     * @param parameters
     * @return
     */
    private static MashSketch[] mash_sketch(InputParameters parameters) {
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

    /**
     * for input parameters and an already created list of MashSketch objects all pairwise distances get computed. Several threads produce
     * lists of distances that get sorted and combined in the end.
     * @param mashSketches
     * @param parameters
     * @return
     */
    private static MashDistance[] mash_dist(MashSketch[] mashSketches,InputParameters parameters) {
        parameters.mashSketches = mashSketches;
        parameters.mashSketchesSynch = new SynchronizedList<MashSketch>(new ArrayList<MashSketch>(Arrays.asList(mashSketches)));
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
