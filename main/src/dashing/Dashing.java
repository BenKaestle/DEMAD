package dashing;

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
 *  Dashing.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

final class DistTask implements Callable<ArrayList<DashingDistance>> {
    private String name;
    private DashingSketch dashingSketch_1;
    private float harmonicMean_1, harmonicMean_2, harmonicMean_c;
    private int[] combinedRegister;
    private ArrayList<DashingDistance> dashingDistances;
    private float jaccard_index;
    private float mash_distance;
    private InputParameters parameters;
    private CountDownLatch latch;
    private int registerSize;
    final private float ALPHA_M = (float) (1 / (2 * Math.log(2)));

    public DistTask(String name, CountDownLatch latch, InputParameters parameters) {
        this.name = name;
        this.latch = latch;
        this.parameters = parameters;
        this.combinedRegister = new int[(int) Math.pow(2, parameters.prefixSize)];
    }

    public String getName() {
        return name;
    }

    /**
     * estimates cardinalities of sets with the register values by calculating harmonic means
     * @param register
     * @return
     */
    public float harmonicMean(int[] register) {
        registerSize = register.length;
        float denominator = 0;
        for (int i = 0; i < registerSize; i++) {
            denominator += Math.pow(2, -(register[i] + 1));
        }
        return ALPHA_M * registerSize * registerSize / denominator;
    }

    /**
     * callable method for dashing distance tasks:
     * generates a list of DashingDistance objects, as long as the Synchronized list containing DashingSketch elements is
     * not empty. This way only a set number of Dist tasks are generated and reused
     * @return
     */
    @Override
    public ArrayList<DashingDistance> call() {
        dashingDistances = new ArrayList<>();
        while (!parameters.dashingSketchesSynch.isEmpty()) {
            dashingSketch_1 = parameters.dashingSketchesSynch.get();
            if (dashingSketch_1 == null) break;
            for (DashingSketch dashingSketch_2 : parameters.dashingSketches) {
                if (dashingSketch_1.hasNoHarmonicMean()) {
                    dashingSketch_1.setHarmonicMean(harmonicMean(dashingSketch_1.getRegister()));
                    dashingSketch_1.setNoHarmonicMean(false);
                }
                if (dashingSketch_2.hasNoHarmonicMean()) {
                    dashingSketch_2.setHarmonicMean(harmonicMean(dashingSketch_2.getRegister()));
                    dashingSketch_2.setNoHarmonicMean(false);
                }
                harmonicMean_1 = dashingSketch_1.getHarmonicMean();
                harmonicMean_2 = dashingSketch_2.getHarmonicMean();
                for (int i = 0; i < dashingSketch_1.getRegister().length; i++) {
                    combinedRegister[i] = Math.max(dashingSketch_1.getRegister()[i], dashingSketch_2.getRegister()[i]);
                }
                harmonicMean_c = harmonicMean(combinedRegister);
                jaccard_index = (harmonicMean_1 + harmonicMean_2 - harmonicMean_c) / harmonicMean_c;
                if (jaccard_index < 0) jaccard_index = 0;
                mash_distance = (jaccard_index == 0f) ? 1 : -1f / parameters.kmerSize * (float) Math.log((2 * jaccard_index) / (1 + jaccard_index));
                dashingDistances.add(new DashingDistance(dashingSketch_1.getHeader(), dashingSketch_2.getHeader(), jaccard_index, dashingSketch_1.getFilename(), dashingSketch_2.getFilename(),mash_distance));
            }
            if(parameters.dashingSketchesSynch.size()>0) {
                System.out.println("comparison finished by thread " + this.name + "\t" + (parameters.dashingSketchesSynch.size() + parameters.cores) + " comparisons left");
            } else{
                System.out.println("comparison finished by thread " + this.name + "\tless than " + parameters.cores + " comparisons left");
            }
            latch.countDown();
        }
        return dashingDistances;
    }
}


final class KmerTask implements Callable<ArrayList<DashingSketch>> {
    private String name;
    private CountDownLatch latch;
    private DashingSketch dashingSketch;
    private ArrayList<DashingSketch> dashingSketches;
    private BloomFilter bloomFilter;
    private String kmer;
    private String reverseKmer;
    private InputParameters parameters;
    private int sequenceLength;
    private String sequence;
    private String reverseSequence;
    private int registerSize;
    private int registerKey;
    private HashFunction hashFunction;


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

    public static int countZeros(long val, int prefixSize) {
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
            return Math.min(n - 2, 64 - prefixSize);
        return Math.min(n - (int) val, 64 - prefixSize);
    }

    /**
     * adds a hash value into the sketch of hashes if all conditions apply. additionally a bloom filter can be used
     * @param hash
     * @param register
     */
    public void addHash(long hash, int[] register) {
        if ((parameters.bloomFilter && bloomFilter.contains(kmer)) || !parameters.bloomFilter) {
            registerKey = (int) (hash % (long) registerSize);
            if (registerKey < 0) registerKey = registerKey + registerSize;
            register[registerKey] = Math.max(register[registerKey], countZeros(hash, parameters.prefixSize));
        } else {
            bloomFilter.add(kmer);
        }
    }

    /**
     * callable method for mash kmer tasks:
     * generates a list of MashSketch objects, as long as the Synchronized list containing input fasta filepaths is
     * not empty. This way only a set number of kmer tasks are generated and reused
     * @return
     */
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
                sequenceInfo = FastaParser.parseFasta(new File(path), parameters.alphabet,parameters.reduceAlphabet);
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
            if (parameters.bloomFilter)
                bloomFilter = new BloomFilter(parameters.bloomFilterSize, parameters.bloomFilterHashes);
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
            for (int i = 0; i <= sequenceLength - parameters.kmerSize; i++) {
                kmer = createKmer(i,sequence,reverseSequence);
                hash = hashFunction.hash(kmer);
                count++;
                addHash(hash, register);
            }

            if(parameters.sequences.size()>0) {
                System.out.println("kmers hashed: " + count + " by thread " + this.name + "\t" + (parameters.cores + parameters.sequences.size()) + " genomes left");
            } else{
                System.out.println("kmers hashed: " + count + " by thread " + this.name + "\t less than " + parameters.cores + " genomes left");
            }
            latch.countDown();
            dashingSketches.add(dashingSketch);
        }
        return dashingSketches;
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

    public String getName() {
        return this.name;
    }
}

public class Dashing {
    /**
     * main function of dashing. parses inputs and performs all calculations and outputs
     * @param args
     */
    public static void dashing(String[] args) {

        InputParameters parameters = new InputParameters();
        parameters.parseInput(args,false);
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

                Arrays.sort(dashingDistances, Comparator.comparing(a -> a.getFilePath1()));
                Arrays.sort(dashingDistances, Comparator.comparing(a -> a.getFilePath2()));
                WriteReadObject.writeTxtFile(dashingDistances, parameters.outputFile);
                for (DashingDistance d : dashingDistances) {
                    System.out.println(d.toString());
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

    }
    /**
     * for input parameters and an already created list of DashingSketch objects all pairwise distances get computed.
     * Several threads produce lists of distances that get sorted and combined in the end.
     * @param dashingSketches
     * @param parameters
     * @return
     */
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

    /**
     * for the parsed input parameters several threads are created and organized
     * the resulting DashingSketch elements are sorted into one list of Sketches and returned
     * @param parameters
     * @return
     */
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

    /**
     * prints for a 2d array tableOutput the corresponding table into the command line
     * @param tableOutput
     */
    private static void printTable(String[][] tableOutput) {
        for (int i = 0; i < tableOutput.length; i++) {
            for (int j = 0; j < tableOutput[0].length; j++) {
                System.out.print(tableOutput[i][j] + "\t");
            }
            System.out.println();
        }
    }

    /**
     * for all calculated DashingDistances and the corresponding sequence file names a distance matrix gets generated and
     * returned
     * @param dashingDistances
     * @param sequenceFiles
     * @return
     */
    private static String[][] tableOutput(DashingDistance[] dashingDistances, ArrayList<String> sequenceFiles) {
        String[][] table = new String[sequenceFiles.size() + 1][sequenceFiles.size() + 1];
        table[0][0] = "Mash-Distance";
        for (int i = 1; i < sequenceFiles.size() + 1; i++) {
            table[i][0] = sequenceFiles.get(i - 1);
            table[0][i] = sequenceFiles.get(i - 1);
            table[i][i] = "0";
        }
        for (DashingDistance dashingDistance : dashingDistances) {
            int i, j;
            for (i = 0; i < sequenceFiles.size(); i++) {
                if (dashingDistance.getFilePath1().equals(sequenceFiles.get(i))) break;
            }
            for (j = 0; j < sequenceFiles.size(); j++) {
                if (dashingDistance.getFilePath2().equals(sequenceFiles.get(j))) break;
            }
            table[i + 1][j + 1] = String.format(Locale.ROOT,"%f", dashingDistance.getMash_distance());
            table[j + 1][i + 1] = String.format(Locale.ROOT,"%f", dashingDistance.getMash_distance());
        }
        return table;
    }


}