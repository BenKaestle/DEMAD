package utility;

import dashing.DashingSketch;
import mash.MashSketch;
import org.apache.commons.cli.*;

import java.io.File;
import java.util.ArrayList;

public class InputParameters {
    public SynchronizedList<String> sequences;
    public SynchronizedList<MashSketch[]> pairsOfSketches;
    public int sketchSize = 1000;
    public int kmerSize = 21;
    public int cores = 4;
    public boolean bloomFilter = false;
    public int hashFunction = 1;
    public int seed=42;
    public int bloomFilterSize = 23; //bis 31
    public int bloomFilterHashes = 7; //bis 8
    public String type; //sketch, dist, info
    public MashSketch[] mashSketches;
    public DashingSketch[] dashingSketches;
    public String outputFile;
    public boolean tableOutput;
    public ArrayList<String> sequenceFiles;
    public float lowKThreshold;
    public boolean outputPhylip;
    public int prefixSize;


    public void parseInputDashing(String[] args){
        this.type = args[0].trim();
        if (!type.equals("sketch") && !type.equals("dist") && !type.equals("info")){
            System.out.println("Choose one of: sketch, dist, info");
            System.exit(1);
        }
        Options options = new Options();

        Option input=null;
        if (this.type.equals("dist"))
            input = new Option("i", "input", true, "input file path (either *.dsk file or all fasta files in directory)");
        else if (this.type.equals("sketch"))
            input = new Option("i", "input", true, "input file path (all fasta files in directory will get sketched)");
        else if (this.type.equals("info"))
            input = new Option("i", "input", true, "*.dsk file of interest");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "output file");
        if (this.type.equals("sketch")) output.setRequired(true);
        options.addOption(output);

        Option prefixSize = new Option("ps", "prefixSize", true, "length of the prefix (log2 of register size) (default = 8)");
        options.addOption(prefixSize);

        Option kmerSize = new Option("k", "kmerSize", true, "length of kmers (default = 21)");
        options.addOption(kmerSize);

        Option cores = new Option("p", "parellism", true, "this many threads will be spawned for processing (default = 4)");
        options.addOption(cores);

        Option bloomFilter = new Option("b", "bloomFilter", false, "a bloom filter is used for sketching (default = false)");
        options.addOption(bloomFilter);

        Option tableOut = new Option("t", "table", false, "output a distance matrix (only mash dist)");
        options.addOption(tableOut);

        Option hashFunction = new Option("h", "hashFunction", true, "which hash-function should be used");
        options.addOption(hashFunction);

        Option seed = new Option("S", "seed", true, "seed used for hash function (default = 42)");
        options.addOption(seed);

        Option lowKThreshold = new Option("w", "lowKThreshold", true, " Probability threshold for warning about low k-mer size. (0-1) (default = 0.01)");
        options.addOption(lowKThreshold);

        Option outputPhylip = new Option("P", "outputPhylip", false, " output the distance Matrix as a *.phylip file");
        options.addOption(outputPhylip);

        Option bloomFilterParameters = new Option("bp", "bloomFilterParameters", true, "size (log2(#bits)) and number of hash functions of the Used BloomFilter");
        bloomFilterParameters.setArgs(2);
        options.addOption(bloomFilterParameters);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("utility-name", options);

            System.exit(1);
        }

        this.update(cmd);
        this.updateDashing(cmd);
    }


    public void parseInputMash(String[] args){
        this.type = args[0].trim();
        if (!type.equals("sketch") && !type.equals("dist") && !type.equals("info")){
            System.out.println("Choose one of: sketch, dist, info");
            System.exit(1);
        }
        Options options = new Options();

        Option input=null;
        if (this.type.equals("dist"))
            input = new Option("i", "input", true, "input file path (either *.msk file or all fasta files in directory)");
        else if (this.type.equals("sketch"))
            input = new Option("i", "input", true, "input file path (all fasta files in directory will get sketched)");
        else if (this.type.equals("info"))
            input = new Option("i", "input", true, "*.msk file of interest");
        input.setRequired(true);
        options.addOption(input);

        Option output = new Option("o", "output", true, "output file");
        if (this.type.equals("sketch")) output.setRequired(true);
        options.addOption(output);

        Option sketchSize = new Option("s", "sketchSize", true, "size of the generated sketches (default = 1000)");
        options.addOption(sketchSize);

        Option kmerSize = new Option("k", "kmerSize", true, "length of kmers (default = 21)");
        options.addOption(kmerSize);

        Option cores = new Option("p", "parellism", true, "this many threads will be spawned for processing (default = 4)");
        options.addOption(cores);

        Option bloomFilter = new Option("b", "bloomFilter", false, "a bloom filter is used for sketching (default = false)");
        options.addOption(bloomFilter);

        Option tableOut = new Option("t", "table", false, "output a distance matrix (only mash dist)");
        options.addOption(tableOut);

        Option hashFunction = new Option("h", "hashFunction", true, "which hash-function should be used");
        options.addOption(hashFunction);

        Option seed = new Option("S", "seed", true, "seed used for hash function (default = 42)");
        options.addOption(seed);

        Option lowKThreshold = new Option("w", "lowKThreshold", true, " Probability threshold for warning about low k-mer size. (0-1) (default = 0.01)");
        options.addOption(lowKThreshold);

        Option outputPhylip = new Option("P", "outputPhylip", false, " output the distance Matrix as a *.phylip file");
        options.addOption(outputPhylip);

        Option bloomFilterParameters = new Option("bp", "bloomFilterParameters", true, "size (log2(#bits)) and number of hash functions of the Used BloomFilter");
        bloomFilterParameters.setArgs(2);
        options.addOption(bloomFilterParameters);


        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("utility-name", options);

            System.exit(1);
        }

        this.update(cmd);
        this.updateMash(cmd);
    }

    private void updateMash(CommandLine cmd){
        sequenceFiles = new ArrayList<>();
        String input = cmd.getOptionValue("input");
        if (input.endsWith(".msk")&&(this.type.equals("dist")||this.type.equals("info"))){
            this.mashSketches = (MashSketch[])WriteReadObject.readObjectFromFile(input);
            for (MashSketch sketch : mashSketches){
                sequenceFiles.add(sketch.getFilename());
            }
        }
        else{
            File folder = new File(input);
            File[] listOfFiles = folder.listFiles();
            if (listOfFiles==null){
                System.out.println("Input file/path: \""+ input + "\" not found");
                System.exit(1);
            }
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()&& (listOfFiles[i].getPath().endsWith(".fasta")||listOfFiles[i].getPath().endsWith(".fna"))) {
                    sequenceFiles.add(listOfFiles[i].getPath());
                }
            }
            if (sequenceFiles.size()<2){
                System.out.println("There was only one or none fasta files at given directory. Please retry with valid input.");
                System.exit(1);
            }
        }
        if (this.mashSketches ==null&&this.type.equals("info")){
            System.out.println("To get information on a sketch, you need to input a valid *.msk file");
            System.exit(1);
        }
    }

    private void updateDashing(CommandLine cmd){
        sequenceFiles = new ArrayList<>();
        String input = cmd.getOptionValue("input");
        if (input.endsWith(".dsk")&&(this.type.equals("dist")||this.type.equals("info"))){
            this.dashingSketches = (DashingSketch[])WriteReadObject.readObjectFromFile(input);
            for (DashingSketch sketch : dashingSketches){
                sequenceFiles.add(sketch.getFilename());
            }
        }
        else{
            File folder = new File(input);
            File[] listOfFiles = folder.listFiles();
            if (listOfFiles==null){
                System.out.println("Input file/path: \""+ input + "\" not found");
                System.exit(1);
            }
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()&& (listOfFiles[i].getPath().endsWith(".fasta")||listOfFiles[i].getPath().endsWith(".fna"))) {
                    sequenceFiles.add(listOfFiles[i].getPath());
                }
            }
            if (sequenceFiles.size()<2){
                System.out.println("There was only one or none fasta files at given directory. Please retry with valid input.");
                System.exit(1);
            }
        }
        if (this.dashingSketches ==null&&this.type.equals("info")){
            System.out.println("To get information on a sketch, you need to input a valid *.dsk file");
            System.exit(1);
        }
    }

    private void update(CommandLine cmd) {
        this.tableOutput = cmd.hasOption("table");
        this.bloomFilter = cmd.hasOption("bloomFilter");
        this.outputPhylip = cmd.hasOption("outputPhylip");
        if (this.bloomFilter){
            if (cmd.hasOption("bloomFilterParameters")){
                this.bloomFilterSize = Integer.parseInt(cmd.getOptionValues("bloomFilterParameters")[0]);
                this.bloomFilterHashes = Integer.parseInt(cmd.getOptionValues("bloomFilterParameters")[1]);

            } else{
                calculateOptimalBloomFilter(0.01f,1000000);
            }
            if (this.bloomFilterHashes>8){
                System.out.println("WARNING: maximum number of bloom filter hashes is 8 -> your value got reduced to 8");
                this.bloomFilterHashes=8;
            }
            if (this.bloomFilterHashes<1){
                System.out.println("WARNING: minimum number of bloom filter hashes is 1 -> your value got increased to 1");
                this.bloomFilterHashes=1;
            }
            if (this.bloomFilterSize<1){
                System.out.println("WARNING: minimum size of bloom filter is 1 -> your value got increased to 1");
                this.bloomFilterSize=1;
            }
            if (this.bloomFilterSize>31){
                System.out.println("WARNING: maximum size of bloom filter is 31 -> your value got reduced to 31");
                this.bloomFilterSize=31;
            }
            if (calculatePvalueBloomFilter(1000000)>0.01f){
                System.out.println("WARNING: possibility of false positive in bloom filter higher than 1%");
            }

        }
        this.sketchSize = Integer.parseInt(cmd.getOptionValue("sketchSize","1000"));
        if (sketchSize<1){
            System.out.println("WARNING: minimum sketchsize  is 1 -> your value got increased to 1000 (default)");
            this.sketchSize = 1000;
        }

        this.prefixSize = Integer.parseInt(cmd.getOptionValue("prefixSize","8"));
        if (prefixSize<1){
            System.out.println("WARNING: minimum prefixSize is 1 -> your value got increased to 8 (default)");
            this.prefixSize = 21;
        }

        this.kmerSize = Integer.parseInt(cmd.getOptionValue("kmerSize","21"));
        if (kmerSize<1){
            System.out.println("WARNING: minimum kmerSize is 1 -> your value got increased to 21 (default)");
            this.kmerSize = 21;
        }
        if (kmerSize>32){
            System.out.println("WARNING: maximum kmerSize is 32 -> your value got decreased to 21 (default)");
            this.kmerSize = 21;
        }

        if (cmd.hasOption("hashFunction")){
            hashFunction = Integer.parseInt(cmd.getOptionValue("hashFunction"));
        } else{
            calculateHashFunction();
        }

        this.cores = Integer.parseInt(cmd.getOptionValue("parellism","4"));
        if (cores<1){
            System.out.println("WARNING: minimum number of threads is 1 -> your value got increased to 4 (default)");
            this.cores = 4;
        }

        this.seed = Integer.parseInt(cmd.getOptionValue("seed","42"));
        if (seed<1){
            System.out.println("WARNING: minimum seed is 1 -> your value got increased to 42 (default)");
            this.seed = 42;
        }

        this.lowKThreshold = Float.parseFloat(cmd.getOptionValue("lowKThreshold","0.01"));
        if (lowKThreshold>1 || lowKThreshold<0){
            System.out.println("WARNING: lowK-threshold has to be a probability (between 0 and 1) -> your value has been set to default(0.1)");
            this.lowKThreshold=0.01f;
        }

        this.outputFile = cmd.getOptionValue("output","test");

        sequences = new SynchronizedList<>((ArrayList<String>)sequenceFiles.clone());
    }
    public void calculateHashFunction(){
        this.hashFunction = (int) this.kmerSize/16;
        System.out.println("estimated hash function: "+this.hashFunction);
    }

    public void calculateOptimalBloomFilter(float p, int numOfItems){
        bloomFilterSize = (int)Math.ceil((numOfItems * Math.log(p)) / Math.log(1 / Math.pow(2, Math.log(2))));
        bloomFilterHashes = (int) Math.round((this.bloomFilterSize / (float)numOfItems) * Math.log(2));
        bloomFilterSize = (int)Math.ceil(Math.log(bloomFilterSize)/Math.log(2));

        System.out.println("Bloom Filter:");
        System.out.println("  estimated size:                     "+this.bloomFilterSize);
        System.out.println("  estimated number of hash functions: "+this.bloomFilterHashes);
    }
    public float calculatePvalueBloomFilter(int numOfItems){
        return (float)Math.pow(1 - Math.exp(-this.bloomFilterHashes / (Math.pow(2,this.bloomFilterSize) / numOfItems)), this.bloomFilterHashes);
    }


}
