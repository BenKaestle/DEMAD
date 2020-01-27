package utility;

import dashing.DashingSketch;
import mash.MashSketch;
import org.apache.commons.cli.*;

import java.io.File;
import java.util.ArrayList;

public class InputParameters {
    public SynchronizedList<String> sequences;
    public SynchronizedList<MashSketch> mashSketchesSynch;
    public SynchronizedList<DashingSketch> dashingSketchesSynch;
    public int sketchSize = 1000;
    public int kmerSize = 21;
    public int cores = 4;
    public boolean bloomFilter = false;
    public int hashFunction = 1;
    public int seed = 42;
    public int bloomFilterSize = 23; //max 31
    public int bloomFilterHashes = 7; //max 8
    public String type; //sketch, dist, info
    public MashSketch[] mashSketches;
    public DashingSketch[] dashingSketches;
    public String outputFile;
    public boolean tableOutput;
    public ArrayList<String> sequenceFiles;
    public float lowKThreshold;
    public boolean outputPhylip;
    public int prefixSize;
    public String alphabet;
    public int reduceAlphabet;
    private float bloomThreshold;
    private int bloomFilterGenomeSize;


    public void parseInput(String[] args, boolean mash) {
        if (args.length == 0) {
            System.out.println("Choose one of: sketch, dist, info");
            System.exit(1);
        }
        this.type = args[0].trim();
        if (!type.equals("sketch") && !type.equals("dist") && !type.equals("info")) {
            System.out.println("Choose one of: sketch, dist, info");
            System.exit(1);
        }
        Options options = new Options();


        Option input = null;

        if (mash) {
            if (this.type.equals("dist"))
                input = new Option("i", "input", true, "input file path (either *.msk file or all fasta files in directory)");
            else if (this.type.equals("sketch"))
                input = new Option("i", "input", true, "input file path (all fasta files in directory will get sketched)");
            else if (this.type.equals("info"))
                input = new Option("i", "input", true, "*.msk file of interest");
            input.setRequired(true);
            options.addOption(input);
        } else {
            if (this.type.equals("dist"))
                input = new Option("i", "input", true, "input file path (either *.dsk file or all fasta files in directory)");
            else if (this.type.equals("sketch"))
                input = new Option("i", "input", true, "input file path (all fasta files in directory will get sketched)");
            else if (this.type.equals("info"))
                input = new Option("i", "input", true, "*.dsk file of interest");
            input.setRequired(true);
            options.addOption(input);
        }

        Option output = new Option("o", "output", true, "output file(path)");
        if (this.type.equals("sketch")) output.setRequired(true);
        options.addOption(output);




        Option outputPhylip = new Option("P", "outputPhylip", false, " output the distance Matrix as a *.phylip file");
        options.addOption(outputPhylip);


        if(mash){
            Option tableOut = new Option("t", "table", false, "output a distance matrix (only mash dist)");
            options.addOption(tableOut);

            Option kmerSize = new Option("k", "kmerSize", true, "length of kmers (default = 21)");
            options.addOption(kmerSize);
        }else {
            Option tableOut = new Option("t", "table", false, "output a distance matrix (only jaccard index)");
            options.addOption(tableOut);

            Option kmerSize = new Option("k", "kmerSize", true, "length of kmers (default = 31)");
            options.addOption(kmerSize);
        }

        Option lowKThreshold = new Option("wk", "lowKThreshold", true, " Probability threshold for warning about low k-mer size. (0-1) (default = 0.01)");
        options.addOption(lowKThreshold);

        if(mash) {
            Option sketchSize = new Option("s", "sketchSize", true, "size of the generated sketches (default = 1000)");
            options.addOption(sketchSize);
        }else {
            Option prefixSize = new Option("ps", "prefixSize", true, "length of the prefix (log2 of register size) (default = 10)");
            options.addOption(prefixSize);
        }
        Option cores = new Option("p", "parellism", true, "this many threads will be spawned for processing (default = 4)");
        options.addOption(cores);


        Option hashFunction = new Option("hf", "hashFunction", true, "which hash-function should be used");
        options.addOption(hashFunction);

        Option seed = new Option("S", "seed", true, "seed used for hash function (default = 42)");
        options.addOption(seed);

        Option bloomFilter = new Option("b", "bloomFilter", true, "a bloom filter is used for sketching (optional add estimated genome size)");
        options.addOption(bloomFilter);

        Option bloomFilterParameters = new Option("bp", "bloomFilterParameters", true, "size (log2(#bits)) and number of hash functions of the Used BloomFilter");
        bloomFilterParameters.setArgs(2);
        options.addOption(bloomFilterParameters);

        Option bloomThreshold = new Option("wb", "bloomThreshold", true, " Probability threshold for warning about false positives in bloom filter. (0-1) (default = 0.01)");
        options.addOption(bloomThreshold);

        Option alphabet = new Option("a", "alphabet", true, "input belongs to the alphabet of DNA or AA (default DNA)");
        options.addOption(alphabet);

        Option reduceAA = new Option("r", "reduceAlphabetTo", true, "reduces the alphabet to one of {2,4,8,10,15} cardinalities(requires alphabet=AA)");
        options.addOption(reduceAA);

        Option help = new Option("h", "help", false, "shows help");
        options.addOption(help);

        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        formatter.setOptionComparator(null);
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            formatter.printHelp("utility-name", options);

            System.exit(1);
        }

        if (cmd.hasOption("help")) {
            formatter.printHelp("help", options);
            System.exit(1);
        }

        this.update(cmd);
        if (mash){
            this.updateMash(cmd);
        }else {
            this.updateDashing(cmd);
        }
    }

    private void updateMash(CommandLine cmd) {
        this.kmerSize = Integer.parseInt(cmd.getOptionValue("kmerSize", "21"));
        if (kmerSize < 1) {
            System.out.println("WARNING: minimum kmerSize is 1 -> your value got increased to 21 (default)");
            this.kmerSize = 21;
        }
        if (kmerSize > 32) {
            System.out.println("WARNING: maximum kmerSize is 32 -> your value got decreased to 21 (default)");
            this.kmerSize = 21;
        }
        sequenceFiles = new ArrayList<>();
        String input = cmd.getOptionValue("input");
        if (input.endsWith(".msk") && (this.type.equals("dist") || this.type.equals("info"))) {
            this.mashSketches = (MashSketch[]) WriteReadObject.readObjectFromFile(input);
            for (MashSketch sketch : mashSketches) {
                sequenceFiles.add(sketch.getFilename());
            }
        } else {
            File folder = new File(input);
            File[] listOfFiles = folder.listFiles();
            if (listOfFiles == null) {
                System.out.println("Input file/path: \"" + input + "\" not found");
                System.exit(1);
            }
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile() && (listOfFiles[i].getPath().endsWith(".fasta") || listOfFiles[i].getPath().endsWith(".fna"))) {
                    sequenceFiles.add(listOfFiles[i].getPath());
                }
            }
            if (sequenceFiles.size() < 2) {
                System.out.println("There was only one or none fasta files at given directory. Please retry with valid input.");
                System.exit(1);
            }
        }
        if (this.mashSketches == null && this.type.equals("info")) {
            System.out.println("To get information on a sketch, you need to input a valid *.msk file");
            System.exit(1);
        }
        sequences = new SynchronizedList<>((ArrayList<String>) sequenceFiles.clone());

        if (cmd.hasOption("hashFunction")) {
            hashFunction = Integer.parseInt(cmd.getOptionValue("hashFunction"));
        } else {
            calculateHashFunction();
        }
    }

    private void updateDashing(CommandLine cmd) {
        this.kmerSize = Integer.parseInt(cmd.getOptionValue("kmerSize", "31"));
        if (kmerSize < 1) {
            System.out.println("WARNING: minimum kmerSize is 1 -> your value got increased to 31 (default)");
            this.kmerSize = 31;
        }
        if (kmerSize > 32) {
            System.out.println("WARNING: maximum kmerSize is 32 -> your value got decreased to 31 (default)");
            this.kmerSize = 31;
        }
        sequenceFiles = new ArrayList<>();
        String input = cmd.getOptionValue("input");
        if (input.endsWith(".dsk") && (this.type.equals("dist") || this.type.equals("info"))) {
            this.dashingSketches = (DashingSketch[]) WriteReadObject.readObjectFromFile(input);
            for (DashingSketch sketch : dashingSketches) {
                sequenceFiles.add(sketch.getFilename());
            }
        } else {
            File folder = new File(input);
            File[] listOfFiles = folder.listFiles();
            if (listOfFiles == null) {
                System.out.println("Input file/path: \"" + input + "\" not found");
                System.exit(1);
            }
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile() && (listOfFiles[i].getPath().endsWith(".fasta") || listOfFiles[i].getPath().endsWith(".fna"))) {
                    sequenceFiles.add(listOfFiles[i].getPath());
                }
            }
            if (sequenceFiles.size() < 2) {
                System.out.println("There was only one or none fasta files at given directory. Please retry with valid input.");
                System.exit(1);
            }
        }
        if (this.dashingSketches == null && this.type.equals("info")) {
            System.out.println("To get information on a sketch, you need to input a valid *.dsk file");
            System.exit(1);
        }
        sequences = new SynchronizedList<>((ArrayList<String>) sequenceFiles.clone());

        hashFunction = Integer.parseInt(cmd.getOptionValue("hashFunction", "2"));
    }

    private void update(CommandLine cmd) {
        this.tableOutput = cmd.hasOption("table");
        this.bloomFilter = cmd.hasOption("bloomFilter");
        this.outputPhylip = cmd.hasOption("outputPhylip");
        this.alphabet = cmd.getOptionValue("alphabet", "DNA").trim();
        if (this.alphabet.equals("DNA") && this.alphabet.equals("AA")) {
            System.out.println("WARNING: the given alphabet " + this.alphabet + " does not exist. Please choose one of: \"DNA\", \"AA\"");
            System.exit(1);
        }
        this.reduceAlphabet = Integer.parseInt(cmd.getOptionValue("reduceAlphabetTo", "-1"));
        if (this.reduceAlphabet != -1 && !this.alphabet.equals("AA")) {
            System.out.println("WARNING: reducing the alphabet with -r requires an alphabet consisting of amino acids. Please add \"-a AA\"");
            System.exit(1);
        }
        if (this.alphabet.equals("AA")) {
            switch (this.reduceAlphabet) {
                case 2:
                    this.reduceAlphabet = 4;
                    break;
                case 4:
                    this.reduceAlphabet = 3;
                    break;
                case 8:
                    this.reduceAlphabet = 2;
                    break;
                case 10:
                    this.reduceAlphabet = 1;
                    break;
                case 15:
                    this.reduceAlphabet = 0;
                    break;
                default:
                    System.out.println("WARNING: " + this.reduceAlphabet + " is no valid reduction cardinality. Please choose one of {2,4,8,10,15}");
                    System.exit(1);
                    break;
            }
        }
        if (this.bloomFilter) {
            this.bloomFilterGenomeSize = Integer.parseInt(cmd.getOptionValue("bloomFilter", "1000000"));
            if (cmd.hasOption("bloomFilterParameters")) {
                this.bloomFilterSize = Integer.parseInt(cmd.getOptionValues("bloomFilterParameters")[0]);
                this.bloomFilterHashes = Integer.parseInt(cmd.getOptionValues("bloomFilterParameters")[1]);

            } else {
                calculateOptimalBloomFilter(bloomThreshold, this.bloomFilterGenomeSize);
            }
            if (this.bloomFilterHashes > 8) {
                System.out.println("WARNING: maximum number of bloom filter hashes is 8 -> your value got reduced to 8");
                this.bloomFilterHashes = 8;
            }
            if (this.bloomFilterHashes < 1) {
                System.out.println("WARNING: minimum number of bloom filter hashes is 1 -> your value got increased to 1");
                this.bloomFilterHashes = 1;
            }
            if (this.bloomFilterSize < 1) {
                System.out.println("WARNING: minimum size of bloom filter is 1 -> your value got increased to 1");
                this.bloomFilterSize = 1;
            }
            if (this.bloomFilterSize > 31) {
                System.out.println("WARNING: maximum size of bloom filter is 31 -> your value got reduced to 31");
                this.bloomFilterSize = 31;
            }
            if (calculatePvalueBloomFilter(this.bloomFilterSize) > this.bloomThreshold) {
                System.out.println("WARNING: possibility of false positive in bloom filter higher than" + bloomThreshold);
            }

        }
        this.sketchSize = Integer.parseInt(cmd.getOptionValue("sketchSize", "1000"));
        if (sketchSize < 1) {
            System.out.println("WARNING: minimum sketchsize  is 1 -> your value got increased to 1000 (default)");
            this.sketchSize = 1000;
        }

        this.prefixSize = Integer.parseInt(cmd.getOptionValue("prefixSize", "10"));
        if (prefixSize < 1) {
            System.out.println("WARNING: minimum prefixSize is 1 -> your value got increased to 10 (default)");
            this.prefixSize = 10;
        }

        this.cores = Integer.parseInt(cmd.getOptionValue("parellism", "4"));
        if (cores < 1) {
            System.out.println("WARNING: minimum number of threads is 1 -> your value got increased to 4 (default)");
            this.cores = 4;
        }

        this.seed = Integer.parseInt(cmd.getOptionValue("seed", "42"));
        if (seed < 1) {
            System.out.println("WARNING: minimum seed is 1 -> your value got increased to 42 (default)");
            this.seed = 42;
        }
        this.lowKThreshold = Float.parseFloat(cmd.getOptionValue("lowKThreshold", "0.01"));
        if (lowKThreshold > 1 || lowKThreshold < 0) {
            System.out.println("WARNING: lowK-threshold has to be a probability (between 0 and 1) -> your value has been set to default(0.01)");
            this.lowKThreshold = 0.01f;
        }

        this.bloomThreshold = Float.parseFloat(cmd.getOptionValue("bloomThreshold", "0.01"));
        if (bloomThreshold > 1 || bloomThreshold < 0) {
            System.out.println("WARNING: bloom-Threshold has to be a probability (between 0 and 1) -> your value has been set to default(0.01)");
            this.bloomThreshold = 0.01f;
        }

        this.outputFile = cmd.getOptionValue("output", "test");
    }

    public void calculateHashFunction() {
        this.hashFunction = (int) this.kmerSize / 16;
        System.out.println("estimated hash function: " + this.hashFunction);
    }

    public void calculateOptimalBloomFilter(float p, int numOfItems) {
        bloomFilterSize = (int) Math.ceil((numOfItems * Math.log(p)) / Math.log(1 / Math.pow(2, Math.log(2))));
        bloomFilterHashes = (int) Math.round((this.bloomFilterSize / (float) numOfItems) * Math.log(2));
        bloomFilterSize = (int) Math.ceil(Math.log(bloomFilterSize) / Math.log(2));

        System.out.println("Bloom Filter:");
        System.out.println("  estimated size:                     " + this.bloomFilterSize);
        System.out.println("  estimated number of hash functions: " + this.bloomFilterHashes);
    }

    public float calculatePvalueBloomFilter(int numOfItems) {
        return (float) Math.pow(1 - Math.exp(-this.bloomFilterHashes / (Math.pow(2, this.bloomFilterSize) / numOfItems)), this.bloomFilterHashes);
    }


}
