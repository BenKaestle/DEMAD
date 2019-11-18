package dashing;

import java.io.Serializable;

public class DashingSketch implements Serializable {
    private int[] register;
    private String header;
    private String filename;
    private int kmerSize;
    private int hashFunction;
    private int genome_size;
    private int hashSeed;
    private float harmonicMean;
    private boolean noHarmonicMean = true;

    public DashingSketch(int[] register, String header, String filename, int kmerSize, int hashFunction, int genome_size, int hashSeed) {
        this.register = register;
        this.header = header;
        this.filename = filename;
        this.kmerSize = kmerSize;
        this.hashFunction = hashFunction;
        this.genome_size = genome_size;
        this.hashSeed = hashSeed;
    }

    public boolean hasNoHarmonicMean() {
        return noHarmonicMean;
    }

    public void setNoHarmonicMean(boolean noHarmonicMean) {
        this.noHarmonicMean = noHarmonicMean;
    }

    public float getHarmonicMean() {
        return harmonicMean;
    }

    public void setHarmonicMean(float harmonicMean) {
        this.harmonicMean = harmonicMean;
    }

    public String getFilename() {
        return filename;
    }

    public int getKmerSize() {
        return kmerSize;
    }

    public int getHashFunction() {
        return hashFunction;
    }

    public int[] getRegister() {
        return register;
    }

    public int getGenome_size() {
        return genome_size;
    }

    public String getHeader() {
        return header;
    }



    public String toString() {
        String hashF ="";
        if(hashFunction==0) hashF = "murmur3_x86_32";
        if(hashFunction==1) hashF = "murmur3_x64_128(first half)";
        return  "filename:         "+ filename+"\n"+
                "header:           "+ header+"\n"+
                "hash function:    "+ hashF+"\n"+
                "         seed:    "+ hashSeed+"\n"+
                "kmer size:        "+ kmerSize+"\n"+
                "genome size:      "+ genome_size+"\n"+
                "filename:         "+ filename+"\n"+
                "size of register: "+ register.length+"\n"
//                +
//                hashes[0]+"\n"+
//                hashes[2]+"\n"+
//                hashes[3]+"\n"+
//                hashes[4]+"\n"+
//                hashes[5]+"\n"+
//                hashes[6]+"\n"
                ;
    }
}
