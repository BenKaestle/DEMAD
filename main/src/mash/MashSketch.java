package mash;

import java.io.Serializable;

public class MashSketch implements Serializable {
    private long[] hashes;
    private String header;
    private String filename;
    private int kmerSize;
    private int hashFunction;
    private int genome_size;
    private int hashSeed;

    public MashSketch(long[] hashes, String header, String filename, int kmerSize, int hashFunction, int genome_size, int hashSeed) {
        this.hashes = hashes;
        this.header = header;
        this.filename = filename;
        this.kmerSize = kmerSize;
        this.hashFunction = hashFunction;
        this.genome_size = genome_size;
        this.hashSeed = hashSeed;
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

    public long[] getHashes() {
        return hashes;
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
        if(hashFunction==2) hashF = "Wang-hash";
        if(hashFunction==3) hashF = "hash_code";
        return  "filename:         "+ filename+"\n"+
                "header:           "+ header+"\n"+
                "hash function:    "+ hashF+"\n"+
                "         seed:    "+ hashSeed+"\n"+
                "kmer size:        "+ kmerSize+"\n"+
                "genome size:      "+ genome_size+"\n"+
                "filename:         "+ filename+"\n"+
                "Number of hashes: "+ hashes.length+"\n"
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
