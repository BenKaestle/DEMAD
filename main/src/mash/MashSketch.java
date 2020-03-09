package mash;

import java.io.Serializable;

/*
 *  MashSketch.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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


    /**
     * to string function for command line printing
     * @return
     */
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
