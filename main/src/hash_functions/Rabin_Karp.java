package hash_functions;

/*
 *  Rabin_Karp.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

public class Rabin_Karp{

    public static long hash(String kmer, long modulo, int a, int p) {
        long hash = 0;
        int k = kmer.length()-1;
        for (char c :kmer.toCharArray()){
            switch (c){
                case 'A':
                    hash += 1*Math.pow(a,k+p);
                    break;
                case 'T':
                    hash += 2*Math.pow(a,k+p);
                    break;
                case 'C':
                    hash += 3*Math.pow(a,k+p);
                    break;
                case 'G':
                    hash += 4*Math.pow(a,k+p);
                    break;
            }
            k--;
        }
        return hash%modulo;
    }
    public static long hashNext(char out, char in, long last,long modulo, int a, int p, int kmerSize){
        switch (out){
            case 'A':
                last -= 1*Math.pow(a,kmerSize-1+p)%modulo;
                break;
            case 'T':
                last -= 2*Math.pow(a,kmerSize-1+p)%modulo;
                break;
            case 'C':
                last -= 3*Math.pow(a,kmerSize-1+p)%modulo;
                break;
            case 'G':
                last -= 4*Math.pow(a,kmerSize-1+p)%modulo;
                break;
        }
        last *= a;
        switch (in){
            case 'A':
                last += 1*Math.pow(a,p);
                break;
            case 'T':
                last += 2*Math.pow(a,p);
                break;
            case 'C':
                last += 3*Math.pow(a,p);
                break;
            case 'G':
                last += 4*Math.pow(a,p);
                break;
        }
        return last%modulo;
    }
}
