package hash_functions;

/*
 *  Wang.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

public class Wang implements HashFunction{

    public long timeConsumingHash(String x){
        long value = 0;
        for(int i=0;i<x.length();i++){
            switch (x.charAt(x.length()-i-1)) {
                case 'A':
                    break;
                case 'C':
                    value += Math.pow(2, i);
                    break;
                case 'G':
                    value += 2 * Math.pow(2, i);
                    break;
                case 'T':
                    value += 3 * Math.pow(2, i);
                    break;
                default:
                    break;
            }
        }
        return hash64shift(value);
    }


    public long hash(String kmer){
        return hash64shift(kmer.hashCode());
    }


    public  long hash64shift(long key)
    {
        key = (~key) + (key << 21); // key = (key << 21) - key - 1;
        key = key ^ (key >>> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >>> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >>> 28);
        key = key + (key << 31);
        return key;
    }

    
}
