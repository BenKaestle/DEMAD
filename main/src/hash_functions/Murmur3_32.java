package hash_functions;

/*
 *  Murmur3_32.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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
/**
 *  The MurmurHash3 algorithm was created by Austin Appleby and placed in the public domain.
 *  This java port was authored by Yonik Seeley and also placed into the public domain.
 *  The author hereby disclaims copyright to this source code.
 **/
public class Murmur3_32 implements HashFunction {


    int kmerSize;
    int seed;

    public Murmur3_32(int kmerSize, int seed) {
        this.kmerSize = kmerSize;
        this.seed = seed;
    }

    @Override
    public long hash(String kmer) {
        return murmurhash3_x86_32(kmer, 0, kmerSize, seed);
    }

    /** Returns the MurmurHash3_x86_32 hash of the UTF-8 bytes of the String without actually encoding
     * the string to a temporary buffer.  This is more than 2x faster than hashing the result
     * of String.getBytes().
     */
    public static int murmurhash3_x86_32(CharSequence data, int offset, int len, int seed) {

        final int c1 = 0xcc9e2d51;
        final int c2 = 0x1b873593;

        int h1 = seed;

        int pos = offset;
        int end = offset + len;
        int k1 = 0;
        int k2 = 0;
        int shift = 0;
        int bits = 0;
        int nBytes = 0;   // length in UTF8 bytes


        while (pos < end) {
            int code = data.charAt(pos++);
            if (code < 0x80) {
                k2 = code;
                bits = 8;
            }
            else if (code < 0x800) {
                k2 = (0xC0 | (code >> 6))
                        | ((0x80 | (code & 0x3F)) << 8);
                bits = 16;
            }
            else if (code < 0xD800 || code > 0xDFFF || pos>=end) {
                // we check for pos>=end to encode an unpaired surrogate as 3 bytes.
                k2 = (0xE0 | (code >> 12))
                        | ((0x80 | ((code >> 6) & 0x3F)) << 8)
                        | ((0x80 | (code & 0x3F)) << 16);
                bits = 24;
            } else {
                // surrogate pair
                // int utf32 = pos < end ? (int) data.charAt(pos++) : 0;
                int utf32 = (int) data.charAt(pos++);
                utf32 = ((code - 0xD7C0) << 10) + (utf32 & 0x3FF);
                k2 = (0xff & (0xF0 | (utf32 >> 18)))
                        | ((0x80 | ((utf32 >> 12) & 0x3F))) << 8
                        | ((0x80 | ((utf32 >> 6) & 0x3F))) << 16
                        |  (0x80 | (utf32 & 0x3F)) << 24;
                bits = 32;
            }


            k1 |= k2 << shift;


            shift += bits;
            if (shift >= 32) {
                // mix after we have a complete word

                k1 *= c1;
                k1 = (k1 << 15) | (k1 >>> 17);  // ROTL32(k1,15);
                k1 *= c2;

                h1 ^= k1;
                h1 = (h1 << 13) | (h1 >>> 19);  // ROTL32(h1,13);
                h1 = h1*5+0xe6546b64;

                shift -= 32;
                // unfortunately, java won't let you shift 32 bits off, so we need to check for 0
                if (shift != 0) {
                    k1 = k2 >>> (bits-shift);   // bits used == bits - newshift
                } else {
                    k1 = 0;
                }
                nBytes += 4;
            }

        } // inner

        // handle tail
        if (shift > 0) {
            nBytes += shift >> 3;
            k1 *= c1;
            k1 = (k1 << 15) | (k1 >>> 17);  // ROTL32(k1,15);
            k1 *= c2;
            h1 ^= k1;
        }

        // finalization
        h1 ^= nBytes;

        // fmix(h1);
        h1 ^= h1 >>> 16;
        h1 *= 0x85ebca6b;
        h1 ^= h1 >>> 13;
        h1 *= 0xc2b2ae35;
        h1 ^= h1 >>> 16;

        return h1;
    }
}
