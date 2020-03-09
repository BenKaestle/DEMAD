package utility;

import java.util.BitSet;

/*
 *  BloomFilter.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

public class BloomFilter {
    private static final int MAX_HASHES = 8;
    private static final long[] byteTable;
    private static final long HSTART = 0xBB40E64DA205B064L;
    private static final long HMULT = 7664345821815920749L;

    static {
        byteTable = new long[256 * MAX_HASHES];
        long h = 0x544B2FBACAAF1684L;
        for (int i = 0; i < byteTable.length; i++) {
            for (int j = 0; j < 31; j++)
                h = (h >>> 7) ^ h; h = (h << 11) ^ h; h = (h >>> 10) ^ h;
            byteTable[i] = h;
        }
    }

    /**
     * used hash function number hcNo for a input element s
     * @param s
     * @param hcNo
     * @return
     */
    private long hashCode(String s, int hcNo) {
        long h = HSTART;
        final long hmult = HMULT;
        final long[] ht = byteTable;
        int startIx = 256 * hcNo;
        for (int len = s.length(), i = 0; i < len; i++) {
            char ch = s.charAt(i);
            h = (h * hmult) ^ ht[startIx + (ch & 0xff)];
            h = (h * hmult) ^ ht[startIx + ((ch >>> 8) & 0xff)];
        }
        return h;
    }
    private final BitSet data;
    private final int noHashes;
    private final int hashMask;

    /**
     * constructor for a Bloom filter of size  2^log2noBits with noHashes hash functions
     * @param log2noBits
     * @param noHashes
     */
    public BloomFilter(int log2noBits, int noHashes) {
        if (log2noBits < 1 || log2noBits > 31)
            throw new IllegalArgumentException("Invalid number of bits");
        if (noHashes < 1 || noHashes > MAX_HASHES)
            throw new IllegalArgumentException("Invalid number of hashes");

        this.data = new BitSet(1 << log2noBits);
        this.noHashes = noHashes;
        this.hashMask = (1 << log2noBits) - 1;
    }

    /**
     * adds a new element to the filter
     * @param s
     */
    public void add(String s) {
        for (int n = 0; n < noHashes; n++) {
            long hc = hashCode(s, n);
            int bitNo = (int) (hc) & this.hashMask;
            data.set(bitNo);
        }
    }

    /**
     * checks if the filter contains a certain element s
     * @param s
     * @return
     */
    public boolean contains(String s) {
        for (int n = 0; n < noHashes; n++) {
            long hc = hashCode(s, n);
            int bitNo = (int) (hc) & this.hashMask;
            if (!data.get(bitNo)) return false;
        }
        return true;
    }
}
