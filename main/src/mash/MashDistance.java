package mash;

import java.util.Locale;

/*
 *  MashDistance.java Copyright (C) 2020 Algorithms in Bioinformatics, University of Tuebingen
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

public class MashDistance {
    private String header_1;
    private String header_2;
    private float jaccard_index;
    private double p_value;
    private float mash_distance;
    private String filePath1;
    private String filePath2;
    private int sameHashes;

    public MashDistance(String header_1, String header_2, float jaccard_index, double p_value, float mash_distance, String filePath1, String filePath2, int sameHashes) {
        this.header_1 = header_1;
        this.header_2 = header_2;
        this.jaccard_index = jaccard_index;
        this.p_value = p_value;
        this.mash_distance = mash_distance;
        this.filePath1 = filePath1;
        this.filePath2 = filePath2;
        this.sameHashes = sameHashes;
    }

    public MashDistance(MashDistance mashDistance){
        this.header_1 = mashDistance.header_2;
        this.header_2 = mashDistance.header_1;
        this.jaccard_index = mashDistance.jaccard_index;
        this.p_value = mashDistance.p_value;
        this.mash_distance = mashDistance.mash_distance;
        this.filePath1 = mashDistance.filePath2;
        this.filePath2 = mashDistance.filePath1;
        this.sameHashes = mashDistance.sameHashes;
    }

    public int getSameHashes() {
        return sameHashes;
    }

    public String getFilePath1() {
        return filePath1;
    }

    public String getFilePath2() {
        return filePath2;
    }

    public String getHeader_1() {
        return header_1;
    }

    public String getHeader_2() {
        return header_2;
    }

    public float getJaccard_index() {
        return jaccard_index;
    }

    public double getP_value() {
        return p_value;
    }

    public float getMash_distance() {
        return mash_distance;
    }

    /**
     * to string function for command line printing
     * @return
     */
    public String toString(){
        return this.filePath1 + "\t" + this.filePath2 + "\t" + String.format(Locale.ROOT,"%f", this.mash_distance) + "\t" + String.format(Locale.ROOT,"%f", this.p_value) + "\t" + this.sameHashes;
    }
}
