package mash;

import java.util.Locale;

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

    public void print(){
        System.out.println("header1: "+ this.header_1);
        System.out.println("header2: "+ this.header_2);
        System.out.println("jaccard: "+ this.jaccard_index);
        System.out.println("MashDis: "+ this.mash_distance);
        System.out.println("p_Value: "+ this.p_value);
    }
    public void printShort(){
        System.out.print(this.filePath1 + "\t" + this.filePath2+"\t");
        System.out.printf("%f", this.mash_distance);
        System.out.print("\t");
        System.out.printf("%f", this.p_value);
        System.out.print("\t"+sameHashes+"\n");
    }

    public String toString(){
        return this.filePath1 + "\t" + this.filePath2 + "\t" + String.format(Locale.ROOT,"%f", this.mash_distance) + "\t" + String.format(Locale.ROOT,"%f", this.p_value) + "\t" + this.sameHashes;
    }
}
