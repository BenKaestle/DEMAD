package dashing;

import java.util.Locale;

public class DashingDistance {
    private String header_1;
    private String header_2;
    private float jaccard_index;
    private String filePath1;
    private String filePath2;
    private float mash_distance;

    public DashingDistance(String header_1, String header_2, float jaccard_index, String filePath1, String filePath2, float mash_distance) {
        this.header_1 = header_1;
        this.header_2 = header_2;
        this.jaccard_index = jaccard_index;
        this.filePath1 = filePath1;
        this.filePath2 = filePath2;
        this.mash_distance = mash_distance;
    }

    public DashingDistance(DashingDistance dashingDistance){
        this.header_1 = dashingDistance.header_2;
        this.header_2 = dashingDistance.header_1;
        this.jaccard_index = dashingDistance.jaccard_index;
        this.filePath1 = dashingDistance.filePath2;
        this.filePath2 = dashingDistance.filePath1;
        this.mash_distance = dashingDistance.mash_distance;
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

    public float getMash_distance() {
        return mash_distance;
    }

    public void print(){
        System.out.println("header1: "+ this.header_1);
        System.out.println("header2: "+ this.header_2);
        System.out.println("jaccard: "+ this.jaccard_index);
    }

    public String toString(){
        return this.filePath1 + "\t" + this.filePath2 + "\t" + String.format(Locale.ROOT,"%f", this.jaccard_index) + "\t" + String.format(Locale.ROOT,"%f", this.mash_distance);
    }

    public void printShort(){
        System.out.print(this.filePath1 + "\t" + this.filePath2+"\t");
        System.out.printf("%f", this.jaccard_index);
        System.out.print("\n");
    }
}
