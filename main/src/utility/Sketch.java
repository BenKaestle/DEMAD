package utility;

public class Sketch {
    private long[] hashes;
    private String header;
    private int genome_size;

    public Sketch(long[] hashes, String header, int genome_size) {
        this.hashes = hashes;
        this.header = header;
        this.genome_size = genome_size;
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
}
