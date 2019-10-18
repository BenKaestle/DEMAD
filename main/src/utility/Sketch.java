package utility;

public class Sketch {
    private long[] hashes;
    private String header;

    public Sketch(long[] hashes, String header) {
        this.hashes = hashes;
        this.header = header;
    }

    public long[] getHashes() {
        return hashes;
    }

    public void setHashes(long[] hashes) {
        this.hashes = hashes;
    }

    public String getHeader() {
        return header;
    }

    public void setHeader(String header) {
        this.header = header;
    }
}
