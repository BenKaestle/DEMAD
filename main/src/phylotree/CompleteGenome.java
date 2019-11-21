package phylotree;

public class CompleteGenome {
    public String downloadLink;
    public int year;
    public String name;
    public int taxId;
    public int taxId2;

    public CompleteGenome(String downloadLink, int year, String name, int taxId, int taxId2) {
        this.downloadLink = downloadLink;
        this.year = year;
        this.name = name;
        this.taxId = taxId;
        this.taxId2 = taxId2;
    }

    public String getDownloadLink() {
        return downloadLink;
    }
}
