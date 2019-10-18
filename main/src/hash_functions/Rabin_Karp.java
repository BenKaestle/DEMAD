package hash_functions;

public class Rabin_Karp{

    // c1 * a^k-1+p+...+ck * a^0+p mod x

//    private int a;
//    private int p;
//
//    public Rabin_Karp(int a, int p) {
//        this.a = a;
//        this.p = p;
//    }

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
