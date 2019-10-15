package hash_functions;

public class Rabin_Karp implements hash_function{

    // c1 * a^k-1+p+...+ck * a^0+p mod x

    private int a;
    private int p;

    public Rabin_Karp(int a, int p) {
        this.a = a;
        this.p = p;
    }

    @Override
    public long hash(String kmer, long modulo) {
        long hash = 0;
        int k = kmer.length()-1;
        for (char c :kmer.toCharArray()){
            switch (c){
                case 'A':
                    hash += 1*Math.pow(this.a,k+p)%modulo;
                    break;
                case 'T':
                    hash += 2*Math.pow(this.a,k+p)%modulo;
                    break;
                case 'C':
                    hash += 3*Math.pow(this.a,k+p)%modulo;
                    break;
                case 'G':
                    hash += 4*Math.pow(this.a,k+p)%modulo;
                    break;
            }
            k--;
        }
        return hash%modulo;
    }
}
