package hash_functions;

public class Test_Hash implements hash_function {
    @Override
    public long hash(String kmer, long modulo) {
        return kmer.hashCode();
    }
}
