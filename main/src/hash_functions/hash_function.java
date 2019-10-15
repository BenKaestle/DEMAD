package hash_functions;

public interface hash_function {
    public long hash(String kmer, long modulo);
}
