package hash_functions;

/**
 * @author Ben on 21.01.2020
 * @project k-mer_based_distance_algo
 */
public class JavaHashCode implements HashFunction {
    @Override
    public long hash(String kmer) {
        return kmer.hashCode();
    }
}
