package hash_functions;

public class Wang {

    public static long timeConsumingHash(String x){
        long value = 0;
        for(int i=0;i<x.length();i++){
            switch (x.charAt(x.length()-i-1)) {
                case 'A':
                    break;
                case 'C':
                    value += Math.pow(2, i);
                    break;
                case 'G':
                    value += 2 * Math.pow(2, i);
                    break;
                case 'T':
                    value += 3 * Math.pow(2, i);
                    break;
                default:
                    break;
            }
        }
        return hash64shift(value);
    }
    public static long hash(String x){
        return hash64shift(x.hashCode());
    }


    public static long hash64shift(long key)
    {
        key = (~key) + (key << 21); // key = (key << 21) - key - 1;
        key = key ^ (key >>> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >>> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >>> 28);
        key = key + (key << 31);
        return key;
    }

    
}
