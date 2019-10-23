package utility;

import java.util.HashMap;
import java.util.Map;

public class SketchSet {
    private final HashMap<Long,Integer> sketches;
    private final int count;

    public SketchSet(int count) {
        this.sketches = new HashMap<>();
        this.count = count;
    }

    public void add (long sketch_hash){
        sketches.putIfAbsent(sketch_hash,0);
        sketches.put(sketch_hash,sketches.get(sketch_hash)+1);
    }
    public boolean contains (long sketch_hash){
        return (sketches.get(sketch_hash)==this.count);
    }
}
