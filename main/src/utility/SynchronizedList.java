package utility;

import java.util.ArrayList;

public class SynchronizedList<T> {
    ArrayList<T> list;

    public SynchronizedList(ArrayList<T> list) {
        this.list = list;
    }

    public synchronized T get(){
        if (list.isEmpty()) return null;
        return list.remove(0);
    }
    public synchronized boolean isEmpty(){
        return list.isEmpty();
    }
    public synchronized int size(){
        return list.size();
    }
}
