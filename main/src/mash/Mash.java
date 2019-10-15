package mash;

import hash_functions.Rabin_Karp;
import utility.FastaParser;

import java.util.Arrays;
import java.util.concurrent.*;

class Task implements Runnable
{
    private String name;
    private CountDownLatch latch;

    public Task(String name, CountDownLatch latch)
    {
        this.name = name;
        this.latch = latch;
    }

    public String getName() {
        return name;
    }

    @Override
    public void run()
    {
        try
        {
            Long duration = (long) (Math.random() * 50);
            System.out.println("Doing a task during : " + name);
            TimeUnit.SECONDS.sleep(duration);
            latch.countDown();
        }
        catch (InterruptedException e)
        {
            e.printStackTrace();
        }
    }
}



public class Mash {
    public static void main(String[] args){
        Rabin_Karp rabin_karp = new Rabin_Karp(11,0);
        System.out.println(rabin_karp.hash("ATCAT", (long) Math.pow(2,32)));
        System.out.println(rabin_karp.hash("ATGTGCGATTCGATTCGATTA", (long) Math.pow(2,32)));
        System.out.println(rabin_karp.hash("GGGGGGGGGGGGGGGGGGGGG", (long) Math.pow(2,32)));
        System.out.println( (long) Math.pow(2,64));

    }
    public static int[] mash_sketch (String sequence, int sketchSize, int kmerSize, int cores){
        String[] split = new String[cores];
        split = splitBy(cores, sequence, kmerSize);


        int[] sketch = new int[sketchSize];
        Arrays.fill(sketch, Integer.MAX_VALUE);
        for (int i=0; i<=sequence.length()-kmerSize; i++){
            String kmer = canonical_kmer(sequence.substring(i,i+kmerSize));
            if (kmer.hashCode()<sketch[0]){
                boolean found = false;
                for (int j=1; j<sketchSize;j++){
                    if (kmer.hashCode()<sketch[j]){
                        sketch[j-1] = sketch[j];
                    }
                    else{
                        sketch[j-1] = kmer.hashCode();
                        found = true;
                        break;
                    }
                }
                if (!found){
                    sketch[sketchSize-1]=kmer.hashCode();
                }
            }

        }
        return sketch;
    }

    public static String[] splitBy(int cores, String sequence, int kmersize){
        String[] split = new String[cores];
        float splitHere= sequence.length()/(float)cores;
        int addDown = (kmersize-1)/2;
        int addUp = (int) Math.ceil((kmersize-1)/2);

        for (int i=0;i<cores;i++){
            int start=(int)(splitHere*i)-addDown;
            int end =(int)(splitHere*(i+1))+addUp;
            if (start<0)
                start=0;
            if (end>sequence.length())
                end = sequence.length();
            split[i] = sequence.substring(start,end);
        }
        return split;
    }

    public static void runInParallel(int cores){
        int threads=10;
        ExecutorService pool = Executors.newFixedThreadPool(cores);
        CountDownLatch latch = new CountDownLatch(threads);
        for (int i = 0; i < threads; i++)
        {
            Task task = new Task("Task " + i, latch);
            System.out.println("A new task has been added : " + task.getName());
            pool.submit(task);
        }
        try {
            latch.await();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("test");
        pool.shutdown();
    }


    //todo is it faster to just reverse the whole sequence and take the smaller one?
    public static String canonical_kmer (String kmer){
        String reverse = "";
        for (int i=kmer.length()-1; i>=0;i--){
            switch (kmer.charAt(i)){
                case 'A':
                    reverse+="T";
                    break;
                case 'T':
                    reverse+="A";
                    break;
                case 'G':
                    reverse+="C";
                    break;
                case 'C':
                    reverse+="G";
                    break;
            }
        }
        return kmer.compareTo(reverse) > 0 ? reverse : kmer;
    }
}
