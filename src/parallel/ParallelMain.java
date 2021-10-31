package parallel;

import qut.Sigma70Consensus;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class ParallelMain {
    private final HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();

    // rename this method to main to run genbank records in parallel stream
    public static void main1(String[] args) throws IOException {
        long startTime = System.nanoTime();
        ParallelStreaming main = new ParallelStreaming();
        main.run("referenceGenes.list", "Ecoli");

        long runTime = System.nanoTime() - startTime;

        for (Map.Entry<String, Sigma70Consensus> entry : main.consensus.entrySet()) {
            System.out.println(entry.getKey() + " " + entry.getValue());
        }

        System.out.println("Total runtime: " + runTime / 1000000000f + " seconds");
    }

    // rename this method to main to run parallel threads for reference genes
    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        long startTime = System.nanoTime();
        ParallelMain main = new ParallelMain();
        Thread parallelThread1 = new ParallelRunner(2, 0, main.consensus);
        Thread parallelThread2 = new ParallelRunner(2, 1, main.consensus);

        parallelThread1.start();
        parallelThread2.start();

        long runTime = System.nanoTime() - startTime;

        for (Map.Entry<String, Sigma70Consensus> entry : main.consensus.entrySet()) {
            System.out.println(entry.getKey() + " " + entry.getValue());
        }

        System.out.println("Total runtime: " + runTime / 1000000000f + " seconds");
    }
}
