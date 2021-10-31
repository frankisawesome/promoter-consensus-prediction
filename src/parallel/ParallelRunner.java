package parallel;

import qut.Sigma70Consensus;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;



public class ParallelRunner extends Thread {
    private int totalThreads = 1;
    private int threadIndex = 0;
    private HashMap<String, Sigma70Consensus> consensus;

    public ParallelRunner (int totalThreads, int threadIndex, HashMap<String, Sigma70Consensus> consensus) {
        this.totalThreads = totalThreads;
        this.threadIndex = threadIndex;
        this.consensus = consensus;
    }

    @Override
    public void run() {
        try {
            StaticRunner.execute("referenceGenes.list", "Ecoli", this.totalThreads, this.threadIndex, this.consensus);
        } catch (Exception e) {
            System.out.println("Error" + e.toString());
        }
    }
}
