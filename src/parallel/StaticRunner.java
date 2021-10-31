package parallel;

import edu.au.jacobi.pattern.Match;
import qut.*;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class StaticRunner extends Sequential {
    private static int[] findStartAndFinishIndex(int count, int totalThreads, int threadIndex) {
        int start = count / totalThreads * threadIndex;
        int end = -1;
        if (totalThreads == threadIndex + 1) {
            end = count - 1;
        } else {
            end = count / totalThreads * (threadIndex + 1);
        }
        return new int[]{start, end};
    }

    private static List<Gene> ParseReferenceGenes(String referenceFile, HashMap<String, Sigma70Consensus> consensus) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true) {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    public static void execute(String referenceFile, String dir, int total, int index, HashMap<String, Sigma70Consensus> consensus) throws IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile, consensus);

        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);

            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                int[] startFinish = findStartAndFinishIndex(record.genes.size(), total, index);

                for (int i = startFinish[0]; i <= startFinish[1]; i++) {
                    Gene gene = record.genes.get(i);

                    if (Homologous(gene.sequence, referenceGene.sequence)) {
                        NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                        Match prediction = PredictPromoter(upStreamRegion);

                        if (prediction != null) {
                            consensus.get(referenceGene.name).addMatch(prediction);
                            consensus.get("all").addMatch(prediction);
                        }
                    }
                }
            }
        }
    }
}
