package parallel;

import edu.au.jacobi.pattern.Match;
import qut.*;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

public class ParallelStreaming extends Sequential {
    private class ComputationData {
        public final Gene gene;
        public final Gene referenceGene;
        public final NucleotideSequence nucleotides;
        public final String name;

        public ComputationData(Gene gene, Gene referenceGene, NucleotideSequence nucleotides, String name) {
            this.gene = gene;
            this.referenceGene = referenceGene;
            this.nucleotides = nucleotides;
            this.name = name;
        }
    }

    public HashMap<String, Sigma70Consensus> consensus = new HashMap<>();

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

    public void run(String referenceFile, String dir) throws IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile, consensus);
        List<ComputationData> computationDataSet = new LinkedList<>();
        List<GenbankRecord> records = new LinkedList<>();

        // we might want to parallelise this for maxmimum performance
        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            records.add(Parse(filename));
        }

        // have all dependencies ready to run in parallel stream
        for (GenbankRecord record : records) {
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    computationDataSet.add(new ComputationData(gene, referenceGene, record.nucleotides, referenceGene.name));
                }
            }
        }

        computationDataSet.parallelStream()
                .filter(computationData -> Homologous(computationData.gene.sequence, computationData.referenceGene.sequence))
                .collect(Collectors.toList())
                .forEach(computationData -> {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(computationData.nucleotides, computationData.gene);
                    Match prediction = PredictPromoter(upStreamRegion);

                    if (prediction != null) {
                        consensus.get(computationData.name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                    }
                });
    }
}
