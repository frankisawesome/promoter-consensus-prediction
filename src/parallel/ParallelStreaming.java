package parallel;

import edu.au.jacobi.pattern.Match;
import qut.*;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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

//    protected static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
//        int upStreamDistance = 250;
//        if (gene.location < upStreamDistance)
//            upStreamDistance = gene.location - 1;
//
//        if (gene.strand == 1)
//            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
//        else {
//            byte[] result = new byte[upStreamDistance];
//            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
//            IntStream.rangeClosed(0, upStreamDistance).parallel()
//                    .map(i -> result[i] = complement[dna.bytes[reverseStart - i]])
//                    .boxed()
//                    .collect(Collectors.toList());
//
//            return new NucleotideSequence(result);
//        }
//    }

    public void run(String referenceFile, String dir) throws IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile, consensus);
        List<GenbankRecord> records = new LinkedList<>();

        // we might want to parallelise this for maxmimum performance
        records = ListGenbankFiles(dir).parallelStream()
                .map(Sequential::Parse)
                .collect(Collectors.toList());

        // have all dependencies ready to run in parallel stream
        final List<ComputationData> computationDataSet = records.parallelStream()
                .map(record -> {
                    List<ComputationData> compList= new ArrayList<>();
                     referenceGenes.forEach(referenceGene -> {
                         record.genes.forEach(gene -> {
                             compList.add(new ComputationData(gene, referenceGene, record.nucleotides, referenceGene.name));
                         });
                     });
                     return compList;
                })
                .flatMap(List::stream)
                .collect(Collectors.toList());

        ForkJoinPool customThreadPool = new ForkJoinPool(11);

        try {
            customThreadPool.submit(
                    () ->  {
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
                    }).get();
        } catch (Exception e) {
            throw new RuntimeException(e);
        } finally {
            customThreadPool.shutdown();
        }
    }
}
