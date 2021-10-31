package qut;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;

import java.io.*;
import java.util.*;

public class Sequential {
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    protected static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
    protected static final Matrix BLOSUM_62 = BLOSUM62.Load();
    protected static byte[] complement = new byte['z'];

    static {
        complement['C'] = 'G';
        complement['c'] = 'g';
        complement['G'] = 'C';
        complement['g'] = 'c';
        complement['T'] = 'A';
        complement['t'] = 'a';
        complement['A'] = 'T';
        complement['a'] = 't';
    }


    protected static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException {
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

    protected static boolean Homologous(PeptideSequence A, PeptideSequence B) {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    protected static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location - 1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
        else {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i = 0; i < upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart - i]];
            return new NucleotideSequence(result);
        }
    }

    protected static Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    protected static void ProcessDir(List<String> list, File dir) {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    protected static List<String> ListGenbankFiles(String dir) {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    protected static GenbankRecord Parse(String file) {
        GenbankRecord record = new GenbankRecord();
        try {
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
            record.Parse(reader);
            reader.close();
        } catch (Exception e) {
            //
        }
        return record;
    }

    public void run(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes)
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

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

    public static void main(String[] args) throws IOException {
        long startTime = System.nanoTime();
        new Sequential().run("referenceGenes.list", "Ecoli");
        long runTime = System.nanoTime() - startTime;
        System.out.println("Total runtime: " + runTime / 1000000000f + " seconds");
    }
}
