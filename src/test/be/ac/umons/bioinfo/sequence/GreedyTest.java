package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.sequence.SequenceAlignment;


import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static junit.framework.TestCase.assertEquals;


/**
 * Created by aline on 8/04/16.
 */
public class GreedyTest {

    @Test
    @Ignore
    public void greedyTroisSequencesTest() {

        Greedy greed = new Greedy();
        String path2 = "src/test/be/ac/umons/bioinfo/sequence/Collections/test1/collection1.fasta";
        File file = new File(path2);

        FastaReader fasta = new FastaReader();

        try
        {
            List<Sequence> list = fasta.readFromFile(file);
            List<SequenceAlignment> result = greed.computePath(list, 1, -1, -2);
            Iterator<SequenceAlignment> iterator = result.iterator();

            while (iterator.hasNext()) {
                SequenceAlignment alignment = iterator.next();
                System.out.println(alignment.s1);
                System.out.println(alignment.s2);

                System.out.println("****************");
            }


        } catch (IOException e) {
            System.out.println("ça ne fonctionne pas ");
            System.out.println(e);
        }
    }

    @Test
    public void testArcGenerationNumber()
    {
        List<Sequence> sequences = Arrays.asList(
                new Sequence("actttacg"),
                new Sequence("ttgcacgat"),
                new Sequence("ttgcg"),
                new Sequence("ggaatctgcgagtta"),
                new Sequence("tact"),
                new Sequence("gaccgat")
        );

        List<Arc> arcs = Greedy.generateArcs(sequences, 1, -1, -2);

        assertEquals(8*(6*5)/2, arcs.size());
    }

    @Test
    public void testArcFilteringNumber()
    {
        List<Sequence> sequences = Arrays.asList(
                new Sequence("actttacg"),
                new Sequence("ttgcacgat"),
                new Sequence("ttgcg"),
                new Sequence("ggaatctgcgagtta"),
                new Sequence("tact"),
                new Sequence("gaccgat")
        );

        List<Arc> arcs = Greedy.generateArcs(sequences, 1, -1, -2);
        Set<Arc> filtered = Greedy.filterArcs(arcs, sequences);

        assertEquals(sequences.size()-1, filtered.size());
    }

    @Test
    public void testArcFilteringNumber2()
    {
        List<Sequence> sequences = Arrays.asList(
            new Sequence("catagtc"),
            new Sequence("taactat"),
            new Sequence("agactatcc")
        );

        List<Arc> arcs = Greedy.generateArcs(sequences, 1, -1, -2);
        Set<Arc> filtered = Greedy.filterArcs(arcs, sequences);

        assertEquals(sequences.size()-1, filtered.size());
    }
}
