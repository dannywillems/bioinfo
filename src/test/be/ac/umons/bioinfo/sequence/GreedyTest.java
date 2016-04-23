package be.ac.umons.bioinfo.sequence;


import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;


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
            List<SequenceAlignment> result = greed.greedy(list, 1, -1, -2);
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

    @Test
    public void filteredArcsMustFormHamiltonianPath()
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

        Set<Sequence> startSequences = new HashSet<>();
        Set<Sequence> endSequences = new HashSet<>();

        for(Arc arc : filtered)
        {
            assertFalse(startSequences.contains(arc.start));
            assertFalse(endSequences.contains(arc.end));

            startSequences.add(arc.start);
            endSequences.add(arc.end);
        }

        assertEquals(sequences.size()-1, startSequences.size());
        assertEquals(sequences.size()-1, endSequences.size());
    }

    @Test
    public void filteredArcsMustFormHamiltonianPath2()
    {
        List<Sequence> sequences = Arrays.asList(
                new Sequence("catagtc"),
                new Sequence("taactat"),
                new Sequence("agactatcc")
        );

        List<Arc> arcs = Greedy.generateArcs(sequences, 1, -1, -2);
        Set<Arc> filtered = Greedy.filterArcs(arcs, sequences);

        Set<Sequence> startSequences = new HashSet<>();
        Set<Sequence> endSequences = new HashSet<>();

        for(Arc arc : filtered)
        {
            assertFalse(startSequences.contains(arc.start));
            assertFalse(endSequences.contains(arc.end));

            startSequences.add(arc.start);
            endSequences.add(arc.end);
        }

        assertEquals(sequences.size()-1, startSequences.size());
        assertEquals(sequences.size()-1, endSequences.size());
    }

    @Test
    public void orderMustNotInfluenceHamiltonianPath()
    {
        Sequence s1 = new Sequence("catagtc");
        Sequence s2 = new Sequence("taactat");
        Sequence s3 = new Sequence("agactatcc");

        List<Sequence> a = Arrays.asList(s1, s2, s3);
        List<Sequence> b = Arrays.asList(s1, s3, s2);
        List<Sequence> c = Arrays.asList(s3, s2, s1);
        List<Sequence> d = Arrays.asList(s2, s1, s3);
        List<Sequence> e = Arrays.asList(s3, s1, s2);
        List<Sequence> f = Arrays.asList(s2, s3, s1);

        Set<Arc> filteredA = Greedy.filterArcs(Greedy.generateArcs(a, 1, -1, -2), a);
        Set<Arc> filteredB = Greedy.filterArcs(Greedy.generateArcs(b, 1, -1, -2), b);
        Set<Arc> filteredC = Greedy.filterArcs(Greedy.generateArcs(c, 1, -1, -2), c);
        Set<Arc> filteredD = Greedy.filterArcs(Greedy.generateArcs(d, 1, -1, -2), d);
        Set<Arc> filteredE = Greedy.filterArcs(Greedy.generateArcs(e, 1, -1, -2), e);
        Set<Arc> filteredF = Greedy.filterArcs(Greedy.generateArcs(f, 1, -1, -2), f);

        List<SequenceAlignment> pathA = Greedy.hamiltonianPath(filteredA);
        List<SequenceAlignment> pathB = Greedy.hamiltonianPath(filteredB);
        List<SequenceAlignment> pathC = Greedy.hamiltonianPath(filteredC);
        List<SequenceAlignment> pathD = Greedy.hamiltonianPath(filteredD);
        List<SequenceAlignment> pathE = Greedy.hamiltonianPath(filteredE);
        List<SequenceAlignment> pathF = Greedy.hamiltonianPath(filteredF);

        assertEquals(pathA, pathB);
        assertEquals(pathB, pathC);
        assertEquals(pathC, pathD);
        assertEquals(pathD, pathE);
        assertEquals(pathE, pathF);
    }

    @Test
    public void testArcs2Path()
    {
        Sequence s1 = new Sequence("actttacg");
        Sequence s2 = new Sequence("ttgcacgat");
        Sequence s3 = new Sequence("ttgcg");
        Sequence s4 = new Sequence("ggaatctgcgagtta");
        Sequence s5 = new Sequence("tact");

        Arc a1 = new Arc(s1, false, s2, false, 42, false, 1, -1, -2);
        Arc a2 = new Arc(s2, false, s3, false, 42, false, 1, -1, -2);
        Arc a3 = new Arc(s3, false, s4, false, 42, false, 1, -1, -2);
        Arc a4 = new Arc(s4, false, s5, false, 42, false, 1, -1, -2);

        Set<Arc> arcs = new HashSet<Arc>(Arrays.asList(a2, a4, a3, a1));

        List<Arc> res = Greedy.arcs2Path(arcs);

        assertEquals(4, res.size());

        assertEquals(a1, res.get(0));
        assertEquals(a2, res.get(1));
        assertEquals(a3, res.get(2));
        assertEquals(a4, res.get(3));
    }
}
