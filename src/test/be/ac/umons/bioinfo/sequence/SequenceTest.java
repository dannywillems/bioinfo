package be.ac.umons.bioinfo.sequence;

/**
 * Created by aline on 31/03/16.
 */
import be.ac.umons.bioinfo.sequence.Sequence;
import be.ac.umons.bioinfo.sequence.SequenceAlignment;
import org.junit.Test;

import java.util.*;

import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;

public class SequenceTest
{
    @Test
    public void complementTest()
    {
        Sequence t = new Sequence("cagcgtgg");
        Sequence res = t.complement();
        Sequence wanted = new Sequence("ccacgctg");

        assertEquals(wanted, res);
    }
    @Test
    public void complementOfComplementMustBeIdentity()
    {
        List<Sequence> sequences = new ArrayList<Sequence>();
        sequences.add(new Sequence("GATTACA"));
        sequences.add(new Sequence("ATCATTAGTGG"));
        sequences.add(new Sequence("A"));
        sequences.add(new Sequence("G"));
        sequences.add(new Sequence("T"));
        sequences.add(new Sequence("C"));

        for(Sequence s : sequences)
            assertEquals(s, s.complement().complement());
    }

    @Test
    public void alignmentScoreTest()
    {
        Sequence s1 = new Sequence("GA-CGGATTAG");
        Sequence s2 = new Sequence("GATCGGAATAG");

        int score = s1.alignmentScore(s2, 1, -1, -2);
        int reverseScore = s2.alignmentScore(s1, 1, -1, -2);

        assertEquals(6, score);
        assertEquals(6, reverseScore);
    }

    @Test
    public void semiGlobalAlignmentTest()
    {
        Sequence s = new Sequence("cagcacttggattctcgg");
        Sequence t = new Sequence("cagcgtgg");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(s,t, 1, -1, -2);


        SequenceAlignment first = result.get(0); //(s,t)
        SequenceAlignment second = result.get(1);//(t,s)

        //System.out.println(first.start);
        //System.out.println(first.end);

        //System.out.println(second.start);
        //System.out.println(second.end);

        assertEquals(new Sequence("cagca-cttggattctcgg"), second.s2);
        assertEquals(new Sequence("---cagcgtgg--------"), second.s1);
    }

    @Test
    public void semiGlobalAlignmentReverseTest()
    {
        Sequence s = new Sequence("cagcacttggattctcgg");
        Sequence t = new Sequence("cagcgtgg");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(t, s, 1, -1, -2);

        SequenceAlignment first = result.get(0); //(s,t)
        SequenceAlignment second = result.get(1);//(t,s)


        assertEquals(new Sequence("cagca-cttggattctcgg"), first.s2);
        assertEquals(new Sequence("---cagcgtgg--------"), first.s1);


    }

    @Test
    public void semiGlobalAlignmentTest2()
    {
        Sequence s = new Sequence("aggagaagaattcaccgctat");
        Sequence t = new Sequence("ttccccttattcaattctaa");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(s, t, 1, -1, -2);


        SequenceAlignment first = result.get(0); //(s,t)
        SequenceAlignment second = result.get(1);//(t,s)


        assertEquals(first.s1.getSize(), first.s2.getSize());

        assertEquals(new Sequence("aggagaagaattcaccgctat----------"), first.s1);
        assertEquals(new Sequence("----------ttc-cccttattcaattctaa"), first.s2);
        assertEquals(11, first.score);

        assertEquals(new Sequence("ttccccttattcaattctaa--------------------"), second.s1);
        assertEquals(new Sequence("-------------------aggagaagaattcaccgctat"), second.s2);
        assertEquals(1, second.score);
    }


    @Test

    public void semiGlobalAlignmentTest3()
    {
        Sequence s = new Sequence("attagaccatgcggc");
        Sequence t = new Sequence("atcggcattcagt");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(s, t, 1, -1, -2);


        SequenceAlignment first = result.get(0); //(s,t)
        SequenceAlignment second = result.get(1);//(t,s)

        assertEquals(first.s1.getSize(),first.s2.getSize());

        assertEquals(new Sequence("attagaccatgcggc-------"), first.s1);
        assertEquals(new Sequence("--------at-cggcattcagt"), first.s2);
        assertEquals(7, first.score);

        assertEquals(new Sequence("atcggcattcagt---------"), second.s1);
        assertEquals(new Sequence("------att-agaccatgcggc"), second.s2);
        assertEquals(7, second.score);

    }

    /*
    @Test
    public void containedSemiGlobalAlignmentTest()
    {
        Sequence s = new Sequence("attagaccatgcggc");
        Sequence t = new Sequence("tagacca");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(s,t, 1, -1, -2);

        assertEquals(result.size(),0);

        //assertEquals(new Sequence("attagaccatgcggc"), result.start);
        //assertEquals(new Sequence("--tagacca------"), result.end);


    }
    */

    @Test
    public void reverseAlignementTest()
    {
        Sequence g = new Sequence("taactat");
        Sequence h = new Sequence("agactatcc");

        List<SequenceAlignment> correct = Sequence.semiGlobalAlignment(g,  h, 1, -1, -2);
        List<SequenceAlignment> reverse = Sequence.semiGlobalAlignment(h, g, 1, -1, -2);

        SequenceAlignment correctPrem = correct.get(0);
        SequenceAlignment correctSec = correct.get(1);

        SequenceAlignment reversePrem = reverse.get(0);
        SequenceAlignment reverseSec = reverse.get(1);

        assertEquals(correctPrem.s1, reverseSec.s1);
        assertEquals(correctPrem.s2, reverseSec.s2);
        assertEquals(correctSec.s1, reversePrem.s1);
        assertEquals(correctSec.s2, reversePrem.s2);


    }

    @Test
    public void insertionOrderDoesntChangeTheResultNoOverlap()
    {
        Greedy greedy = new Greedy();

        Sequence f = new Sequence("catagtc");
        Sequence g = new Sequence("taactat");
        Sequence h = new Sequence("agactatcc");

        List<Sequence> fgh = Arrays.asList(f, g, h);
        List<SequenceAlignment> result_fgh = greedy.computePath(fgh, 1, -1, -2);

        List<Sequence> gfh = Arrays.asList(g, f, h);
        List<SequenceAlignment> result_gfh = greedy.computePath(gfh, 1, -1, -2);

        List<Sequence> ghf = Arrays.asList(g, h, f);
        List<SequenceAlignment> result_ghf = greedy.computePath(ghf, 1, -1, -2);

        List<Sequence> fhg = Arrays.asList(f, h, g);
        List<SequenceAlignment> result_fhg = greedy.computePath(fhg, 1, -1, -2);

        // First, we test the path length. We only change the first argument for the test because of equality transivity.
        assertEquals(result_fgh.size(), result_fhg.size()); // fgh VS fhg
        assertEquals(result_gfh.size(), result_fhg.size()); // gfh VS fhg
        assertEquals(result_ghf.size(), result_fhg.size()); // ghf VS fhg

        // Now, check the path
        for (int i = 0; i < result_fgh.size(); i++)
        {
            SequenceAlignment ref = result_fhg.get(i);

            assertEquals(result_fgh.get(i).s1, ref.s1);
            assertEquals(result_fgh.get(i).s2, ref.s2);
            assertEquals(result_gfh.get(i).s1, ref.s1);
            assertEquals(result_gfh.get(i).s2, ref.s2);
        }
    }

    @Test
    public void orderDoesNotAffectSemiGlobalAlignment()
    {
        Sequence a = new Sequence("tact");
        Sequence b = new Sequence("cgtaaagt");
        Sequence c = new Sequence("catagtc");
        Sequence d = new Sequence("taactat");
        Sequence e = new Sequence("agactatcc");

        List<Sequence> sequences = Arrays.asList(a,b,c,d,e);

        for(Sequence s1 : sequences)
        {
            for(Sequence s2 : sequences)
            {
                Set<SequenceAlignment> x = new HashSet<SequenceAlignment>(Sequence.semiGlobalAlignment(s1, s2, 1, -1, -2));
                Set<SequenceAlignment> y = new HashSet<SequenceAlignment>(Sequence.semiGlobalAlignment(s1, s2, 1, -1, -2));

                assertTrue(x.equals(y));
            }
        }
    }
}
