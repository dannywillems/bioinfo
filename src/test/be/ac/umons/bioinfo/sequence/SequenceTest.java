package be.ac.umons.bioinfo.sequence;

/**
 * Created by aline on 31/03/16.
 */
import be.ac.umons.bioinfo.sequence.Sequence;
import be.ac.umons.bioinfo.sequence.SequenceAlignment;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.TestCase.assertEquals;

public class SequenceTest
{
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
    

        //SequenceAlignment first = result.get(0); //(s,t)
        //SequenceAlignment second = result.get(1);//(t,s)


        assertEquals(result.size(), 0);

        //assertEquals(new Sequence("cagca-cttggattctcgg"), first.s1);
        //assertEquals(new Sequence("---cagcgtgg--------"), first.s2);
    }

    @Test
    public void semiGlobalAlignmentReverseTest()
    {
        Sequence s = new Sequence("cagcacttggattctcgg");
        Sequence t = new Sequence("cagcgtgg");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(t, s, 1, -1, -2);

        //SequenceAlignment first = result.get(0); //(s,t)
        //SequenceAlignment second = result.get(1);//(t,s)

        assertEquals(0, result.size());

        //assertEquals(new Sequence("cagca-cttggattctcgg"), first.s2);
        //assertEquals(new Sequence("---cagcgtgg--------"), first.s1);


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

    }

    /*
    @Test
    public void containedSemiGlobalAlignmentTest()
    {
        Sequence s = new Sequence("attagaccatgcggc");
        Sequence t = new Sequence("tagacca");

        List<SequenceAlignment> result = Sequence.semiGlobalAlignment(s,t, 1, -1, -2);

        assertEquals(result.size(),0);

        //assertEquals(new Sequence("attagaccatgcggc"), result.s1);
        //assertEquals(new Sequence("--tagacca------"), result.s2);


    }
    */
}
