package be.ac.umons.sequence.Test;

/**
 * Created by aline on 31/03/16.
 */
import be.ac.umons.sequence.Sequence;
import be.ac.umons.sequence.SequenceAlignment;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static groovy.util.GroovyTestCase.assertEquals;
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

        SequenceAlignment result = s.semiGlobalAlignment(t, 1, -1, -2);

        assertEquals(result.s1.getSize(), result.s2.getSize());

        assertEquals(new Sequence("cagca-cttggattctcgg"), result.s1);
        assertEquals(new Sequence("---cagcgtgg--------"), result.s2);
    }
    @Test
    public void semiGlobalAlignmentReverseTest()
    {
        Sequence s = new Sequence("cagcacttggattctcgg");
        Sequence t = new Sequence("cagcgtgg");

        SequenceAlignment result = t.semiGlobalAlignment(s, 1, -1, -2);
        assertEquals(result.s1.getSize(), result.s2.getSize());

        assertEquals(new Sequence("cagca-cttggattctcgg"), result.s2);
        assertEquals(new Sequence("---cagcgtgg--------"), result.s1);


    }

    @Test
    public void semiGlobalAlignmentTest2()
    {
        Sequence s = new Sequence("aggagaagaattcaccgctat");
        Sequence t = new Sequence("ttccccttattcaattctaa");

        SequenceAlignment result = s.semiGlobalAlignment(t, 1, -1, -2);

        assertEquals(result.s1.getSize(), result.s2.getSize());

        assertEquals(new Sequence("aggagaagaattcaccgctat----------"), result.s1);
        assertEquals(new Sequence("----------ttcccct-tattcaattctaa"), result.s2);
    }
}