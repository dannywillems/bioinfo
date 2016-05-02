package be.ac.umons.bioinfo.sequence;

import org.junit.Ignore;
import org.junit.Test;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

public class SequencePairSameTest
{
    @Test
    public void findGapsSimpleTest()
    {
        Sequence initial = new Sequence("cgatctg");
        Sequence s = new Sequence("cgat-c-tg");
        Sequence t = new Sequence("c-gat--ctg");
        SequencePairSame a = new SequencePairSame(initial, s, t);

        int gaps_s[] = {1, 5};
        int gaps_t[] = {8};
        int[][] gaps = a.findGaps();
        assertTrue(Arrays.equals(gaps[0], gaps_s));
        assertTrue(Arrays.equals(gaps[1], gaps_t));
    }

    @Test
    public void rebuildWithGapsSimpleTest()
    {
        Sequence initial = new Sequence("cgatctg");
        Sequence s = new Sequence("cgat-c-tg");
        Sequence t = new Sequence("c-gat--ctg");
        SequencePairSame a = new SequencePairSame(initial, s, t);

        a.rebuildWithGaps();
        assertEquals(new Sequence("c-gat--c-tg"), a.s);
        assertEquals(new Sequence("c-gat--c-tg"), a.t);
    }

    @Test
    public void rebuildWithGapsOrderDoesntChangeTest()
    {
        Sequence initial = new Sequence("cgatctg");
        Sequence s = new Sequence("cgat-c-tg");
        Sequence t = new Sequence("c-gat--ctg");
        SequencePairSame a = new SequencePairSame(initial, t, s);

        a.rebuildWithGaps();
        assertEquals(new Sequence("c-gat--c-tg"), a.s);
        assertEquals(new Sequence("c-gat--c-tg"), a.t);
    }
}
