package be.ac.umons.bioinfo.sequence;

import java.util.*;

import org.junit.Test;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;

public class SequenceAbstractTest
{
    @Test
    public void arrayGapsNoGapsTest()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("actg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {0, 0, 0, 0};
        assertArrayEquals(gaps, s.nb_gaps);
    }

    @Test
    public void arrayGapsGapsTest1()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a-c-t-g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {1, 1, 1, 0};
        assertArrayEquals(gaps, s.nb_gaps);
    }

    @Test
    public void arrayGapsGapsTest2()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a---c--t-g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {3, 2, 1, 0};
        assertArrayEquals(gaps, s.nb_gaps);
    }

    @Test
    public void arrayGapsGapsTest3()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a---c--t-g--");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {3, 2, 1, 2};
        assertArrayEquals(gaps, s.nb_gaps);
    }

    @Test
    public void toStringNoGapsTest()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("actg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void toStringGapsTest1()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a--c--tg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void toStringGapsTest2()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a--c--tg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void toStringGapsTest3()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a--c--tg--");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void toStringGapsTest4()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a--ctg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void toStringGapsTestEnd()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("actg--");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void addGapsNoTest()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("actg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(0, 0);
        assertEquals(initial.toString(), s.toString());

        s.addGaps(0, 1);
        assertEquals(initial.toString(), s.toString());

        s.addGaps(0, 2);
        assertEquals(initial.toString(), s.toString());

        s.addGaps(0, 3);
        assertEquals(initial.toString(), s.toString());
    }

    @Test
    public void addGapsTest1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 1);
        assertEquals("a---ttg", s.toString());
    }

    @Test
    public void addGapsTestTwoTimes1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 1);
        s.addGaps(3, 4);
        assertEquals("a------ttg", s.toString());
    }

    @Test
    public void addGapsTestTwoTimes2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 1);
        s.addGaps(3, 5);
        assertEquals("a---t---tg", s.toString());
    }

    @Test
    public void addGapsEndTest1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 4);
        assertEquals("attg---", s.toString());
    }

    @Test
    public void addGapsEndTest2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 7);
        assertEquals("attg---", s.toString());
    }

    @Test
    public void addGapsBeginAndEndTest2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 7);
        s.addGaps(3, 0);
        assertEquals("---attg---", s.toString());
    }

    @Test
    public void addGapsBeginAndEndTest3()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 7);
        s.addGaps(3, 1);
        assertEquals("a---ttg---", s.toString());
    }

    @Test
    public void addGapsWhenInitialGaps1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("a-tt--g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 1);
        assertEquals("a----tt--g", s.toString());
    }

    @Test
    public void addGapsWhenInitialGaps2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("a-tt--g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGaps(3, 3);
        assertEquals("a-t---t--g", s.toString());
    }
}
