package be.ac.umons.bioinfo.sequence;

import org.junit.Ignore;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

import java.lang.Character;
import java.lang.Integer;

public class ConsensusTest
{
    /* BEGIN TEST getBase */
    @Test
    public void getBaseTest()
    {
        Consensus c = new Consensus(null);

        HashMap<Character, Integer> s = new HashMap<Character, Integer>();
        s.put(new Character('c'), 5);
        s.put(new Character('g'), 2);
        s.put(new Character('t'), 4);
        s.put(new Character('-'), 1);
        assertEquals(c.getBase(s, false), 'c');
    }

    @Test
    public void getBaseMultipleTest()
    {
        Consensus c = new Consensus(null);

        HashMap<Character, Integer> s = new HashMap<Character, Integer>();
        s.put(new Character('c'), 5);
        s.put(new Character('g'), 5);
        s.put(new Character('t'), 4);
        s.put(new Character('-'), 1);
        assertEquals(c.getBase(s, false), 'c');
    }

    @Test
    public void getBaseNoGapMultipleTest()
    {
        Consensus c = new Consensus(null);

        HashMap<Character, Integer> s = new HashMap<Character, Integer>();
        s.put(new Character('c'), 5);
        s.put(new Character('g'), 5);
        s.put(new Character('t'), 4);
        s.put(new Character('-'), 6);
        assertEquals(c.getBase(s, false), 'c');
    }

    @Test
    public void getBaseNoGapTest()
    {
        Consensus c = new Consensus(null);

        HashMap<Character, Integer> s = new HashMap<Character, Integer>();
        s.put(new Character('c'), 5);
        s.put(new Character('g'), 2);
        s.put(new Character('t'), 4);
        s.put(new Character('-'), 6);
        assertEquals(c.getBase(s, false), 'c');
    }
    /* END TEST getBase */

    /* BEGIN TEST build with remove_if_max_gap false */
    @Test
    public void buildNoRemoveSimpleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("acttt"));
        l.add(new Sequence("tcttg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build(false).toString(), "acttg");
    }

    @Test
    public void buildNoGapNoRemoveMultipleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("acttt"));
        l.add(new Sequence("tattg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build(false).toString(), "aattg");
    }

    @Test
    public void buildGapNoRemoveMultipleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("a-ttt"));
        l.add(new Sequence("t-ttg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build(false).toString(), "atttg");
    }

    @Test
    public void buildGapSimpleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("a-ttt"));
        l.add(new Sequence("tcttg"));
        l.add(new Sequence("acgag"));

        Consensus c = new Consensus(null); c.setAlignment(l);

        assertEquals(c.build(false).toString(), "acttg");
    }
    /* END TEST build with remove_if_max_gap false */

    /* BEGIN TEST build with remove_if_max_gap true */
    @Test
    public void buildGapRemoveMultipleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("a-ttt"));
        l.add(new Sequence("t-ttg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build(true).toString(), "attg");
    }

    @Test
    public void buildGapRemoveSimpleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("a-ttt"));
        l.add(new Sequence("tcttg"));
        l.add(new Sequence("acgag"));

        Consensus c = new Consensus(null); c.setAlignment(l);

        assertEquals(c.build(true).toString(), "acttg");
    }

    @Test
    public void buildRemoveSimpleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("acttt"));
        l.add(new Sequence("tcttg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build(true).toString(), "acttg");
    }

    @Test
    public void buildNoGapRemoveMultipleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("acttt"));
        l.add(new Sequence("tattg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build(true).toString(), "aattg");
    }
    /* END TEST build with remove_if_max_gap true */

    @Test
    public void computeOffsetSimpleTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("at-cg"), new Sequence("-tcg-"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("t-cg"), new Sequence("acct"), 0, 0));

        Consensus c = new Consensus(l);
        int gaps[][] = c.computeOffset();
        int result[][] = { {0, 1}, {1, 1} };

        for(int i = 0;i < result.length;i++)
        {
            for(int j = 0;j < result[0].length;j++)
                assertEquals(result[i][j], gaps[i][j]);
        }
    }

    @Test
    public void computeOffsetHarderTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact-----"), new Sequence("-actttacg"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("actttacg---"), new Sequence("------cgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("----cgcaa"), new Sequence("atcgtgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("---atcgtgcaa----"), new Sequence("ggaatc-tgcgagtta"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("ggaatctg-cgagtta"), new Sequence("---atcggtc------"), 0, 0));

        Consensus c = new Consensus(l);
        int gaps[][] = c.computeOffset();
        int result[][] = { {0, 1, 7, 3, 0}, {1, 7, 3, 0, 3} };

        for(int j = 0;j < result[0].length;j++)
        {
            assertEquals(result[0][j], gaps[0][j]);
            assertEquals(result[1][j], gaps[1][j]);
        }
    }

    @Test
    public void computeOffsetWithMinChangedTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact-----"), new Sequence("-actttacg"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("actttacg---"), new Sequence("------cgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("----cgcaa"), new Sequence("atcgtgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("-----atcgtgcaa----"), new Sequence("acggaatc-tgcgagtta"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("acggaatctg-cgagtta"), new Sequence("-----atcggtc------"), 0, 0));

        Consensus c = new Consensus(l);
        int gaps[][] = c.computeOffset();
        int result[][] = { {2, 3, 9, 5, 0}, {3, 9, 5, 0, 5} };

        for(int j = 0;j < result[0].length;j++)
        {
            assertEquals(result[0][j], gaps[0][j]);
            assertEquals(result[1][j], gaps[1][j]);
        }
    }
}
