package be.ac.umons.bioinfo.consensus;

import be.ac.umons.bioinfo.sequence.*;
import be.ac.umons.bioinfo.greedy.*;

import org.junit.Ignore;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;

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

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST computeOffset */
    @Test
    public void computeOffsetSimpleTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("at-cg"), new Sequence("-tcg-"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("t-cg"), new Sequence("acct"), 0, 0));

        Consensus c = new Consensus(l);
        c.computeOffset();
        int gaps[][] = c.getOffset();
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
        c.computeOffset();
        int gaps[][] = c.getOffset();
        int result[][] = { {0, 1, 7, 3, 0}, {1, 7, 3, 0, 3} };

        for(int j = 0;j < result[0].length;j++)
        {
            assertEquals(result[0][j], gaps[0][j]);
            assertEquals(result[1][j], gaps[1][j]);
        }
    }

    @Test
    public void computeOffsetNoOverlapTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("--------cgcaa"), new Sequence("cgtaaagt-----"), 0, 0));

        Consensus c = new Consensus(l);
        c.computeOffset();
        int gaps[][] = c.getOffset();
        int result[][] = { {8}, {0} };

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
        c.computeOffset();
        int gaps[][] = c.getOffset();
        int result[][] = { {2, 3, 9, 5, 0}, {3, 9, 5, 0, 5} };

        for(int j = 0;j < result[0].length;j++)
        {
            assertEquals(result[0][j], gaps[0][j]);
            assertEquals(result[1][j], gaps[1][j]);
        }
    }

    @Test
    public void computeOffsetLongerOneNoOverlapTest()
    {
        List<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("tact"));
        l.add(new Sequence("actttacg"));
        l.add(new Sequence("cgcaa"));
        l.add(new Sequence("atcgtgcaa"));
        l.add(new Sequence("ggaatctgcgagtta"));
        l.add(new Sequence("atcggtc"));

        Greedy g = new Greedy();
        Consensus c = new Consensus(g.greedy(l, 1, -1, -2));
        c.computeOffset();
        int gaps[][] = c.getOffset();
        int result[][] = { {16, 8, 8, 11, 3}, {8, 8, 11, 3, 0} };

        for(int j = 0;j < result[0].length;j++)
        {
            assertEquals(result[0][j], gaps[0][j]);
            assertEquals(result[1][j], gaps[1][j]);
        }
    }
    /* END TEST computeOffset */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST showOffset */
    @Test
    public void showOffsetNoOverlapTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("--------cgcaa"), new Sequence("cgtaaagt-----"), 0, 0));

        Consensus c = new Consensus(l);
        c.addOffset();
        Consensus c2 = new Consensus(l);

        assertTrue(c.sameHamiltonianPath(c2));
    }

    @Test
    public void showOffsetLongerOneNoOverlapTest()
    {
        List<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("tact"));
        l.add(new Sequence("actttacg"));
        l.add(new Sequence("cgcaa"));
        l.add(new Sequence("atcgtgcaa"));
        l.add(new Sequence("ggaatctgcgagtta"));
        l.add(new Sequence("atcggtc"));

        Greedy g = new Greedy();
        Consensus c = new Consensus(g.greedy(l, 1, -1, -2));
        c.addOffset();
        c.fillEndWithGaps();

        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("----------------cgcaa"), new Sequence("--------cgtaaagt-----"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("--------cgtaaagt-----"), new Sequence("--------agta---------"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("--------agta---------"), new Sequence("-----------atcggtc---"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("-----------atcggtc---"), new Sequence("---atcgtgcaa---------"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("---atcgtgcaa---------"), new Sequence("ggaatc-tgcgagtta-----"), 0, 0));

        Consensus c2 = new Consensus(l2);
        assertTrue(c.sameHamiltonianPath(c2));
    }

    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST propageGapsDownFrom */
    @Test
    public void propageGapsDownFromSimpleTest()
    {
        /* ------------------------------------------------------------------ */
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact-----"), new Sequence("-actttacg"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("actttacg---"), new Sequence("------cgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("----cgcaa"), new Sequence("atcgtgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("---atcgtgcaa----"), new Sequence("ggaatc-tgcgagtta"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("ggaatctg-cgagtta"), new Sequence("---atcggtc------"), 0, 0));

        Consensus c = new Consensus(l);
        c.addOffset();
        c.fillEndWithGaps();
        c.propageGapsDownFrom(3, 6);
        c.fillEndWithGaps(); // Needed because the last lines has more a nucleotide
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* WAITED OUTPUT WITHOUT OFFSET GAPS AND END */
        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("tact-----"), new Sequence("-actttacg"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actttacg---"), new Sequence("------cgcaa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("----cgcaa"), new Sequence("atcgtgcaa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("---atcgtgcaa----"), new Sequence("ggaatc-tgcgagtta"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("ggaatc-tg-cgagtta"), new Sequence("---atc-ggtc------"), 0, 0));

        Consensus c2 = new Consensus(l2);
        c2.addOffset();
        c2.fillEndWithGaps();
        /* ------------------------------------------------------------------ */

        assertTrue(c.sameHamiltonianPath(c2));
    }
    /* END TEST propageGapsDownFrom */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST propageGapsUpFrom */
    @Test
    public void propageGapsUpFromSimpleTest()
    {
        /* ------------------------------------------------------------------ */
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact-----"), new Sequence("-actttacg"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("actttacg---"), new Sequence("------cgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("----cgcaa"), new Sequence("atcgtgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("---atcgtgcaa----"), new Sequence("ggaatc-tgcgagtta"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("ggaatctg-cgagtta"), new Sequence("---atcggtc------"), 0, 0));

        Consensus c = new Consensus(l);
        c.addOffset();
        c.fillEndWithGaps();
        c.propageGapsUpFrom(4, 8);
        c.fillEndWithGaps(); // Needed because the last lines has more a nucleotide
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* WAITED OUTPUT WITHOUT OFFSET GAPS AND END */
        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("tact------"), new Sequence("-actttac-g"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actttac-g---"), new Sequence("------c-gcaa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("----c-gcaa"), new Sequence("atcgt-gcaa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("---atcgt-gcaa----"), new Sequence("ggaatc-t-gcgagtta"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("ggaatctg-cgagtta"), new Sequence("---atcggtc------"), 0, 0));

        Consensus c2 = new Consensus(l2);
        c2.addOffset();
        c2.fillEndWithGaps();
        /* ------------------------------------------------------------------ */

        assertTrue(c.sameHamiltonianPath(c2));
    }

    @Test
    public void propageGapsDownNextUpFromSimpleTest()
    {
        /* ------------------------------------------------------------------ */
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact-----"), new Sequence("-actttacg"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("actttacg---"), new Sequence("------cgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("----cgcaa"), new Sequence("atcgtgcaa"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("---atcgtgcaa----"), new Sequence("ggaatc-tgcgagtta"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("ggaatctg-cgagtta"), new Sequence("---atcggtc------"), 0, 0));

        Consensus c = new Consensus(l);
        c.addOffset();
        c.fillEndWithGaps();

        // In real cases, the down will be done first because the position (6)
        // is less than 8. We also need to add 1 to 8 because the first
        // propagation move to the left the 4th line. It's what it must happen.
        c.propageGapsDownFrom(3, 6);
        c.propageGapsUpFrom(4, 9);
        c.fillEndWithGaps(); // Needed because the last lines has more a nucleotide
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* WAITED OUTPUT WITHOUT OFFSET GAPS AND END */
        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("tact------"), new Sequence("-actttacg-"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actttacg----"), new Sequence("------cg-caa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("----cg-caa"), new Sequence("atcgtg-caa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("---atcgtg-caa----"), new Sequence("ggaatc-tg-cgagtta"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("ggaatc-tg-cgagtta"), new Sequence("---atc-ggtc------"), 0, 0));

        Consensus c2 = new Consensus(l2);
        c2.addOffset();
        c2.fillEndWithGaps();
        /* ------------------------------------------------------------------ */

        assertTrue(c.sameHamiltonianPath(c2));
    }

    /* END TEST propageGapsUpFrom */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST sameHamiltonianPath */
    @Test
    public void sameHamiltonianPathSimpleTrueTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact"), new Sequence("actt"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("actt"), new Sequence("gcaa"), 0, 0));

        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("tact"), new Sequence("actt"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actt"), new Sequence("gcaa"), 0, 0));

        Consensus c = new Consensus(l);
        Consensus c2 = new Consensus(l2);

        assertTrue(c.sameHamiltonianPath(c2));
    }

    @Test
    public void sameHamiltonianPathSimpleFalseTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact"), new Sequence("actt"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("gctt"), new Sequence("gcaa"), 0, 0));

        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("tact"), new Sequence("actt"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actt"), new Sequence("gcaa"), 0, 0));

        Consensus c = new Consensus(l);
        Consensus c2 = new Consensus(l2);

        assertFalse(c.sameHamiltonianPath(c2));
    }

    @Test
    public void sameHamiltonianPathNotSameSizeTest()
    {
        ArrayList<SequenceAlignment> l = new ArrayList<SequenceAlignment>();
        l.add(new SequenceAlignment(new Sequence("tact"), new Sequence("actt"), 0, 0));
        l.add(new SequenceAlignment(new Sequence("gctt"), new Sequence("gcaa"), 0, 0));

        ArrayList<SequenceAlignment> l2 = new ArrayList<SequenceAlignment>();
        l2.add(new SequenceAlignment(new Sequence("tact"), new Sequence("actt"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actt"), new Sequence("gcaa"), 0, 0));
        l2.add(new SequenceAlignment(new Sequence("actt"), new Sequence("gcaa"), 0, 0));

        Consensus c = new Consensus(l);
        Consensus c2 = new Consensus(l2);

        assertFalse(c.sameHamiltonianPath(c2));
    }
    /* END TEST sameHamiltonianPath */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST computeAlignment */
    @Test
    public void computeAlignmentNoPropagationTest()
    {
        List<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("tact"));
        l.add(new Sequence("actttacg"));
        l.add(new Sequence("cgcaa"));
        l.add(new Sequence("atcgtgcaa"));
        l.add(new Sequence("ggaatctgcgagtta"));
        l.add(new Sequence("atcggtc"));

        Greedy g = new Greedy();
        List<SequenceAlignment> result = g.greedy(l, 1, -1, -2);

        Consensus c = new Consensus(result);
        c.computeAlignment();
        ArrayList<Sequence> align = c.getAlignment();

        ArrayList<Sequence> l2 = new ArrayList<Sequence>();
        l2.add(new Sequence("----------------cgcaa"));
        l2.add(new Sequence("--------cgtaaagt-----"));
        l2.add(new Sequence("--------agta---------"));
        l2.add(new Sequence("-----------atcggtc---"));
        l2.add(new Sequence("---atcgtgcaa---------"));
        l2.add(new Sequence("ggaatc-tgcgagtta-----"));

        for(int i = 0;i < l2.size();i++)
            assertTrue(l2.get(i).equals(align.get(i)));
    }
    /* END TEST computeAlignment */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void computeAlignmentOnePropagationDownAndUpTest()
    {
        Sequence f = new Sequence("actttacg");
        Sequence g = new Sequence("ttgcacgat");
        Sequence h = new Sequence("ttgcg");
        Sequence i = new Sequence("ggaatctgcgagtta");
        Sequence j = new Sequence("tact");
        Sequence k = new Sequence("gaccgat");

        List<Sequence> list = new ArrayList<Sequence>();

        list.add(j);
        list.add(f);
        list.add(k);
        list.add(h);
        list.add(g);
        list.add(i);

        Greedy greed = new Greedy();
        List<SequenceAlignment> result = greed.greedy(list, 1, -1, -2);

        Consensus c = new Consensus(result);
        System.out.println("############################");
        c.computeAlignment();
        ArrayList<Sequence> align = c.getAlignment();

        ArrayList<Sequence> l2 = new ArrayList<Sequence>();
        l2.add(new Sequence("tact-------------"));
        l2.add(new Sequence("-actttacg--------"));
        l2.add(new Sequence("-------cg-caa----"));
        l2.add(new Sequence("---atcgtg-caa----"));
        l2.add(new Sequence("ggaatc-tg-cgagtta"));
        l2.add(new Sequence("---atc-ggtc------"));

        for(int l = 0;l < l2.size();l++)
            assertEquals(l2.get(l), align.get(l));
    }
    /* ---------------------------------------------------------------------- */
}
