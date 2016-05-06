package be.ac.umons.bioinfo.sequence;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;

public class ConsensusAbstractTest
{
    @Test
    public void updateOffsetSimpleTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("at-cg"), new SequenceAbstract("-tcg-")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("t-cg"), new SequenceAbstract("acct")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();

        int result[] = {0, 1, 1, 1};

        for(int i = 0;i < result.length;i++)
        {
            SequenceAlignmentAbstract s = c.getHamiltonianPath().get(i / 2);
            assertEquals(result[i], s.s.getOffset());
            assertEquals(result[++i], s.t.getOffset());
        }
    }

    @Test
    public void updateOffsetHarderTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtgcaa----"), new SequenceAbstract("ggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatctg-cgagtta"), new SequenceAbstract("---atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        int result[][] = { {0, 1, 7, 3, 0}, {1, 7, 3, 0, 3} };

        for(int j = 0;j < result[0].length;j++)
        {
            SequenceAlignmentAbstract sa = c.getHamiltonianPath().get(j);
            assertEquals(result[0][j], sa.s.getOffset());
            assertEquals(result[1][j], sa.t.getOffset());
        }
    }

    @Test
    public void updateOffsetNoOverlapTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("--------cgcaa"), new SequenceAbstract("cgtaaagt-----")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        int result[][] = { {8}, {0} };

        for(int j = 0;j < result[0].length;j++)
        {
            SequenceAlignmentAbstract sa = c.getHamiltonianPath().get(j);
            assertEquals(result[0][j], sa.s.getOffset());
            assertEquals(result[1][j], sa.t.getOffset());
        }
    }

    @Test
    public void updateOffsetWithMinChangedTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("-----atcgtgcaa----"), new SequenceAbstract("acggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("acggaatctg-cgagtta"), new SequenceAbstract("-----atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();

        int result[][] = { {2, 3, 9, 5, 0}, {3, 9, 5, 0, 5} };

        for(int j = 0;j < result[0].length;j++)
        {
            SequenceAlignmentAbstract sa = c.getHamiltonianPath().get(j);
            assertEquals(result[0][j], sa.s.getOffset());
            assertEquals(result[1][j], sa.t.getOffset());
        }
    }

    @Test
    public void updateOffsetSimpleIdempotentTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("at-cg"), new SequenceAbstract("-tcg-")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("t-cg"), new SequenceAbstract("acct")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        c.updateOffset();

        int result[] = {0, 1, 1, 1};

        for(int i = 0;i < result.length;i++)
        {
            SequenceAlignmentAbstract s = c.getHamiltonianPath().get(i / 2);
            assertEquals(result[i], s.s.getOffset());
            assertEquals(result[++i], s.t.getOffset());
        }
    }

    @Test
    public void updateOffsetHarderIdempotentTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtgcaa----"), new SequenceAbstract("ggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatctg-cgagtta"), new SequenceAbstract("---atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        c.updateOffset();

        int result[][] = { {0, 1, 7, 3, 0}, {1, 7, 3, 0, 3} };

        for(int j = 0;j < result[0].length;j++)
        {
            SequenceAlignmentAbstract sa = c.getHamiltonianPath().get(j);
            assertEquals(result[0][j], sa.s.getOffset());
            assertEquals(result[1][j], sa.t.getOffset());
        }
    }

    @Test
    public void updateOffsetNoOverlapIdempotentTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("--------cgcaa"), new SequenceAbstract("cgtaaagt-----")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        c.updateOffset();

        int result[][] = { {8}, {0} };

        for(int j = 0;j < result[0].length;j++)
        {
            SequenceAlignmentAbstract sa = c.getHamiltonianPath().get(j);
            assertEquals(result[0][j], sa.s.getOffset());
            assertEquals(result[1][j], sa.t.getOffset());
        }
    }

    @Test
    public void updateOffsetMinChangedIdempotentTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("-----atcgtgcaa----"), new SequenceAbstract("acggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("acggaatctg-cgagtta"), new SequenceAbstract("-----atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        c.updateOffset();

        int result[][] = { {2, 3, 9, 5, 0}, {3, 9, 5, 0, 5} };

        for(int j = 0;j < result[0].length;j++)
        {
            SequenceAlignmentAbstract sa = c.getHamiltonianPath().get(j);
            assertEquals(result[0][j], sa.s.getOffset());
            assertEquals(result[1][j], sa.t.getOffset());
        }
    }


    @Test
    public void updateOffsetLongerOneNoOverlapTest()
    {
        List<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("tact"));
        l.add(new Sequence("actttacg"));
        l.add(new Sequence("cgcaa"));
        l.add(new Sequence("atcgtgcaa"));
        l.add(new Sequence("ggaatctgcgagtta"));
        l.add(new Sequence("atcggtc"));

        Greedy g = new Greedy();
        ConsensusAbstract c = new ConsensusAbstract(g.greedy(l, 1, -1, -2));
        c.updateOffset();
        int result[][] = { {16, 8, 8, 11, 3}, {8, 8, 11, 3, 0} };

        for(int j = 0;j < result[0].length;j++)
        {
            assertEquals(result[0][j], c.getHamiltonianPath().get(j).s.getOffset());
            assertEquals(result[1][j], c.getHamiltonianPath().get(j).t.getOffset());
        }
    }

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST sameHamiltonianPath */
    @Test
    public void sameHamiltonianPathSimpleTrueTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact"), new SequenceAbstract("actt")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actt"), new SequenceAbstract("gcaa")));

        ArrayList<SequenceAlignmentAbstract> l2 = new ArrayList<SequenceAlignmentAbstract>();
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact"), new SequenceAbstract("actt")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actt"), new SequenceAbstract("gcaa")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        ConsensusAbstract c2 = new ConsensusAbstract(l2);

        assertTrue(c.sameHamiltonianPath(c2));
    }

    @Test
    public void sameHamiltonianPathSimpleFalseTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact"), new SequenceAbstract("actt")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("gctt"), new SequenceAbstract("gcaa")));

        ArrayList<SequenceAlignmentAbstract> l2 = new ArrayList<SequenceAlignmentAbstract>();
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact"), new SequenceAbstract("actt")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actt"), new SequenceAbstract("gcaa")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        ConsensusAbstract c2 = new ConsensusAbstract(l2);

        assertFalse(c.sameHamiltonianPath(c2));
    }

    @Test
    public void sameHamiltonianPathNotSameSizeTest()
    {
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact"), new SequenceAbstract("actt")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("gctt"), new SequenceAbstract("gcaa")));

        ArrayList<SequenceAlignmentAbstract> l2 = new ArrayList<SequenceAlignmentAbstract>();
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact"), new SequenceAbstract("actt")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actt"), new SequenceAbstract("gcaa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actt"), new SequenceAbstract("gcaa")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        ConsensusAbstract c2 = new ConsensusAbstract(l2);

        assertFalse(c.sameHamiltonianPath(c2));
    }
    /* END TEST sameHamiltonianPath */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST propageGapsDownFrom */
    @Test
    public void propageGapsDownFromSimpleTest()
    {
        /* ------------------------------------------------------------------ */
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtgcaa----"), new SequenceAbstract("ggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatctg-cgagtta"), new SequenceAbstract("---atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        c.propageGapsDownFrom(3);
        c.addEndGaps(); // Needed because the last lines has more a nucleotide
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* WAITED OUTPUT WITHOUT OFFSET GAPS AND END */
        ArrayList<SequenceAlignmentAbstract> l2 = new ArrayList<SequenceAlignmentAbstract>();
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtgcaa----"), new SequenceAbstract("ggaatc-tgcgagtta")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatc-tg-cgagtta"), new SequenceAbstract("---atc-ggtc------")));

        ConsensusAbstract c2 = new ConsensusAbstract(l2);
        c2.updateOffset();
        c2.addEndGaps();
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
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtgcaa----"), new SequenceAbstract("ggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatctg-cgagtta"), new SequenceAbstract("---atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();
        c.propageGapsUpFrom(4);
        c.addEndGaps(); // Needed because the last lines has more a nucleotide
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* WAITED OUTPUT WITHOUT OFFSET GAPS AND END */
        ArrayList<SequenceAlignmentAbstract> l2 = new ArrayList<SequenceAlignmentAbstract>();
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact------"), new SequenceAbstract("-actttacg")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cg-caa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cg-caa"), new SequenceAbstract("atcgtg-caa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtg-caa----"), new SequenceAbstract("ggaatc-tg-cgagtta")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatctg-cgagtta"), new SequenceAbstract("---atcggtc------")));

        ConsensusAbstract c2 = new ConsensusAbstract(l2);
        c2.updateOffset();
        c2.addEndGaps();
        /* ------------------------------------------------------------------ */

        assertTrue(c.sameHamiltonianPath(c2));
    }

    @Test
    public void propageGapsDownNextUpFromSimpleTest()
    {
        /* ------------------------------------------------------------------ */
        ArrayList<SequenceAlignmentAbstract> l = new ArrayList<SequenceAlignmentAbstract>();
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact-----"), new SequenceAbstract("-actttacg")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg---"), new SequenceAbstract("------cgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cgcaa"), new SequenceAbstract("atcgtgcaa")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtgcaa----"), new SequenceAbstract("ggaatc-tgcgagtta")));
        l.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatctg-cgagtta"), new SequenceAbstract("---atcggtc------")));

        ConsensusAbstract c = new ConsensusAbstract(l);
        c.updateOffset();

        // In real cases, the down will be done first because the position (6)
        // is less than 8. We also need to add 1 to 8 because the first
        // propagation move to the left the 4th line. It's what it must happen.
        c.propageGapsDownFrom(3);
        c.propageGapsUpFrom(4);
        c.addEndGaps(); // Needed because the last lines has more a nucleotide
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* WAITED OUTPUT WITHOUT OFFSET GAPS AND END */
        ArrayList<SequenceAlignmentAbstract> l2 = new ArrayList<SequenceAlignmentAbstract>();
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("tact------"), new SequenceAbstract("-actttacg-")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("actttacg----"), new SequenceAbstract("------cg-caa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("----cg-caa"), new SequenceAbstract("atcgtg-caa")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("---atcgtg-caa----"), new SequenceAbstract("ggaatc-tg-cgagtta")));
        l2.add(new SequenceAlignmentAbstract(new SequenceAbstract("ggaatc-tg-cgagtta"), new SequenceAbstract("---atc-ggtc------")));

        ConsensusAbstract c2 = new ConsensusAbstract(l2);
        c2.updateOffset();
        c2.addEndGaps();
        /* ------------------------------------------------------------------ */

        assertTrue(c.sameHamiltonianPath(c2));
    }

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

        ConsensusAbstract c = new ConsensusAbstract(result);
        c.computeAlignment();
        ArrayList<SequenceAbstract> align = c.getAlignment();

        ArrayList<SequenceAbstract> l2 = new ArrayList<SequenceAbstract>();
        l2.add(new SequenceAbstract("----------------cgcaa"));
        l2.add(new SequenceAbstract("--------cgtaaagt-----"));
        l2.add(new SequenceAbstract("--------agta---------"));
        l2.add(new SequenceAbstract("-----------atcggtc---"));
        l2.add(new SequenceAbstract("---atcgtgcaa---------"));
        l2.add(new SequenceAbstract("ggaatc-tgcgagtta-----"));

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

        ConsensusAbstract c = new ConsensusAbstract(result);
        c.computeAlignment();
        ArrayList<SequenceAbstract> align = c.getAlignment();

        ArrayList<SequenceAbstract> l2 = new ArrayList<SequenceAbstract>();
        l2.add(new SequenceAbstract("tact-------------"));
        l2.add(new SequenceAbstract("-actttacg--------"));
        l2.add(new SequenceAbstract("-------cg-caa----"));
        l2.add(new SequenceAbstract("---atcgtg-caa----"));
        l2.add(new SequenceAbstract("ggaatc-tg-cgagtta"));
        l2.add(new SequenceAbstract("---atc-ggtc------"));

        for(int l = 0;l < l2.size();l++)
            assertEquals(l2.get(l), align.get(l));
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST getBase */
    @Test
    public void getBaseTest()
    {
        ConsensusAbstract c = new ConsensusAbstract();

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
        ConsensusAbstract c = new ConsensusAbstract();

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
        ConsensusAbstract c = new ConsensusAbstract();

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
        ConsensusAbstract c = new ConsensusAbstract();

        HashMap<Character, Integer> s = new HashMap<Character, Integer>();
        s.put(new Character('c'), 5);
        s.put(new Character('g'), 2);
        s.put(new Character('t'), 4);
        s.put(new Character('-'), 6);
        assertEquals(c.getBase(s, false), 'c');
    }
    /* END TEST getBase */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST build with remove_if_max_gap false */
    @Test
    public void buildNoRemoveSimpleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("acttt"));
        l.add(new SequenceAbstract("tcttg"));
        l.add(new SequenceAbstract("atgag"));

        ConsensusAbstract c = new ConsensusAbstract();
        c.setAlignment(l);

        assertEquals(c.build(false).toString(), "acttg");
    }

    @Test
    public void buildNoGapNoRemoveMultipleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("acttt"));
        l.add(new SequenceAbstract("tattg"));
        l.add(new SequenceAbstract("atgag"));

        ConsensusAbstract c = new ConsensusAbstract();
        c.setAlignment(l);

        assertEquals(c.build(false).toString(), "aattg");
    }

    @Test
    public void buildGapNoRemoveMultipleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("a-ttt"));
        l.add(new SequenceAbstract("t-ttg"));
        l.add(new SequenceAbstract("atgag"));

        ConsensusAbstract c = new ConsensusAbstract();
        c.setAlignment(l);

        assertEquals(c.build(false).toString(), "atttg");
    }

    @Test
    public void buildGapSimpleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("a-ttt"));
        l.add(new SequenceAbstract("tcttg"));
        l.add(new SequenceAbstract("acgag"));

        ConsensusAbstract c = new ConsensusAbstract(); c.setAlignment(l);

        assertEquals(c.build(false).toString(), "acttg");
    }
    /* END TEST build with remove_if_max_gap false */
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /* BEGIN TEST build with remove_if_max_gap true */
    @Test
    public void buildGapRemoveMultipleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("a-ttt"));
        l.add(new SequenceAbstract("t-ttg"));
        l.add(new SequenceAbstract("atgag"));

        ConsensusAbstract c = new ConsensusAbstract();
        c.setAlignment(l);

        assertEquals(c.build(true).toString(), "attg");
    }

    @Test
    public void buildGapRemoveSimpleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("a-ttt"));
        l.add(new SequenceAbstract("tcttg"));
        l.add(new SequenceAbstract("acgag"));

        ConsensusAbstract c = new ConsensusAbstract(); c.setAlignment(l);

        assertEquals(c.build(true).toString(), "acttg");
    }

    @Test
    public void buildRemoveSimpleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("acttt"));
        l.add(new SequenceAbstract("tcttg"));
        l.add(new SequenceAbstract("atgag"));

        ConsensusAbstract c = new ConsensusAbstract();
        c.setAlignment(l);

        assertEquals(c.build(true).toString(), "acttg");
    }

    @Test
    public void buildNoGapRemoveMultipleTest()
    {
        ArrayList<SequenceAbstract> l = new ArrayList<SequenceAbstract>();
        l.add(new SequenceAbstract("acttt"));
        l.add(new SequenceAbstract("tattg"));
        l.add(new SequenceAbstract("atgag"));

        ConsensusAbstract c = new ConsensusAbstract();
        c.setAlignment(l);

        assertEquals(c.build(true).toString(), "aattg");
    }

    @Test
    public void buildRealTestNoRemove()
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

        ConsensusAbstract c = new ConsensusAbstract(result);
        c.computeAlignment();
        SequenceAbstract s_final = c.build(false);

        SequenceAbstract s = new SequenceAbstract("gacatcacgtcaagtta");
        assertEquals(s, s_final);
    }

    @Test
    public void buildRealTestRemove()
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

        ConsensusAbstract c = new ConsensusAbstract(result);
        c.computeAlignment();
        SequenceAbstract s_final = c.build(true);

        SequenceAbstract s = new SequenceAbstract("atccgc");
        assertEquals(s, s_final);
    }
    /* END TEST build with remove_if_max_gap true */
    /* ---------------------------------------------------------------------- */
}
