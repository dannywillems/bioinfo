package be.ac.umons.bioinfo.sequence;

import static junit.framework.TestCase.assertEquals;

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

}
