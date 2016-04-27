package be.ac.umons.bioinfo.sequence;

import org.junit.Ignore;
import static junit.framework.TestCase.assertEquals;

import org.junit.Test;
import java.util.HashMap;
import java.util.ArrayList;
import java.lang.Character;
import java.lang.Integer;

public class ConsensusTest
{
    @Test
    public void getBaseTest()
    {
        Consensus c = new Consensus(null);

        HashMap<Character, Integer> s = new HashMap<Character, Integer>();
        s.put(new Character('c'), 5);
        s.put(new Character('g'), 2);
        s.put(new Character('t'), 4);
        s.put(new Character('-'), 1);
        assertEquals(c.getBase(s), 'c');
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
        assertEquals(c.getBase(s), 'c');
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
        assertEquals(c.getBase(s), 'c');
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
        assertEquals(c.getBase(s), 'c');
    }

    @Test
    public void buildSimpleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("acttt"));
        l.add(new Sequence("tcttg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build().toString(), "acttg");
    }

    @Test
    public void buildNoGapMultipleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("acttt"));
        l.add(new Sequence("tattg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build().toString(), "aattg");
    }

    @Test
    public void buildGapMultipleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("a-ttt"));
        l.add(new Sequence("t-ttg"));
        l.add(new Sequence("atgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build().toString(), "atttg");
    }

    @Test
    public void buildGapSimpleTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("a-ttt"));
        l.add(new Sequence("tcttg"));
        l.add(new Sequence("acgag"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build().toString(), "acttg");
    }

    @Test
    public void buildLargerTest()
    {
        ArrayList<Sequence> l = new ArrayList<Sequence>();
        l.add(new Sequence("-------------------ttgcg"));
        l.add(new Sequence("--------------atcggtc---"));
        l.add(new Sequence("--------------atcgtgcaa-"));
        l.add(new Sequence("----taaccgcagattcc------"));
        l.add(new Sequence("--------------atcgtgcaa-"));
        l.add(new Sequence("actttacg----------------"));
        l.add(new Sequence("tact--------------------"));

        Consensus c = new Consensus(null);
        c.setAlignment(l);

        assertEquals(c.build().toString(), "aacttaaccgcagaatcgttcaag");
    }
}
