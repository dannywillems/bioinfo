package be.ac.umons.bioinfo.sequence;

import org.junit.Ignore;
import static junit.framework.TestCase.assertEquals;

import org.junit.Test;
import java.util.HashMap;
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
}
