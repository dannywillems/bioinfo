package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.greedy.*;

import java.util.*;

import org.junit.Test;
import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertArrayEquals;

public class SequenceAbstractTest
{
    /* ---------------------------------------------------------------------- */
    @Test
    public void constructorTest()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("--actg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(s.toString(), aligned.toString());
    }

    @Test
    public void constructorTestBig1()
    {
        Sequence aligned = new Sequence("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------cttgggcctcctggtcgagg-ccatcagccgac-cggcctcgtg--gacggttg-gg-accggtaacgggcccccagcgcggcccgggatcccgagggggaggggtatgtggccagcgccgacgggccaggagagcccatcgtggttcgggcctccacgaagcgagaggcctaccgggagggccaggaaggcgtgggtgcggggcttctcgagtagctgggcccgtgccctggggcactatctcctgcccctcttggactaccggatccgcg-ggagtgg-aacgtggaccggggaagg-ct-cctttcctgggc--gt-ac-ag-c-tg-gaatggctggacc-ctg-atgaagggc-g-ct-t-c-g-t-c---t---g-g---g---t-g-g-ccc-ggt-aga-caccgcagga-g--g--c--c--t--a--c--c-t-a-----c----a----c-----c----------g---------g----g-----a-----c----c----g----t----g-----a----a--c----c---g---a---a--t--c---g---c-----c--c--g---t-----c-----a--t--g--a-----g-----c--c--c--a-----t----g--c--g--c--c----a--c---c-----t--t----t--c---c-----c--c---c---a---c--c---t------c---g---g--g---g--g-----a---g--a----c---c--t-----c--a--a----c--c----g--t--g--c--c--c----g--g--g--a-----c--g--a--c----c----t-----c--c--t--t----g----c---t----c----g---g---g---c---c----t----g---g---g---a---g---t-----a-----c---a---t---t----g-----a-----g--------g--------c---c-----t-----a--------a--------a--------g-----g-----c-----c-----g-----c---g-----a-----g-----g---c-----c---------g-----a----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
        SequenceAbstract s = new SequenceAbstract(aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void constructorTestBig2()
    {
        Sequence aligned = new
            Sequence("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------gtcctgcccaagga-ctacgccccggg-ggctaaccgctggc-t-g-g-ga-caacggcctgccggtgcgggtgacgggcg-cct-a-c-c-t---g-c-g-cgacttctgacgcctccttcctcct-g--a--a--c---c--c-c-t-----a----c----a--t----c-c----c----c--a----t-----c--a--c-----c---g---g-------c---c--c--c--a--t--c-----g--c--c--t-----g--g--g--t--g--a--g--c--t--a--c-----c--c--c--a--g--t---a-----c--a---c-----c--t--g---g---a--c------c------a--c--c---a--a--g--cg-g----g--c--t--g--c--t--c-a--a-g--g--a--c--a-c--c--a--t-c---g-g-a----a---t--a-g---g--g---a---t--c-t--t-c---c---g---g---c---g--g-t-------g-----t--t--c---t--c----t-a---c--------g--a-----c----a-g------c-----g-c--a-----c---a---g---a-----t---a--c---------g-----t---a--g-----g-----g-----g---t---c-g---t---c---c--t--a-----c--t-----a-----c-----c-----t--c--c---a-c-----t---------a--------g---c-----c--------c-c----g--g---c----c---t----a--------c--c-----g-----c-------t----a--c-----c-g-----g-------g-----t---------g----g----t---g-----g-----------a--g--t-----a--c------t---g--g----g-----a-----------c-------t-g-------g----t--c--c-g---c-----c-------a----a----t---t----g--c----c--a--c---g-----g-----t----g--g--a---c----a----a--c---g---g-----g-----a-----a-------g----g-----a---g-----g--g----t-c----t--t--t----c----c--c----------c---c--c----g---t----c-------c---------t-----t----t---g----g------t------c----c----a---t----t--------g---a-----------g------a-----------t-----------t--t----a---c------------c--------t-c----a----g---c-----g-------a------c----t----a------c-------t---------a---------c-------g---------a-------t------c----t-----c---t--c----c--a---a----c----g-----g-------------g---g-----t-------g---------t----g-g-----c-----------a----a-t---c--a--g-----------a-----c----c---------g---c----t------a-------g--------g---g----g-----g------------t--------g----------------g----t------------------g---------a--------a--------g-------g-----g-------a----------t----g-----c------g-------t---------a----a--a----g-------------------c------t----c------g-----t----------c-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------^C-----------------------------------t----t--c---c-t-c-ag---c-------------a----c---t--a------c-------t---------a---------cg---------a-------t--c-t----c---tc----c--a---a----c----g-----g-------------g---g-----t-------g---------t----g-g----c-------g---c----a-t------a-g----t-----c----c-g---ct------aa--------g---g---g-----g------------t--------g----------------g----t------------------g---------a--------a--------g-------g----g------------------t---g-----c------g-------t-------------a--a----g-------------------c------t----c------g---t----------c--g----g-------c-----------c----t----t------g---g----------c-----------c-------t----g----c--------t---g---g-------c--------g-----------t-------c------g----c--------t---g---g-----------c--------g----------g----------t-----g---c-----t---g-------t--t--------t----------t----------c----g----------g----c-----c---t---------g--t------g-------g-------t---g----------a----a----t-----g---------g----g------t-------c-------t--t----t-----g------t-----t------c----t-----t----a-------c----t----c-----c------c-----a------t----a-----g----a------t-----------c---t-----t----c----a------t----g----t-----g-------t-----------t---g-----t----------a----c---t----a--------c-----c-------c------c---c----------a-------------g---a-------a--------t------c----c------c--a----------g-------c--c--t-------c----------t----g----c----c--c-------c----g-----g----a---------g------c----------a----g-------c----g--g-------g----a--g----c---------g------g-------t--g-------a------g--c---a----c-----c------a----a---c---a--a--c------c------c-------a-----g------g----g---t---t----c--t----a------c--a----g-----c-----t----g---c--c--c----c---t-g------a--t---a-----g----c--c----cc--c--a-------c--c------t-----c----c----c--c-----tc--c--c-------a--c--g--c-------c----c-c--c---------t---g-------g--g-------c---c--g---g---g-------g----g--a--g---------a---t-------c---a-------t---t----t---c------t---t---t---t--c--c---g--a---c--t-----c----c--t------a---g--t---t--a--t---g--g--g---c-----c---c---g---t---c---a-----c---g---c---t--g--c-------t-------a---c---c--------t-------t---t-----c---t---a--------c----a---c-------g-------a---g---c--------g-------c---c----g---c---c--------c----a---g-------a-------g---a---a--------g----t---a--c-------g---a----c---c---t--------g-------a---c---a--------c----t-t--c-------t--a----c---a---g--------c----a--a--c---g--c-c-g-g---g-g--a--g-g--g--c----cg--c-------c-a----c---c-----t--------g-----ac----c-a----c-------c-c-------c---c-----g--------g-----a--------c-----c---c----c-c---g-------t-g----g---c-----c--------g-----c---c--a---c---gg--g-a---t---t--a---t--t--c---g----t-c--g-a---a---c---a-----a--t-----c-a----g-a----a-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
        SequenceAbstract s = new SequenceAbstract(aligned);
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void constructorStringTest1()
    {
        SequenceAbstract s1 = new SequenceAbstract("tact-----");
        SequenceAbstract s2 = new SequenceAbstract(new Sequence("tact-----"));
        SequenceAbstract s3 = new SequenceAbstract(new Sequence("tact"), new Sequence("tact-----"));

        assertEquals(9, s1.getSize());
        assertEquals(9, s2.getSize());
        assertEquals(9, s3.getSize());
        assertEquals(s1, s2);
        assertEquals(s1, s3);
        assertEquals(s2, s3);
    }

    @Test
    public void constructorStringTest2()
    {
        SequenceAbstract s1 = new SequenceAbstract("actttacg---");
        SequenceAbstract s2 = new SequenceAbstract(new Sequence("actttacg---"));
        SequenceAbstract s3 = new SequenceAbstract(new Sequence("actttacg"), new Sequence("actttacg---"));

        assertEquals(11, s1.getSize());
        assertEquals(11, s2.getSize());
        assertEquals(11, s3.getSize());
        assertEquals(s1, s2);
        assertEquals(s1, s3);
        assertEquals(s2, s3);
    }

    @Test
    public void constructorStringTest3()
    {
        SequenceAbstract s1 = new SequenceAbstract("------cgcaa");
        SequenceAbstract s2 = new SequenceAbstract(new Sequence("------cgcaa"));
        SequenceAbstract s3 = new SequenceAbstract(new Sequence("cgcaa"), new Sequence("------cgcaa"));

        assertEquals(11, s1.getSize());
        assertEquals(11, s2.getSize());
        assertEquals(11, s3.getSize());
        assertEquals(s1, s2);
        assertEquals(s1, s3);
        assertEquals(s2, s3);
    }

    @Test
    public void constructorStringTest4()
    {
        SequenceAbstract s1 = new SequenceAbstract("ggaatc-tgcgagtta");
        SequenceAbstract s2 = new SequenceAbstract(new Sequence("ggaatc-tgcgagtta"));
        SequenceAbstract s3 = new SequenceAbstract(new Sequence("ggaatctgcgagtta"), new Sequence("ggaatc-tgcgagtta"));

        assertEquals(s1, s2);
        assertEquals(s1, s3);
        assertEquals(s2, s3);
    }

    @Test
    public void constructorStringTest5()
    {
        SequenceAbstract s1 = new SequenceAbstract("---atcgtgcaa----");
        SequenceAbstract s2 = new SequenceAbstract(new Sequence("---atcgtgcaa----"));
        SequenceAbstract s3 = new SequenceAbstract(new Sequence("atcgtgcaa"), new Sequence("---atcgtgcaa----"));

        assertEquals(s1, s2);
        assertEquals(s1, s3);
        assertEquals(s2, s3);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void complementTest()
    {
        SequenceAbstract t = new SequenceAbstract("cagcgtgg");
        SequenceAbstract res = t.complement();
        SequenceAbstract wanted = new SequenceAbstract("ccacgctg");

        assertEquals(wanted, res);
    }
    @Test
    public void complementOfComplementMustBeIdentity()
    {
        List<SequenceAbstract> sequences = new ArrayList<SequenceAbstract>();
        sequences.add(new SequenceAbstract("GATTACA"));
        sequences.add(new SequenceAbstract("ATCATTAGTGG"));
        sequences.add(new SequenceAbstract("A"));
        sequences.add(new SequenceAbstract("G"));
        sequences.add(new SequenceAbstract("T"));
        sequences.add(new SequenceAbstract("C"));

        for(SequenceAbstract s : sequences)
            assertEquals(s, s.complement().complement());
    }

    @Test
    public void complementWithGaps()
    {
        SequenceAbstract t = new SequenceAbstract("--ca--gc--gtgg---");
        SequenceAbstract res = t.complement();
        SequenceAbstract wanted = new SequenceAbstract("---ccac--gc--tg--");

        assertEquals(wanted, res);
    }

    @Test
    public void complementWithGapsComplementOfComplementMustBeIdentidyTest()
    {
        List<SequenceAbstract> sequences = new ArrayList<SequenceAbstract>();
        sequences.add(new SequenceAbstract("GA--TT--A-CA---"));
        sequences.add(new SequenceAbstract("-----ATCA--TTA--GTGG--"));
        sequences.add(new SequenceAbstract("A----"));
        sequences.add(new SequenceAbstract("---G---"));
        sequences.add(new SequenceAbstract("---T--"));
        sequences.add(new SequenceAbstract("----C----"));

        for(SequenceAbstract s : sequences)
            assertEquals(s, s.complement().complement());
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void arrayGapsNoGapsTest()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("actg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {0, 0, 0, 0};
        assertArrayEquals(gaps, s.getNbGaps());
    }

    @Test
    public void arrayGapsGapsTest1()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a-c-t-g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {1, 1, 1, 0};
        assertArrayEquals(gaps, s.getNbGaps());
    }

    @Test
    public void arrayGapsGapsTest2()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a---c--t-g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {3, 2, 1, 0};
        assertArrayEquals(gaps, s.getNbGaps());
    }

    @Test
    public void arrayGapsGapsTest3()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("a---c--t-g--");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        int[] gaps = {3, 2, 1, 2};
        assertArrayEquals(gaps, s.getNbGaps());
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
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
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
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
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void addGapsAndReturnIndiceNoTest()
    {
        Sequence initial = new Sequence("actg");
        Sequence aligned = new Sequence("actg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(-1, s.addGapsAndReturnIndice(0, 0));
        assertEquals(initial.toString(), s.toString());

        assertEquals(0, s.addGapsAndReturnIndice(0, 1));
        assertEquals(initial.toString(), s.toString());

        assertEquals(1, s.addGapsAndReturnIndice(0, 2));
        assertEquals(initial.toString(), s.toString());

        assertEquals(2, s.addGapsAndReturnIndice(0, 3));
        assertEquals(initial.toString(), s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceTest1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGapsAndReturnIndice(3, 1);
        assertEquals("a---ttg", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceTestTwoTimes1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(0, s.addGapsAndReturnIndice(3, 1));
        assertEquals(0, s.addGapsAndReturnIndice(3, 4));
        assertEquals("a------ttg", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceTestTwoTimes2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(0, s.addGapsAndReturnIndice(3, 1));
        assertEquals(1, s.addGapsAndReturnIndice(3, 5));
        assertEquals("a---t---tg", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceEndTest1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(3, s.addGapsAndReturnIndice(3, 4));
        assertEquals("attg---", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceEndTest2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(3, s.addGapsAndReturnIndice(3, 7));
        assertEquals("attg---", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceBeginAndEndTest2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(3, s.addGapsAndReturnIndice(3, 7));
        assertEquals(-1, s.addGapsAndReturnIndice(3, 0));
        assertEquals("---attg---", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceBeginAndEndTest3()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("attg");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(3, s.addGapsAndReturnIndice(3, 7));
        assertEquals(0, s.addGapsAndReturnIndice(3, 1));
        assertEquals("a---ttg---", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceWhenInitialGaps1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("a-tt--g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(0, s.addGapsAndReturnIndice(3, 1));
        assertEquals("a----tt--g", s.toString());
    }

    @Test
    public void addGapsAndReturnIndiceWhenInitialGaps2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("a-tt--g");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(1, s.addGapsAndReturnIndice(3, 3));
        assertEquals("a-t---t--g", s.toString());
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void addGapsAfterIndiceTest1()
    {
        Sequence initial = new Sequence("attgc");
        SequenceAbstract s = new SequenceAbstract(initial);

        s.addGapsAfterIndice(1, 1);
        assertEquals("at-tgc", s.toString());
    }

    @Test
    public void addGapsAfterIndiceTest2()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGapsAfterIndice(1, 2);
        assertEquals("at--t-gc", s.toString());
    }

    @Test
    public void addGapsAfterIndiceTest3()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGapsAfterIndice(3, 5);
        assertEquals("at--tgc---", s.toString());
    }

    @Test
    public void addGapsAfterIndiceTest4()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGapsAfterIndice(3, 0);
        assertEquals("a---t--tgc", s.toString());
    }

    @Test
    public void addGapsAfterIndiceTest5()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        s.addGapsAfterIndice(3, -5);
        assertEquals("---at--tgc", s.toString());
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void addGapsAfterIndiceAndReturnPositionTest1()
    {
        Sequence initial = new Sequence("attgc");
        SequenceAbstract s = new SequenceAbstract(initial);

        assertEquals(2, s.addGapsAfterIndiceAndReturnPosition(1, 1));
        assertEquals("at-tgc", s.toString());
    }

    @Test
    public void addGapsAfterIndiceAndReturnPositionTest2()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(5, s.addGapsAfterIndiceAndReturnPosition(1, 2));
        assertEquals("at--t-gc", s.toString());
    }

    @Test
    public void addGapsAfterIndiceAndReturnPositionTest3()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(7, s.addGapsAfterIndiceAndReturnPosition(3, 5));
        assertEquals("at--tgc---", s.toString());
    }

    @Test
    public void addGapsAfterIndiceAndReturnPositionTest4()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(1, s.addGapsAfterIndiceAndReturnPosition(3, 0));
        assertEquals("a---t--tgc", s.toString());
    }

    @Test
    public void addGapsAfterIndiceAndReturnPositionTest5()
    {
        Sequence initial = new Sequence("attgc");
        Sequence aligned = new Sequence("at--tgc");
        SequenceAbstract s = new SequenceAbstract(initial, aligned);

        assertEquals(0, s.addGapsAfterIndiceAndReturnPosition(3, -5));
        assertEquals("---at--tgc", s.toString());
    }
    /* ---------------------------------------------------------------------- */


    /* ---------------------------------------------------------------------- */
    @Test
    public void buildFromAlignedSequence1()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("a--ttg");
        SequenceAbstract s = new SequenceAbstract(aligned);

        assertEquals(aligned.getPosFirstNucleotide(), s.getOffset());
        assertEquals(initial.toString(), s.initial.toString());
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void buildFromAlignedSequence2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("--a--ttg");
        SequenceAbstract s = new SequenceAbstract(aligned);

        assertEquals(aligned.getPosFirstNucleotide(), s.getOffset());
        assertEquals(initial.toString(), s.initial.toString());
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void buildFromAlignedSequence3()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("--a--t---tg");
        SequenceAbstract s = new SequenceAbstract(aligned);

        assertEquals(aligned.getPosFirstNucleotide(), s.getOffset());
        assertEquals(initial.toString(), s.initial.toString());
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void buildFromAlignedSequence4()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("--a--t---tg----");
        SequenceAbstract s = new SequenceAbstract(aligned);

        assertEquals(aligned.getPosFirstNucleotide(), s.getOffset());
        assertEquals(initial.toString(), s.initial.toString());
        assertEquals(aligned.toString(), s.toString());
    }

    @Test
    public void buildFromAlignedSequenceSimplest()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = initial;
        SequenceAbstract s = new SequenceAbstract(aligned);

        assertEquals(aligned.getPosFirstNucleotide(), s.getOffset());
        assertEquals(initial.toString(), s.initial.toString());
        assertEquals(aligned.toString(), s.toString());
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void equals1()
    {
        Sequence initial = new Sequence("attg");
        SequenceAbstract s1 = new SequenceAbstract(initial, initial);
        SequenceAbstract s2 = new SequenceAbstract(initial, initial);

        assertTrue(s1.equals(s2));
    }

    @Test
    public void equals2()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("--attg");

        SequenceAbstract s1 = new SequenceAbstract(initial, initial);
        SequenceAbstract s2 = new SequenceAbstract(initial, aligned);
        assertFalse(s1.equals(s2));

        s1 = new SequenceAbstract(initial, aligned);
        assertTrue(s1.equals(s2));
    }

    @Test
    public void equals3()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("at--t--g");

        SequenceAbstract s1 = new SequenceAbstract(initial, initial);
        SequenceAbstract s2 = new SequenceAbstract(initial, aligned);
        assertFalse(s1.equals(s2));

        s1 = new SequenceAbstract(initial, aligned);
        assertTrue(s1.equals(s2));
    }

    @Test
    public void equals4()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("at--t--g");
        SequenceAbstract s1 = new SequenceAbstract(initial, aligned);
        SequenceAbstract s2 = new SequenceAbstract(aligned);

        assertTrue(s1.equals(s2));
    }

    @Test
    public void equals5()
    {
        Sequence initial = new Sequence("attg");
        Sequence aligned = new Sequence("---at--t--g---");
        SequenceAbstract s1 = new SequenceAbstract(initial, aligned);
        SequenceAbstract s2 = new SequenceAbstract(aligned);

        assertTrue(s1.equals(s2));
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void size1()
    {
        Sequence initial = new Sequence("attg");
        SequenceAbstract s1 = new SequenceAbstract(initial, initial);

        assertEquals(4, s1.getSize());
    }

    @Test
    public void size2()
    {
        Sequence aligned = new Sequence("att--g");
        SequenceAbstract s1 = new SequenceAbstract(aligned);

        assertEquals(6, s1.getSize());
    }

    @Test
    public void size3()
    {
        Sequence aligned = new Sequence("--att--g");
        SequenceAbstract s1 = new SequenceAbstract(aligned);

        assertEquals(8, s1.getSize());
    }

    @Test
    public void size4()
    {
        Sequence aligned = new Sequence("--att--g----");
        SequenceAbstract s1 = new SequenceAbstract(aligned);

        assertEquals(12, s1.getSize());
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Test
    public void iteratorTest()
    {
        SequenceAbstract s = new SequenceAbstract("--at--g----");

        StringBuilder s_n = new StringBuilder();
        Iterator<Character> i = s.iterator();
        while (i.hasNext())
            s_n.append(i.next().charValue());
        assertEquals(s.toString(), s_n.toString());
    }

    @Test
    public void iteratorNoGapsTest()
    {
        SequenceAbstract s = new SequenceAbstract("atg");

        StringBuilder s_n = new StringBuilder();
        Iterator<Character> i = s.iterator();
        while (i.hasNext())
            s_n.append(i.next().charValue());
        assertEquals(s.toString(), s_n.toString());
    }

    @Test
    public void iteratorRandomTest()
    {
        StringBuilder str = new StringBuilder();
        Random rdm = new Random();
        int size = rdm.nextInt(1500);
        for(int i = 0;i < size;i++)
            str.append(Nucleotide.base2letter((byte) rdm.nextInt(5)));
        SequenceAbstract s = new SequenceAbstract(str.toString());

        StringBuilder s_n = new StringBuilder();
        Iterator<Character> i = s.iterator();
        while (i.hasNext())
            s_n.append(i.next().charValue());

        assertEquals(s.toString(), s_n.toString());
    }
    /* ---------------------------------------------------------------------- */
}
