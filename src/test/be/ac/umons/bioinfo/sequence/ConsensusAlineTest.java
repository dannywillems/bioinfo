package be.ac.umons.bioinfo.sequence;

import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

/**
 * Created by aline on 6/05/16.
 */
public class ConsensusAlineTest
{
    @Test
    public void findMaximumVotesTest()
    {

        Map<Integer, MyCounter> data = new HashMap();
        MyCounter counter1 = new MyCounter();
        counter1.vote((byte) Nucleotide.A);
        counter1.vote((byte) Nucleotide.A);
        counter1.vote((byte) Nucleotide.C);
        counter1.vote((byte) Nucleotide.GAP);
        counter1.vote((byte) Nucleotide.GAP);
        counter1.vote((byte) Nucleotide.GAP);

        MyCounter counter2 = new MyCounter();
        counter2.vote((byte) Nucleotide.GAP);
        counter2.vote((byte) Nucleotide.T);
        counter2.vote((byte) Nucleotide.T);
        counter2.vote((byte) Nucleotide.GAP);
        counter2.vote((byte) Nucleotide.GAP);
        counter2.vote((byte) Nucleotide.G);

        MyCounter counter3 = new MyCounter();
        counter3.vote((byte) Nucleotide.A);
        counter3.vote((byte) Nucleotide.T);
        counter3.vote((byte) Nucleotide.T);
        counter3.vote((byte) Nucleotide.GAP);
        counter3.vote((byte) Nucleotide.GAP);
        counter3.vote((byte) Nucleotide.C);

        data.put(-1,counter1);
        data.put(0,counter2);
        data.put(1,counter3);

        Sequence result = ConsensusAline.produceConsensusSequence(data);
        assertEquals(new Sequence("att"),result);

    }
}
