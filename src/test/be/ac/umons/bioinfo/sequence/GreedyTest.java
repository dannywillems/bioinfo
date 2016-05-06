package be.ac.umons.bioinfo.sequence;


import org.junit.Ignore;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;



/**
 * Created by aline on 8/04/16.
 */
public class GreedyTest {

    @Test
    //@Ignore
    public void greedyCollection1Test() {

        Greedy greed = new Greedy();
        String path2 = "src/test/be/ac/umons/bioinfo/sequence/Collections/test4/collection4.fasta";
        File file = new File(path2);

        FastaReader fasta = new FastaReader();

        try
        {
            List<Sequence> list = fasta.readFromFile(file);
            List<SequenceAlignment> result = greed.greedy(list, 1, -1, -2);
            Iterator<SequenceAlignment> iterator = result.iterator();
            /*
            while (iterator.hasNext()) {
                SequenceAlignment alignment = iterator.next();
                System.out.println(alignment.s1);
                System.out.println(alignment.s2);

                System.out.println("****************");
            }
            */

            ConsensusAline consensus = new ConsensusAline(result);
            Sequence seq = consensus.buildSequence();
            FastaWriter.write("helloworld", seq,new File("./output4.fasta"),80);
            System.out.println(seq);
        } catch (IOException e) {
            System.out.println("ça ne fonctionne pas ");
            System.out.println(e);
        }
    }
}
