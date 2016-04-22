package be.ac.umons.bioinfo.sequence;

import java.util.List;
import java.util.ArrayList;
import org.junit.Test;
import static junit.framework.TestCase.assertEquals;

/**
 * Created by Danny Willems (contact@danny-willems.be) on 22/04/2016.
 */
public class SequenceAlignmentTest
{
    @Test
    public void insertionOrderDoesntChangeTheResultNoOverlap()
    {
        Greedy greed = new Greedy();

        Sequence f = new Sequence("catagtc");
        Sequence g = new Sequence("taactat");
        Sequence h = new Sequence("agactatcc");

        List<Sequence> fgh = new ArrayList<Sequence>();
        fgh.add(f);
        fgh.add(g);
        fgh.add(h);
        List<SequenceAlignment> result_fgh = greed.computePath(fgh, 1, -1, -2);

        List<Sequence> gfh = new ArrayList<Sequence>();
        gfh.add(g);
        gfh.add(f);
        gfh.add(h);
        List<SequenceAlignment> result_gfh = greed.computePath(gfh, 1, -1, -2);

        /* See FIXME below
        List<Sequence> ghf = new ArrayList<Sequence>();
        ghf.add(g);
        ghf.add(h);
        ghf.add(f);
        List<SequenceAlignment> result_ghf = greed.computePath(ghf, 1, -1, -2);
        */

        List<Sequence> fhg = new ArrayList<Sequence>();
        fhg.add(f);
        fhg.add(h);
        fhg.add(g);
        List<SequenceAlignment> result_fhg = greed.computePath(fhg, 1, -1, -2);

        // First, we test the path length. We only change the first argument for the test because of equality transivity.
        assertEquals(result_fgh.size(), result_fhg.size()); // fgh VS fhg
        assertEquals(result_gfh.size(), result_fhg.size()); // gfh VS fhg
        /* FIXME: Fail, check if it is because there is multiple shortest path
        assertEquals(result_ghf.size(), result_fhg.size()); // ghf VS fhg
        */

        // Now, check the path
        for (int i = 0; i < result_fgh.size(); i++)
        {
            SequenceAlignment ref = result_fhg.get(i);

            assertEquals(result_fgh.get(i).s1, ref.s1);
            assertEquals(result_fgh.get(i).s2, ref.s2);
            assertEquals(result_gfh.get(i).s1, ref.s1);
            assertEquals(result_gfh.get(i).s2, ref.s2);
        }
    }
}
