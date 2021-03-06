package be.ac.umons.bioinfo.greedy;

import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import be.ac.umons.bioinfo.sequence.*;

import static org.junit.Assert.assertTrue;

/**
 * Created by aline on 23/04/16.
 */
public class ArcComparatorTest
{
    @Test
    public void simpleComparison()
    {
        Sequence s = new Sequence("");

        Arc a1 = new Arc(s, false, s, false, 1, 42, 42, 42, false);
        Arc a2 = new Arc(s, false, s, false, 2, 42, 42, 42, false);
        Arc a3 = new Arc(s, false, s, false, 3, 42, 42, 42, false);

        ArcComparator comparator = new ArcComparator();

        assertTrue(comparator.compare(a1, a2) > 0);
        assertTrue(comparator.compare(a2, a1) < 0);
        assertTrue(comparator.compare(a1, a1) == 0);

        assertTrue(comparator.compare(a1, a3) > 0);
        assertTrue(comparator.compare(a3, a1) < 0);
    }

    @Test
    public void sortComparison()
    {
        Sequence s = new Sequence("");

        Arc a1 = new Arc(s, false, s, false, 1, 42, 42, 42, false);
        Arc a2 = new Arc(s, false, s, false, 2, 42, 42, 42, false);
        Arc a3 = new Arc(s, false, s, false, 3, 42, 42, 42, false);
        Arc a4 = new Arc(s, false, s, false, 4, 42, 42, 42, false);
        Arc a5 = new Arc(s, false, s, false, 5, 42, 42, 42, false);

        List<Arc> list = Arrays.asList(a3, a5, a2, a1, a4);
        Collections.sort(list, new ArcComparator());

        double[] scores = list.stream().mapToDouble(a -> a.score).toArray();

        assert(Arrays.equals(scores, new double[]{5,4,3,2,1}));
    }
}
