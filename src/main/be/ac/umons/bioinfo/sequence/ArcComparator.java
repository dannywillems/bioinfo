package be.ac.umons.bioinfo.sequence;

import java.util.Comparator;

/**
 * Compares two arcs according to their respective scores,
 * so that the arcs are ordered by decreasing score.
 */
public class ArcComparator implements Comparator<Arc>
{
    @Override
    public int compare(Arc a1, Arc a2)
    {
        if (a1.score < a2.score)
            return -1;
        else
            return (a1.score > a2.score ? 1 : 0);
    }
}
