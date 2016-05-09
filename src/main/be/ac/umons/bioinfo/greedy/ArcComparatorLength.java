package be.ac.umons.bioinfo.greedy;

import java.util.Comparator;

/**
 * Compares two arcs according to their respective scores,
 * so that the arcs are ordered by decreasing score.
 */
public class ArcComparatorLength implements Comparator<Arc>
{
    @Override
    public int compare(Arc a1, Arc a2)
    {
        if(a1.equals(a2)) return 0;

        if(a1.score > a2.score) return -1;
        if(a1.score < a2.score) return 1;

        if (a1.longestCommon > a2.longestCommon)
            return (-1);
        else if (a1.longestCommon < a2.longestCommon)
            return (1);
        else
        {
            if (a1.start.hashCode() > a2.start.hashCode())
                return -1;
            if (a1.start.hashCode() < a2.start.hashCode())
                return 1;

            if (a1.end.hashCode() > a2.end.hashCode())
                return -1;
            if (a1.end.hashCode() < a2.end.hashCode())
                return 1;

            return (a2.start.toString().hashCode() - a1.start.toString().hashCode());
        }
    }
}
