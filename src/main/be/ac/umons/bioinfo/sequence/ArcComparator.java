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
        else{
            if (a1.score > a2.score)
            return 1;

            else
                return 0;

        }


        //return a2.score - a1.score;
    }
}
