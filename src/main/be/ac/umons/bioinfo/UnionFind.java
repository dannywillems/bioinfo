package be.ac.umons.bioinfo;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A structure for efficiently merging sets.
 * See https://fr.wikipedia.org/wiki/Union-Find for details about the implementation.
 */
public class UnionFind<T>
{
    /** Maps an element to his parent which is the representative. */
    private Map<T, T> parent;
    /** Maps an element to his rank which is the height in the tree */
    private Map<T, Integer> rank;
    /** Number of equivalence classes */
    private int nbSets;

    /**
     * @param elements the elements that can potentially be merged.
     */
    public UnionFind(List<T> elements)
    {
        this.parent = new HashMap<>();
        this.rank = new HashMap<>();

        for(T e : elements)
        {
            this.parent.put(e, e);
            this.rank.put(e, 0);
        }

        this.nbSets = elements.size();
    }

    /**
     * Merges the sets that contain two elements.
     * @param x an element.
     * @param y an other element.
     */
    public void union(T x, T y)
    {
        T xRoot = find(x);
        T yRoot = find(y);

        if(!xRoot.equals(yRoot))
        {
            final Integer xRank = this.rank.get(xRoot);
            final Integer yRank = this.rank.get(yRoot);

            if(xRank < yRank)
            {
                this.parent.put(xRoot, yRoot);
            }
            else
            {
                if(xRank > yRank)
                {
                    this.parent.put(yRoot, xRoot);
                }
                else
                {
                    this.parent.put(yRoot, xRoot);
                    this.rank.put(xRoot, xRank + 1);
                }
            }

            this.nbSets--;
        }
    }

    public boolean sameSet(T x, T y)
    {
        return find(x).equals(find(y));
    }

    /**
     * Finds the element representing the set of an element.
     * @param e an element.
     * @return the element representing the set containing e.
     */
    private T find(T e)
    {
        final T parent = this.parent.get(e);

        if(!parent.equals(e))
            this.parent.put(e, find(parent));

        return this.parent.get(e);
    }

    public int size()
    {
        return nbSets;
    }
}
