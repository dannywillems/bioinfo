package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.UnionFind;

import java.util.*;

/**
 * An algorithm that proposes a list of sequence alignments that forms an Hamiltonian path.
 */
public class Greedy
{
    /**
     * Determine a list of aligment, that corresponds to the succession of alignements that
     * have been ordered aligned in order to optimize the alignment.
     * @param sequences The primary sequences to work on.
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return A sequence of alignement.
     */
    public List<SequenceAlignment> computePath(List<Sequence> sequences, int match, int mismatch, int gap)
    {
        LinkedList<Arc> arcs = new LinkedList<>();
        UnionFind<Sequence> groups = new UnionFind<>(sequences);
        Set<Sequence> entered = new HashSet<>();
        Set<Sequence> exited = new HashSet<>();
        Set<Arc> accepted = new HashSet<>();

        //first computes all the arc of the graph
        for(int i=0 ; i < sequences.size()-1 ; i++)
        {
            for(int j=i+1 ; j < sequences.size() ; j++)
            {
                Sequence s1 = sequences.get(i);
                Sequence s2 = sequences.get(j);

                List<Arc> result = s1.arcGenerator(s2, match, mismatch, gap);
                arcs.addAll(result);
            }
        }


        //orders the arcs by their score
        Collections.sort(arcs, new ArcComparator());


        Map<Sequence, Arc> right = new HashMap<>(); // Arcs at the right of a sequence
        Map<Sequence, Arc> left = new HashMap<>(); // Arcs at the left of a sequence
        Map<Sequence, Boolean> comp = new HashMap<>();

        //computes the greedy algorithm
        while(groups.size() > 1 )
        {
            Arc candidate = arcs.pop();


            Boolean bool = isAcceptable(candidate, entered, exited, groups,comp);

            if(bool)
            {
                accepted.add(candidate);
                exited.add(candidate.s1);
                entered.add(candidate.s2);
                groups.union(candidate.s1, candidate.s2);
                right.put(candidate.s1, candidate);
                left.put(candidate.s2, candidate);
                comp.put(candidate.s1, candidate.s1Comp);
                comp.put(candidate.s2, candidate.s2Comp);

            }
        }

        //Map<Sequence, Arc> right = new HashMap<>(); // Arcs at the right of a sequence
        //Map<Sequence, Arc> left = new HashMap<>(); // Arcs at the left of a sequence
        /*
        for(Arc a : accepted)
        {
            //right.put(a.s1, a);
            //left.put(a.s2, a);


        }
       */



        // Select a pivot arc, and look for all the arcs at the right of the pivot
        LinkedList<Arc> path = new LinkedList<>();

        Arc pivot = accepted.iterator().next();
        path.add(pivot);

        Sequence current = pivot.s2;

        while(right.containsKey(current))
        {
            Arc next = right.get(current);
            path.addLast(next);
            current = next.s2;
        }

        // Now, let's have a look at the left of the pivot
        current = pivot.s1;

        while(left.containsKey(current))
        {
            Arc previous = left.get(current);
            path.addFirst(previous);
            current = previous.s1;
        }

        List<SequenceAlignment> ret = new ArrayList<>(path.size()+1);

        for(Arc a : path) {
            SequenceAlignment seq = new SequenceAlignment(a.s1Aligned, a.s2Aligned, a.score);
            ret.add(seq);
        }


        return ret;
    }

    /**
     * Determines if an arc is acceptable for building an Hamiltonien path.
     * @param a an arc
     * @param entered the sequences that have been entered
     * @param exited the exited that have been exited
     * @param groups the sequences are stand together
     * @return true if a is acceptable, false otherwise
     */
    public static boolean isAcceptable(Arc a,
                                       Set<Sequence> entered,
                                       Set<Sequence> exited,
                                       UnionFind<Sequence> groups,
                                       Map<Sequence,Boolean> comp)
    {

        Boolean _s1Comp = a.s1Comp;
        Sequence s = a.s1;

        Boolean _s2Comp = a.s2Comp;
        Sequence t = a.s2;

        Boolean test1 = !comp.containsKey(s)||comp.get(s) == _s1Comp;
        Boolean test2 = !comp.containsKey(t)||comp.get(t) == _s2Comp;

        return !exited.contains(a.s1)
                && !entered.contains(a.s2)
                && !groups.sameSet(a.s1, a.s2)
                && test1
                && test2;
    }
}
