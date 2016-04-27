package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.Pair;
import be.ac.umons.bioinfo.UnionFind;

import java.util.*;
import java.util.stream.Collectors;

/**
 * An algorithm that proposes a list of sequence alignments that forms an Hamiltonian path.
 */
public class Greedy
{

    public static List<Pair<Sequence>> generateAllPairs(List<Sequence> sequences)
    {
        List<Pair<Sequence>> ret = new ArrayList<>(sequences.size()*sequences.size()/2);

        for(int i=0 ; i<sequences.size()-1 ; i++)
            for(int j=i+1 ; j<sequences.size() ; j++)
                ret.add(new Pair<Sequence>(sequences.get(i), sequences.get(j)));

        return ret;
    }

    /**
     * Determine a list of aligment, that corresponds to the succession of alignements that
     * have been ordered aligned in order to optimize the alignment.
     * @param sequences The primary sequences to work on.
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return A sequence of alignement.
     */
    public static List<SequenceAlignment> greedy(List<Sequence> sequences, int match, int mismatch, int gap)
    {
        List<Pair<Sequence>> pairs = generateAllPairs(sequences);

        List<Arc> arcs = pairs.parallelStream()
                                    .map(pair -> pair.a.arcGenerator(pair.b, match, mismatch, gap))
                                    .flatMap(Collection::stream)
                                    .collect(Collectors.toList());

        Collections.sort(arcs, new ArcComparator());

        List<Arc> path = hamiltonianPath(filterArcs(arcs, sequences));

        return path.parallelStream().map(arc -> arc.getAlignment()).collect(Collectors.toList());
    }

    /**
     * @param arcs Candidates arcs that can form an hamiltonian path.
     * @return The hamiltonian path made of these arcs.
     */
    public static List<Arc> hamiltonianPath(Set<Arc> arcs)
    {
        Map<Sequence, Arc> right = right(arcs);
        Map<Sequence, Arc> left = left(arcs);

        // Select a pivot arc, and look for all the arcs at the right of the pivot
        LinkedList<Arc> path = new LinkedList<>();

        Arc pivot = arcs.iterator().next();
        path.add(pivot);

        Sequence current = pivot.end;

        while(right.containsKey(current))
        {
            Arc next = right.get(current);
            path.addLast(next);
            current = next.end;
        }

        // Now, let's have a look at the left of the pivot
        current = pivot.start;

        while(left.containsKey(current))
        {
            Arc previous = left.get(current);
            path.addFirst(previous);
            current = previous.start;
        }

        return path;
    }
    /**
     * Filters the arcs that must be retained for the Hamiltonian path.
     * @param arcs  the candidate arcs
     * @param sequences the canonical representation of the sequences contained in the arc.
     * @return the candidate arcs that must belong to the hamiltonian path
     */
    public static Set<Arc> filterArcs(List<Arc> arcs, List<Sequence> sequences)
    {
        LinkedList<Arc> lArcs = new LinkedList<Arc>();
        lArcs.addAll(arcs);

        UnionFind<Sequence> groups = new UnionFind<>(sequences);
        Set<Sequence> entered = new HashSet<>();
        Set<Sequence> exited = new HashSet<>();
        Set<Arc> accepted = new HashSet<>();

        Map<Sequence, Boolean> comp = new HashMap<>();

        //computes the greedy algorithm
        while(groups.size() > 1)
        {
            Arc candidate = lArcs.pop();

            if(isAcceptable(candidate, entered, exited, groups,comp))
            {
                accepted.add(candidate);
                exited.add(candidate.start);
                entered.add(candidate.end);
                groups.union(candidate.start, candidate.end);
                comp.put(candidate.start, candidate.startComp);
                comp.put(candidate.end, candidate.endComp);
            }
        }

        return accepted;
    }

    /**
     * Computes the arc at the right of each sequence (expect the extremum)
     * @param arcs the arcs in the hamiltonian path
     * @return the arc at the right of each
     */
    public static Map<Sequence, Arc> right(Set<Arc> arcs)
    {
        Map<Sequence, Arc> ret = new HashMap<Sequence, Arc>(arcs.size());

        for(Arc arc : arcs)
            ret.put(arc.start, arc);

        return ret;
    }

    /**
     * Computes the arc at the left of each sequence (expect the extremum)
     * @param arcs the arcs in the hamiltonian path
     * @return the arc at the left of each
     */
    public static Map<Sequence, Arc> left(Set<Arc> arcs)
    {
        Map<Sequence, Arc> ret = new HashMap<Sequence, Arc>(arcs.size());

        for(Arc arc : arcs)
            ret.put(arc.end, arc);

        return ret;
    }

    /**
     * Determines if an arc is acceptable for building an Hamiltonien path.
     * @param arc an arc
     * @param entered the canonical sequences that have been entered
     * @param exited the canonical sequences that have been exited
     * @param groups the sequences that stand together
     * @param comp The complementary with which sequences are used in the hamiltonian path
     * @return true if arc is acceptable, false otherwise
     */
    public static boolean isAcceptable(Arc arc,
                                       Set<Sequence> entered,
                                       Set<Sequence> exited,
                                       UnionFind<Sequence> groups,
                                       Map<Sequence,Boolean> comp)
    {
        Boolean test1 = (!comp.containsKey(arc.start)) || comp.get(arc.start) == (Boolean) arc.startComp;
        Boolean test2 = (!comp.containsKey(arc.end))   || comp.get(arc.end) == (Boolean) arc.endComp;

        return     !exited.contains(arc.start)
                && !entered.contains(arc.end)
                && !groups.sameSet(arc.start, arc.end)
                && test1
                && test2;
    }
}
