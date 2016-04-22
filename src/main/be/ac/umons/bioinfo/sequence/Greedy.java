package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.UnionFind;
import be.ac.umons.bioinfo.parallel.Pair;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * An algorithm that proposes a list of sequence alignments that forms an Hamiltonian path.
 */
public class Greedy
{
    /**
     * @param sequences the sequences for which the arcs must be generated.
     * @param match match cost
     * @param mismatch mismatch cost
     * @param gap gap cost
     * @return a list of sequence arcs.
     */
    private List<Arc> generateArcs(List<Sequence> sequences, int match, int mismatch, int gap)
    {
        List<Arc> arcs = new ArrayList<>();

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

        return arcs;
    }

    /**
     * @param sequences the sequences for which the arcs must be generated.
     * @param match match cost
     * @param mismatch mismatch cost
     * @param gap gap cost
     * @return a list of sequence arcs.
     */
    private List<Arc> parallelGenerateArcs(List<Sequence> sequences,
                                           int match, int mismatch, int gap)
    {
        List<Pair<Sequence>> pairs = new ArrayList<Pair<Sequence>>();

        for(int i=0 ; i<sequences.size()-1 ; i++)
            for( int j=i+1 ; j<sequences.size() ; j++)
                pairs.add(new Pair<Sequence>(sequences.get(i), sequences.get(j)));

        return pairs.parallelStream()
                    .map(p -> p.a.arcGenerator(p.b, match, mismatch, gap))
                    .flatMap(Collection::stream)
                    .collect(Collectors.toList());
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
    public List<SequenceAlignment> computePath(List<Sequence> sequences, int match, int mismatch, int gap)
    {
        //List<Arc> arcs = generateArcs(sequences, match, mismatch, gap);
        List<Arc> arcs = parallelGenerateArcs(sequences, match, mismatch, gap);

        Collections.sort(arcs, new ArcComparator());

        return greedy(sequences, arcs, match, mismatch, gap);
    }

    private List<SequenceAlignment> greedy(List<Sequence> sequences,
                                           List<Arc> arcs,
                                           int match,
                                           int mismatch,
                                           int gap)
    {
        LinkedList<Arc> lArcs = new LinkedList<Arc>();
        lArcs.addAll(arcs);

        UnionFind<Sequence> groups = new UnionFind<>(sequences);
        Set<Sequence> entered = new HashSet<>();
        Set<Sequence> exited = new HashSet<>();
        Set<Arc> accepted = new HashSet<>();

        Map<Sequence, Arc> right = new HashMap<>(); // Arcs at the right of a sequence
        Map<Sequence, Arc> left = new HashMap<>(); // Arcs at the left of a sequence
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
                right.put(candidate.start, candidate);
                left.put(candidate.end, candidate);
                comp.put(candidate.start, candidate.startComp);
                comp.put(candidate.end, candidate.endComp);
            }
        }

        // Select a pivot arc, and look for all the arcs at the right of the pivot
        LinkedList<Arc> path = new LinkedList<>();

        Arc pivot = accepted.iterator().next();
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

        return path.parallelStream()
                   .map(arc -> arc.getAlignment(match, mismatch, gap))
                   .collect(Collectors.toList());
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

        Boolean _s1Comp = a.startComp;
        Sequence s = a.start;

        Boolean _s2Comp = a.endComp;
        Sequence t = a.end;

        Boolean test1 = !comp.containsKey(s)||comp.get(s) == _s1Comp;
        Boolean test2 = !comp.containsKey(t)||comp.get(t) == _s2Comp;

        return !exited.contains(a.start)
                && !entered.contains(a.end)
                && !groups.sameSet(a.start, a.end)
                && test1
                && test2;
    }
}
