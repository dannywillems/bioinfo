package be.ac.umons.bioinfo.sequence;


import java.util.List;

public class Arc
{
    /** The initial source sequence */
    public final Sequence start;
    /** True if we need to use the complementary of start */
    public final boolean startComp;
    /** The initial target sequence */
    public final Sequence end;
    /** True if we need to use the complementary of end */
    public final boolean endComp;
    /** The alignment score */
    public final int score;
    private final boolean startBefore;

    /**
     * A representation of a transition between two sequences.
     *
     * @param start        The canonical representation of the sequence at the origin of the arc.
     * @param startComp    true if the complement of start is used instead of the sequence itself; false otherwise.
     * @param end          The canonical representation of the sequence at the end of the arc.
     * @param endComp      true if the complement of end is used instead of the sequence itself; false otherwise.
     * @param score        length of the longest suffix-prefix or prefix-suffix.
     * @param startBefore
     */
    public Arc(Sequence start,
               boolean startComp,
               Sequence end,
               boolean endComp,
               int score,
               boolean startBefore)
    {
        this.start = start;
        this.startComp = startComp;
        this.end = end;
        this.endComp = endComp;
        this.score = score;
        this.startBefore = startBefore;
    }

    /**
     *
     * @param match
     * @param mismatch
     * @param gap
     * @return
     */
    public SequenceAlignment getAlignment(int match, int mismatch, int gap)
    {
        Sequence a, b;

        if(startComp) a = start.complement();
        else a = start;

        if(endComp) b = end.complement();
        else b = end;

        List<SequenceAlignment> candidates = Sequence.semiGlobalAlignment(a, b, match, mismatch, gap);

        if(startBefore) return candidates.get(0);
        else return new SequenceAlignment(  candidates.get(1).s2,
                                            candidates.get(1).s1,
                                            candidates.get(1).score);
    }
}
