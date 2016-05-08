package be.ac.umons.bioinfo.greedy;

import java.util.List;

import be.ac.umons.bioinfo.sequence.*;

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

    public final int match;
    public final int mismatch;
    public final int gap;

    public final boolean inside;

    /**
     * A representation of a transition between two sequences.
     *
     * @param start        The canonical representation of the sequence at the origin of the arc.
     * @param startComp    true if the complement of start is used instead of the sequence itself; false otherwise.
     * @param end          The canonical representation of the sequence at the end of the arc.
     * @param endComp      true if the complement of end is used instead of the sequence itself; false otherwise.
     * @param score        length of the longest suffix-prefix or prefix-suffix.
     */
    public Arc(Sequence start,
               boolean startComp,
               Sequence end,
               boolean endComp,
               int score,
               int match,
               int mismatch,
               int gap,
               boolean inside)
    {
        this.start = start;
        this.startComp = startComp;
        this.end = end;
        this.endComp = endComp;
        this.score = score;
        this.match = match;
        this.mismatch = mismatch;
        this.gap = gap;

        this.inside = inside;
    }

    /**
     * @return
     */
    public SequenceAlignment getAlignment()
    {

        Sequence a, b;

        if(startComp) a = start.complement();
        else a = start;

        if(endComp) b = end.complement();
        else b = end;

        List<SequenceAlignment> candidates = Sequence.semiGlobalAlignment(a, b, match, mismatch, gap);

        return candidates.get(0);
    }

    public String toString()
    {
        SequenceAlignment alignment = getAlignment();
        return "Arc(" +
                start + ", " +
                startComp + ", " +
                end + " , " +
                endComp + ", " +
                alignment.score + ")";
    }
}
