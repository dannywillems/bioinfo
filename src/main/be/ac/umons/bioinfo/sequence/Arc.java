package be.ac.umons.bioinfo.sequence;


public class Arc {

    public final Sequence start; // The initial source sequence
    public final boolean startComp; //True if we need to use the complementary of start
    public final Sequence end; // The initial target sequence
    public final boolean endComp; //True if we need to use the complementary of end
    public final Sequence alignedStart; // The start sequence when aligned with end. Can be the complementary.
    public final Sequence alignedEnd; //he end sequence when aligned with start. Can be the complementary
    public final int score; //The alignment score

    /**
     * A representation of a transition between two sequences.
     *
     * @param start        The canonical representation of the sequence at the origin of the arc.
     * @param startComp    true if the complement of start is used instead of the sequence itself; false otherwise.
     * @param end          The canonical representation of the sequence at the end of the arc.
     * @param endComp      true if the complement of end is used instead of the sequence itself; false otherwise.
     * @param alignedStart the actually used start sequence after its alignment with the actual end.
     * @param alignedEnd   the actually used end sequence after its alignement with the actual start.
     * @param score        length of the longest suffix-prefix or prefix-suffix.
     */
    public Arc(Sequence start,
               boolean startComp,
               Sequence end,
               boolean endComp,
               Sequence alignedStart,
               Sequence alignedEnd,
               int score) {
        this.start = start;
        this.startComp = startComp;
        this.end = end;
        this.endComp = endComp;
        this.alignedStart = alignedStart;
        this.alignedEnd = alignedEnd;
        this.score = score;
    }

}
