package be.ac.umons.bioinfo.sequence;


public class Arc
{
    /**
     * The initial source sequence
     */
    public final Sequence s1;
    /**
     * True if we need to use the complementary of s1
     */
    public final boolean s1Comp;

    /**
     * The initial target sequence
     */
    public final Sequence s2;
    /**
     * True if we need to use the complementary of s2
     */
    public final boolean s2Comp;

    /**
     * The s1 sequence when aligned with s2. Can be the complementary depending
     * on the value of s1Comp
     */
    public final Sequence s1Aligned;
    /**
     * The s2 sequence when aligned with s1. Can be the complementary depending
     * on the value of s2Comp
     */
    public final Sequence s2Aligned;

    /**
     * The alignment score
     */
    public final int score;

    public Arc(Sequence s1,
               boolean s1Comp,
               Sequence s2,
               boolean s2Comp,
               Sequence s1Aligned,
               Sequence s2Aligned,
               int score)
    {
        this.s1 = s1;
        this.s1Comp = s1Comp;
        this.s2 = s2;
        this.s2Comp = s2Comp;
        this.s1Aligned = s1Aligned;
        this.s2Aligned = s2Aligned;
        this.score = score;
    }
}
