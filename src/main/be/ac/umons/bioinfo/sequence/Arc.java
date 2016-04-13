package be.ac.umons.bioinfo.sequence;


public class Arc
{

    public final Sequence s1; // The initial source sequence
    public final boolean s1Comp; //True if we need to use the complementary of s1
    public final Sequence s2; // The initial target sequence
    public final boolean s2Comp; //True if we need to use the complementary of s2
    public final Sequence s1Aligned; // The s1 sequence when aligned with s2. Can be the complementary.
    public final Sequence s2Aligned; //he s2 sequence when aligned with s1. Can be the complementary
    public final int score; //The alignment score

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
