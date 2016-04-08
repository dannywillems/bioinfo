package be.ac.umons.bioinfo.sequence;


public class Arc
{
    public final Sequence s1;
    public final boolean s1Comp;

    public final Sequence s2;
    public final boolean s2Comp;

    public final Sequence s1Aligned;
    public final Sequence s2Aligned;

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
