package be.ac.umons.bioinfo.sequence;

import java.util.ArrayList;
import java.util.List;

public class SequenceAlignmentAbstract
{
    public SequenceAbstract s;
    public SequenceAbstract t;

    public SequenceAlignmentAbstract(SequenceAbstract s, SequenceAbstract t)
    {
        this.s = s;
        this.t = t;
    }

    public SequenceAlignmentAbstract(SequenceAlignment sa)
    {
        this.s = new SequenceAbstract(sa.initial_s1, sa.s1);
        this.t = new SequenceAbstract(sa.initial_s2, sa.s2);
    }

    public SequenceAlignmentAbstract(Sequence s, Sequence t)
    {
        this.s = new SequenceAbstract(s);
        this.t = new SequenceAbstract(t);
    }

    public void showAlignment()
    {
        System.out.println(this.s.toString());
        System.out.println(this.t.toString());
    }
}
