package be.ac.umons.bioinfo.sequence;

import java.util.ArrayList;

public class SequencePairSame
{
    public Sequence initial;
    public Sequence s;
    public Sequence t;

    public SequencePairSame(Sequence initial, Sequence s, Sequence t)
    {
        this.initial = initial;
        this.s = s;
        this.t = t;
    }

    public int[][] findGaps()
    {
        int i_s = 0;
        int i_t = 0;
        int i = 0;
        ArrayList<Integer> gaps_s = new ArrayList<Integer>();
        ArrayList<Integer> gaps_t = new ArrayList<Integer>();

        while (i < this.initial.getSize())
        {
            if (this.s.getBaseAsByte(i_s) == Nucleotide.GAP && this.t.getBaseAsByte(i_t) == Nucleotide.GAP)
            {
                i_s++;
                i_t++;
            }
            else if (this.s.getBaseAsByte(i_s) == Nucleotide.GAP)
            {
                if (this.t.getBaseAsByte(i_t) != Nucleotide.GAP)
                    gaps_t.add(i_t);
                i_s++;
            }
            else if (this.t.getBaseAsByte(i_t) == Nucleotide.GAP)
            {
                gaps_s.add(i_s);
                i_t++;
            }
            else
            {
                i_s++;
                i_t++;
                i++;
            }
        }

        int[][] gaps = new int[2][];

        gaps[0] = new int[gaps_s.size()];
        for(int j = 0;j < gaps_s.size();j++)
            gaps[0][j] = gaps_s.get(j);

        gaps[1] = new int[gaps_t.size()];
        for(int j = 0;j < gaps_t.size();j++)
            gaps[1][j] = gaps_t.get(j);

        return (gaps);
    }

    public void rebuildWithGaps()
    {
        this.rebuildWithGaps(this.findGaps());
    }

    public void rebuildWithGaps(int[][] gaps)
    {
        this.s = this.s.rebuildAddingGaps(gaps[0]);
        this.t = this.t.rebuildAddingGaps(gaps[1]);
    }

    public String toString()
    {
        return ("Initial sequence: " + this.initial + "\nAligned with s: " + this.s + "\nAligned with t: " + this.t);
    }
}
