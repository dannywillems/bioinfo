package be.ac.umons.bioinfo.sequence;

import java.lang.Integer;
import java.util.ArrayList;

public class SequenceAbstract
{
    public Sequence initial;
    public int[] nb_gaps;
    private int offset;

    /* initial sequence is useless, if no mutation, the initial sequence is the
     * aligned sequence where gaps has been removed
     */
    public SequenceAbstract(Sequence initial, int[] nb_gaps)
    {
        this.initial = initial;
        this.nb_gaps = nb_gaps;
        this.offset = 0;
    }

    /* initial sequence is useless, if no mutation, the initial sequence is the
     * aligned sequence where gaps has been removed. Only used to defined the
     * size of the gaps array but can be replaced by an arraylist.
     */
    public SequenceAbstract(Sequence initial, Sequence aligned)
    {
        this.initial = initial;
        this.nb_gaps = new int[initial.getSize()];
        int j, i = 0;
        int count = 0;
        String s = aligned.toString();
        while (i < aligned.getSize())
        {
            j = 0;
            while (i + j + 1 < aligned.getSize() && aligned.getBaseAsByte(i + j + 1) == (byte) Sequence.GAP)
                j++;
            nb_gaps[count] = j;
            count++;
            i = i + j + 1;
        }

        this.offset = 0;
    }

    /**
     * Create a new abstract sequence from an aligned sequence (or not) eg
     * derive the initial sequence and compute the gaps array.
     */
    public SequenceAbstract(Sequence s)
    {
        this.offset = 0;
        while (s.getBaseAsByte(this.offset) == (byte) Sequence.GAP)
            this.offset++;

        StringBuilder str = new StringBuilder();
        int j;
        int i = this.offset;
        ArrayList<Integer> gaps = new ArrayList<Integer>();
        while (i < s.getSize())
        {
            str.append(s.getLetter(i));
            j = 0;
            while (i + j + 1 < s.getSize() && s.getBaseAsByte(i + j + 1) == (byte) Sequence.GAP)
                j++;
            gaps.add(j);
            i = i + j + 1;
        }

        this.nb_gaps = new int[gaps.size()];
        for(int k = 0;k < gaps.size();k++)
            this.nb_gaps[k] = gaps.get(k);

        this.initial = new Sequence(str.toString());
    }

    /**
     * Add [nb] gaps before the element at position [pos]. If position is
     * greater than the number of element in the sequence, gaps are added at the
     * end, ie in the nb_gaps[sequence_size - 1].
     * If position is before the offset or at 0, gaps are added as offset.
     * Examples:
     *  - addGaps(1, 1) on the sequence
     *      attgc ==> a-ttgc
     *  - addGaps(1, 2) on the sequence
     *      at--tgc ==> at---tgc
     *  - addGaps(3, 5) on the sequence
     *      at--tgc ==> at--t---gc
     *  - addGaps(3, 0) on the sequence
     *      at--tgc ==> ---at--tgc
     * @param nb number of gaps to add
     * @param pos where the gaps must be added
     */
    public void addGaps(int nb, int pos)
    {
        int i = 0;
        int real_pos = this.offset;
        // if we need to add gaps before the beginning of the sequence
        if (pos <= real_pos)
            this.offset += nb;
        else // else, we looking for the position where we must add
        {
            while (i < this.nb_gaps.length)
            {
                // if the current real position (which is a nucleotide, not a
                // gap) is exactly where we must add the gaps, we must add gaps
                // before this nucleotide, so add gaps at the position i - 1.
                if (real_pos == pos)
                {
                    this.nb_gaps[i - 1] += nb;
                    break;
                }
                // else if the next real position (ie real_pos + nb_gaps[i] + 1)
                // is before where we need to add gaps, we go to the next
                // position.
                else if (real_pos + this.nb_gaps[i] < pos)
                {
                    real_pos += this.nb_gaps[i] + 1;
                    i++;
                }
                // else, we add the number of gaps we need.
                else
                {
                    this.nb_gaps[i] += nb;
                    break;
                }
            }
            if (i == this.nb_gaps.length)
                this.nb_gaps[i - 1] += nb;
        }
    }

    public String toString()
    {
        StringBuilder s = new StringBuilder();
        for(int i = 0;i < this.offset;i++)
            s.append("-");
        for(int i = 0;i < this.initial.getSize();i++)
        {
            s.append(this.initial.getLetter(i));
            for(int j = 0;j < this.nb_gaps[i];j++)
                s.append("-");
        }
        return (s.toString());
    }

    public int getOffset()
    {
        return (this.offset);
    }

    public void setOffset(int offset)
    {
        this.offset = offset;
    }
}
