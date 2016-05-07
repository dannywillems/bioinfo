package be.ac.umons.bioinfo.sequence;

import java.lang.Integer;
import java.util.ArrayList;
import java.util.Iterator;

public class SequenceAbstract //implements Iterable<Character>
{
    public Sequence initial;
    public int[] nb_gaps;
    private int offset;

    /* ---------------------------------------------------------------------- */
    // CONSTRUCTORS
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
        int j, i = aligned.getPosFirstNucleotide();;
        int count = 0;
        String s = aligned.toString();
        while (i < aligned.getSize())
        {
            j = 0;
            while (i + j + 1 < aligned.getSize() && aligned.getBaseAsByte(i + j + 1) == (byte) Nucleotide.GAP)
                j++;
            nb_gaps[count] = j;
            count++;
            i = i + j + 1;
        }

        this.offset = aligned.getPosFirstNucleotide();;
    }

    /**
     * Create a new abstract sequence from an aligned sequence (or not) eg
     * derive the initial sequence and compute the gaps array.
     */
    public SequenceAbstract(Sequence s)
    {
        this.offset = 0;
        while (s.getBaseAsByte(this.offset) == (byte) Nucleotide.GAP)
            this.offset++;

        StringBuilder str = new StringBuilder();
        int j;
        int i = this.offset;
        ArrayList<Integer> gaps = new ArrayList<Integer>();
        while (i < s.getSize())
        {
            str.append(s.getLetter(i));
            j = 0;
            while (i + j + 1 < s.getSize() && s.getBaseAsByte(i + j + 1) == (byte) Nucleotide.GAP)
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
     * Create a new abstract sequence from a string.
     * FIXME: lazy to reimplement an independent method from Sequence class.
     * Must be independent!
     */
    public SequenceAbstract(String s)
    {
        this(new Sequence(s));
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /**
     * TODO
     */
    public SequenceAbstract complement()
    {
        Sequence s = new Sequence(this.toString());
        return (new SequenceAbstract(s.complement()));
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
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
     * @param pos where the gaps must be added. Absolute position eg including
     * the offset.
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

    /**
     * Add [nb] gaps before the element at position [pos] and return the indice
     * of the nucleotide which gaps has been inserted after, -1 if it's in the
     * offset. If position is greater than the number of element in the
     * sequence, gaps are added at the end, ie in the nb_gaps[sequence_size -
     * 1].
     * If position is before the offset or at 0, gaps are added as offset.
     * Examples:
     *  - addGapsAndReturnIndice(1, 1) on the sequence
     *      attgc ==> a-ttgc and returns 0
     *  - addGapsAndReturnIndice(1, 2) on the sequence
     *      at--tgc ==> at---tgc and returns 1
     *  - addGapsAndReturnIndice(3, 5) on the sequence
     *      at--tgc ==> at--t---gc and returns 2
     *  - addGapsAndReturnIndice(3, 0) on the sequence
     *      at--tgc ==> ---at--tgc and returns -1
     * @param nb number of gaps to add
     * @param pos where the gaps must be added. Absolute position eg including
     * the offset.
     * @return indice of the nucleotide which gaps has been inserted after. -1
     * if it's in the offset.
     */
    public int addGapsAndReturnIndice(int nb, int pos)
    {
        int i = 0;
        int real_pos = this.offset;
        // if we need to add gaps before the beginning of the sequence
        if (pos <= real_pos)
        {
            this.offset += nb;
            return (-1);
        }
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
                    return (i - 1);
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
                    return (i);
                }
            }
            this.nb_gaps[i - 1] += nb;
            return (i - 1);
        }
    }
    /* ---------------------------------------------------------------------- */

    /**
     * Add [nb] gaps after the nucleotide at the position [pos] in the initial
     * sequence.
     * If [pos] is negative, it inserts in the offset. To add gaps at the end,
     * use [pos] greater than or equal to nb_gaps.length - 1.
     * Examples:
     *  - addGapsAfterIndice(1, 1) on the sequence
     *      attgc ==> at-tgc
     *  - addGapsAfterIndice(1, 2) on the sequence
     *      at--tgc ==> at--t-gc
     *  - addGapsAfterIndice(3, 5) on the sequence
     *      at--tgc ==> at--tgc---
     *  - addGapsAfterIndice(3, 0) on the sequence
     *      at--tgc ==> a---t--tgc
     *  - addGapsAfterIndice(3, -5) on the sequence
     *      at--tgc ==> ---at--tgc
     * @param nb number of gaps to add
     * @param indice which the gaps must be added after.
     */
    public void addGapsAfterIndice(int nb, int indice)
    {
        if (indice < 0)
            this.offset += nb;
        else if (indice < nb_gaps.length - 1)
            this.nb_gaps[indice] += nb;
        else
            this.nb_gaps[nb_gaps.length - 1] += nb;
    }

    /**
     * Add [nb] gaps after the nucleotide at the position [pos] in the initial
     * sequence and returns the absolute position where gaps were inserted.
     * If [pos] is negative, it inserts in the offset. To add gaps at the end,
     * use [pos] greater than or equal to nb_gaps.length - 1.
     * Examples:
     *  - addGapsAfterIndice(1, 1) on the sequence
     *      attgc ==> at-tgc and returns 2
     *  - addGapsAfterIndice(1, 2) on the sequence
     *      at--tgc ==> at--t-gc and returns 5
     *  - addGapsAfterIndice(3, 5) on the sequence
     *      at--tgc ==> at--tgc--- and returns 7
     *  - addGapsAfterIndice(3, 0) on the sequence
     *      at--tgc ==> a---t--tgc and returns 1
     *  - addGapsAfterIndice(3, -5) on the sequence
     *      at--tgc ==> ---at--tgc and returns 0
     * @param nb number of gaps to add
     * @param indice which the gaps must be added after.
     * @return absolute position (eg including offset) where the gaps has been
     * added.
     */
    public int addGapsAfterIndiceAndReturnPosition(int nb, int indice)
    {
        int real_pos = this.getOffset();
        int i = 0;
        if (indice < 0)
        {
            this.offset += nb;
            return (0);
        }
        while (i < indice && i < nb_gaps.length)
            real_pos += nb_gaps[i++] + 1;
        if (i == nb_gaps.length)
        {
            this.nb_gaps[nb_gaps.length - 1] += nb;
            return (real_pos);
        }
        else if (i == indice)
            nb_gaps[i] += nb;
        return (real_pos + 1);
    }

    /**
     * Add [nb] gaps after the nucleotide at the position [pos] in the initial
     * sequence and returns the absolute position where gaps were inserted.
     * If [pos] is negative, it inserts in the offset. To add gaps at the end,
     * use [pos] greater than or equal to nb_gaps.length - 1.
     * Examples:
     *  - addGapsAfterIndice(1, 1) on the sequence
     *      attgc ==> at-tgc and returns 2
     *  - addGapsAfterIndice(1, 2) on the sequence
     *      at--tgc ==> at--t-gc and returns 5
     *  - addGapsAfterIndice(3, 5) on the sequence
     *      at--tgc ==> at--tgc--- and returns 7
     *  - addGapsAfterIndice(3, 0) on the sequence
     *      at--tgc ==> a---t--tgc and returns 1
     *  - addGapsAfterIndice(3, -5) on the sequence
     *      at--tgc ==> ---at--tgc and returns 0
     * @param nb number of gaps to add
     * @param indice which the gaps must be added after.
     * @return absolute position (eg including offset) where the gaps has been
     * added.
     */
    public int addGapsAfterIndiceEndAndReturnPosition(int nb, int indice)
    {
        int real_pos = this.getOffset();
        int i = 0;
        if (indice < 0)
        {
            this.offset += nb;
            return (this.offset - nb);
        }
        while (i < indice && i < nb_gaps.length)
            real_pos += nb_gaps[i++] + 1;
        if (i == nb_gaps.length)
        {
            //this.nb_gaps[nb_gaps.length - 1] += nb;
            return (real_pos + nb);
        }
        else if (i == indice)
            nb_gaps[i] += nb;
        return (real_pos + 1 + nb_gaps[i] - nb);
    }


    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /**
     * @return the string representation of the sequence.
     */
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
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    // GETTERS
    /**
     * Get the number of gaps before the first nucleotide eg the offset.
     * @return offset
     */
    public int getOffset()
    {
        return (this.offset);
    }

    /**
     * @return size of the sequence.
     */
    public int getSize()
    {
        int size = this.getOffset();
        for(int i = 0;i < this.nb_gaps.length;i++)
            size += 1 + this.nb_gaps[i];
        return (size);
    }

    public int getSizeWithoutEndGaps()
    {
        int size = this.getOffset();
        for(int i = 0;i < this.nb_gaps.length - 1;i++)
            size += 1 + this.nb_gaps[i];
        return (size);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    // SETTERS
    /**
     * Set the number of gaps before the first nucleotide eg the offset. Could
     * be useful to dynamically create sequences during an alignment or compute
     * the offset dynamically on a sequence.
     */
    public void setOffset(int offset)
    {
        this.offset = offset;
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Override
    /**
     * @return true if the sequences string representation are equals else
     * false.
     */
    public boolean equals(Object o)
    {
        if (o instanceof SequenceAbstract)
        {
            SequenceAbstract other = (SequenceAbstract) o;
            boolean equal = this.getOffset() == other.getOffset() && this.nb_gaps.length == other.nb_gaps.length;
            int i = 0;
            while (equal && i < this.nb_gaps.length)
                equal &= (this.nb_gaps[i] == other.nb_gaps[i] && this.initial.getLetter(i) == other.initial.getLetter(i++));
            return (equal);
        }
        return (false);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    public boolean hasSameGapsNumber(SequenceAbstract other)
    {
        int i = 0;
        boolean equal = this.nb_gaps.length == other.nb_gaps.length;
        while (equal && i < this.nb_gaps.length - 1)
            equal &= this.nb_gaps[i] == other.nb_gaps[i++];
        return (equal);
    }

    /* ---------------------------------------------------------------------- */
    /* ---------------------------------------------------------------------- */
    /*
    @Override
    public Iterator<Character> iterator()
    {
        SequenceAbstract s = this;
        Iterator<Character> it = new Iterator<Character>()
        {
            private int size = s.getSize();
            private int index = 0;
            private int gaps_i = -1;

            @Override
            public boolean hasNext()
            {
                return (index < s);
            }

            public Character next()
            {
                if (index < s.getOffset())
                    return (Character("-"));
                else if (index
                index++;
            }
        }
    }
    */
    /* ---------------------------------------------------------------------- */
}
