package be.ac.umons.bioinfo.sequence;

import java.lang.Integer;
import java.util.ArrayList;
import java.util.Iterator;

public class SequenceAbstract implements Iterable<Character>
{
    public Sequence initial;
    private int[] nb_gaps;
    private int offset;
    private int size;
    /* ---------------------------------------------------------------------- */
    // CONSTRUCTORS
    /**
     * Build a SequenceAbstract object from an initial sequence and the array of
     * gaps.
     *
     * Initial sequence is useless, if no mutation, the initial sequence is the
     * aligned sequence where gaps has been removed
     *
     * Complexity: O(1)
     */
    public SequenceAbstract(Sequence initial, int[] nb_gaps)
    {
        this.initial = initial;
        this.nb_gaps = nb_gaps;
        this.offset = 0;
        this.updateSize();
    }

    /**
     * Build a SequenceAbstract object from an initial sequence, the array of
     * gaps and the offset.
     *
     * Initial sequence is useless, if no mutation, the initial sequence is the
     * aligned sequence where gaps has been removed
     *
     * Complexity: O(1)
     */
    public SequenceAbstract(Sequence initial, int[] nb_gaps, int offset)
    {
        this.initial = initial;
        this.nb_gaps = nb_gaps;
        this.offset = offset;
        this.updateSize();
    }

    /**
     * Build SequenceAbtract object from an aligned sequence and his initial
     * sequence.
     *
     * Initial sequence is useless, if no mutation, the initial sequence is the
     * aligned sequence where gaps has been removed. Only used to defined the
     * size of the gaps array but can be replaced by an arraylist.
     *
     * Complexity: O(n) where n is the size of the aligned sequence.
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
        this.updateSize();
    }

    /**
     * Create a new abstract sequence from an aligned sequence (or not) eg
     * derive the initial sequence and compute the gaps array.
     *
     * Complexity: O(n) where n is the size of the
     * aligned sequence.
     */
    public SequenceAbstract(Sequence s)
    {
        this.offset = 0;
        while (s.getBaseAsByte(this.offset) == (byte) Nucleotide.GAP)
            this.offset++;

        StringBuilder str = new StringBuilder(s.getSize());
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
        this.updateSize();
    }

    /**
     * Create a new abstract sequence from a string.
     *
     * Complexity: O(n * k^2) where n is the size of the aligned sequence and k
     * the number of nucleotides in the initial sequence. k^2 is due to
     * ArrayList.add method. charAt is in contant time.
     */
    public SequenceAbstract(String s)
    {
        this.offset = 0;
        while (Nucleotide.letter2Base(s.charAt(this.offset)) == (byte) Nucleotide.GAP)
            this.offset++;

        StringBuilder str = new StringBuilder(s.length());
        int j;
        int i = this.offset;
        ArrayList<Integer> gaps = new ArrayList<Integer>(s.length());
        while (i < s.length())
        {
            str.append(s.charAt(i));
            j = 0;
            while (i + j + 1 < s.length() && Nucleotide.letter2Base(s.charAt(i + j + 1)) == (byte) Nucleotide.GAP)
                j++;
            gaps.add(j);
            i = i + j + 1;
        }

        this.nb_gaps = new int[gaps.size()];
        for(int k = 0;k < gaps.size();k++)
            this.nb_gaps[k] = gaps.get(k);

        this.initial = new Sequence(str.toString());
        this.updateSize();
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /**
     * Compute and return the complement inverse of the sequence. Faster than
     * complement in Sequence if lots of gaps.
     * @return the complement inverse of the sequence.
     *
     * Complexity: O(2n) = O(n) where n is the length of the initial sequence.
     */
    public SequenceAbstract complement()
    {
        Sequence init = this.initial.complement();
        int[] gaps = new int[this.nb_gaps.length];

        gaps[this.nb_gaps.length - 1] = this.getOffset();

        for(int i = 0;i < this.nb_gaps.length - 2;i++)
            gaps[i] = this.nb_gaps[this.nb_gaps.length - i - 2];

        return (new SequenceAbstract(init, gaps, this.nb_gaps[this.nb_gaps.length - 1]));
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
     *
     * Complexity: O(n) where n is the length of the initial sequence.
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
        this.updateSize(nb);
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
     *
     * Complexity: O(n) where n is the length of the initial sequence.
     */
    public int addGapsAndReturnIndice(int nb, int pos)
    {
        int i = 0;
        int real_pos = this.offset;
        // if we need to add gaps before the beginning of the sequence
        if (pos <= real_pos)
        {
            this.updateSize(nb);
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
                    this.updateSize(nb);
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
                    this.updateSize(nb);
                    this.nb_gaps[i] += nb;
                    return (i);
                }
            }
            this.updateSize(nb);
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
     *
     * Complexity: O(1)
     */
    public void addGapsAfterIndice(int nb, int indice)
    {
        if (indice < 0)
            this.offset += nb;
        else if (indice < nb_gaps.length - 1)
            this.nb_gaps[indice] += nb;
        else
            this.nb_gaps[nb_gaps.length - 1] += nb;
        this.updateSize(nb);
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
     *
     * Complexity: O(n) where n is the length of the initial sequence. This is
     * due to the necessity to return the position.
     */
    public int addGapsAfterIndiceAndReturnPosition(int nb, int indice)
    {
        int real_pos = this.getOffset();
        int i = 0;
        if (indice < 0)
        {
            this.updateSize(nb);
            this.offset += nb;
            return (0);
        }
        while (i < indice && i < nb_gaps.length)
            real_pos += nb_gaps[i++] + 1;
        if (i == nb_gaps.length)
        {
            this.updateSize(nb);
            this.nb_gaps[nb_gaps.length - 1] += nb;
            return (real_pos);
        }
        else if (i == indice)
            nb_gaps[i] += nb;
        this.updateSize(nb);
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
     *
     * Complexity: O(n) where n is the length of the initial sequence. This is
     * due to the necessity to return the position.
     */
    public int addGapsAfterIndiceEndAndReturnPosition(int nb, int indice)
    {
        int real_pos = this.getOffset();
        int i = 0;
        if (indice < 0)
        {
            this.updateSize(nb);
            this.offset += nb;
            return (this.offset - nb);
        }
        while (i < indice && i < nb_gaps.length)
            real_pos += nb_gaps[i++] + 1;
        if (i == nb_gaps.length)
        {
            this.updateSize(nb);
            return (real_pos + nb);
        }
        else if (i == indice)
            nb_gaps[i] += nb;

        this.updateSize(nb);
        return (real_pos + 1 + nb_gaps[i] - nb);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /**
     * @return the string representation of the sequence.
     *
     * Complexity: O(m) where m is the length of the aligned sequence.
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
        return (this.size);
    }

    public int getSizeWithoutEndGaps()
    {
        int size = this.getOffset();
        for(int i = 0;i < this.nb_gaps.length - 1;i++)
            size += 1 + this.nb_gaps[i];
        return (size);
    }

    public int[] getNbGaps()
    {
        return (this.nb_gaps);
    }

    public int getNbGapsLength()
    {
        return (this.nb_gaps.length);
    }

    public int getNbGaps(int i)
    {
        return (this.nb_gaps[i]);
    }

    private void updateSize()
    {
        this.size = this.getOffset();
        for(int i = 0;i < this.initial.getSize();i++)
            this.size += 1 + this.nb_gaps[i];
    }

    private void updateSize(int nb)
    {
        this.size += nb;
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
        this.updateSize();
    }

    public void addNbGaps(int i, int nb)
    {
        this.nb_gaps[i] += nb;
        this.updateSize(nb);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    @Override
    /**
     * @return true if the sequences string representation are equals else
     * false.
     *
     * Complexity: O(n) where n is the length of the initial sequences. O(1) if
     * the the sequences doesn't have the same number of nucleotide in their
     * initial sequences.
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
    /**
     * @return true if the sequence has the same number of gaps between its
     * nucleotides than the [other] sequence abstract passed in parameter.
     *
     * Complexity: O(n) where n is the length of the initial sequence. O(1) if
     * the the sequences doesn't have the same number of nucleotide in their
     * initial sequences.
     */
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
    @Override
    public Iterator<Character> iterator()
    {
        SequenceAbstract s = this;
        Iterator<Character> it = new Iterator<Character>()
        {
            private int real_pos = s.getOffset();
            private int size = s.getSize();
            private int index = 0;
            private int gaps_i = 0;

            @Override
            public boolean hasNext()
            {
                return (index < size);
            }

            public Character next()
            {
                Character c;
                if (index < s.getOffset())
                    c = new Character('-');
                else if (index == s.getOffset())
                    c = new Character(s.initial.getLetter(gaps_i));
                else if (gaps_i == s.nb_gaps.length)
                    c = new Character('-');
                else if (index == real_pos + s.nb_gaps[gaps_i] + 1)
                {
                    real_pos += s.nb_gaps[gaps_i] + 1;
                    gaps_i++;
                    c = new Character(s.initial.getLetter(gaps_i));
                }
                else
                    c = new Character('-');
                index++;
                return (c);
            }

            public void remove()
            {
                throw new UnsupportedOperationException();
            }
        };
        return (it);
    }
    /* ---------------------------------------------------------------------- */
}
