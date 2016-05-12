package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.greedy.*;

/**
 * Created by aline on 31/03/16.
 */

import be.ac.umons.bioinfo.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List; /** Represents a DNA sequence
 */


public class Sequence
{
    public byte[] content;

    /**
     * @param seq A DNA sequence to encode
     */

    public Sequence(String seq)
    {
        this.content = letter2Base(seq);
    }

    public Sequence(byte[] seq)
    {
        this.content = seq;
    }

    public Sequence(Byte[] seq)
    {
        byte[] content = new byte[seq.length];

        for(int i=0; i<seq.length ; i++)
            content[i] = seq[i];

        this.content = content;
    }

    public void setContent(String s)
    {
        this.content = letter2Base(s);
    }

    public static byte[] letter2Base(String text)
    {
        byte[] ret = new byte[text.length()];

        for(int i = 0; i<text.length() ; i++)
            ret[i] = Nucleotide.letter2Base(text.charAt(i));

        return ret;
    }

    public static String base2letter(byte[] bases)
    {
        StringBuilder builder = new StringBuilder();

        for(byte c : bases)
            builder.append(Nucleotide.base2letter(c));

        return builder.toString();
    }

    public char getLetter(int i)
    {
        return (Nucleotide.base2letter(this.content[i]));
    }

    public Pair<Byte> getPairBase(Sequence that, int i)
    {
        return new Pair<>(this.content[i], that.content[i]);
    }

    public byte getBaseAsByte(int i)
    {
        return (this.content[i]);
    }

    @Override
    public String toString()
    {
        StringBuilder builder = new StringBuilder();

        for(byte c : this.content)
            builder.append(Nucleotide.base2letter(c));

        return builder.toString();
    }

    /**
     * FIXME: strange behavior when changing bytes... Cast doesn't work ? How
     * does the cast is done?
     */
    public void addByteAtPos(byte b, int i)
    {
        StringBuilder s = new StringBuilder();
        int j = 0;
        while (j < this.content.length)
        {
            if (i == j)
            {
                s.append(Nucleotide.base2letter(b));
                i = -1;
            }
            else
            {
                s.append(this.getLetter(j));
                j++;
            }
        }
        this.content = this.letter2Base(s.toString());
    }

    public void addGapAtPos(int i)
    {
        this.addByteAtPos((byte) Nucleotide.GAP, i);
    }

    /**
     * Get the size of the sequence which is the number of bases.
     * @return The number of bases in this DNA sequence.
     */
    public int getSize()
    {
        return this.content.length;
    }

    /**
     * Reverses the order of the bases contained in this sequence,
     * and replaces their values by their binded base.
     *
     * Ex: AATG becomes CATT
     *
     * @return The DNA sequence being the complement to this sequence.
     */
    public Sequence complement()
    {
        byte[] comp = new byte[this.content.length];

        for(int i = 0; i<this.content.length ; i++)
            comp[i] = Nucleotide.complement(this.content[this.content.length - 1 - i]);

        return new Sequence(comp);
    }

    /**
     * Computes the score associated an alignment of this sequence and an other one.
     * This score is based on the fact that matches, mismatches, and gaps have different values.
     * @param that the other sequence to consider
     * @param match the score associated to a match
     * @param mismatch the score associated to a mismatch
     * @param gap the score associated to a gap
     * @return the score associated to the alignment of this sequence and an other one, according to the score rules.
     */
    public int alignmentScore(Sequence that, int match, int mismatch, int gap)
    {
        int minLength = Math.min(this.getSize(), that.getSize());
        int maxLength = Math.max(this.getSize(), that.getSize());

        // We fill with gap if sequences don't have the same size.
        int score = (maxLength - minLength) * gap;

        for(int i = 0; i<minLength ; i++)
        {
            if (this.content[i] == Nucleotide.GAP || that.content[i] == Nucleotide.GAP)
                score += gap;
            else if (this.content[i] == that.content[i])
                score += match;
            else
                score += mismatch;
        }

        return score;
    }

    /**
     * Two sequences are equals if they
     * - have the same size
     * - contains the same nucleotides.
     */
    @Override
    public boolean equals(Object that)
    {
        if(that instanceof Sequence)
        {
            Sequence s = (Sequence) that;
            if(this.getSize() != s.getSize()) return false;

            for(int i = 0; i<this.getSize() ; i++)
                if(this.content[i] != s.content[i]) return false;

            return true;
        }
        else return false;
    }

    @Override
    public int hashCode()
    {
        return Arrays.hashCode(content);
    }

    /**
     * Cost function for a match or a mismatch between two bases.
     * @param a a base
     * @param b an other base
     * @param match cost of a match between a and b
     * @param mismatch cost of a mismatch between a and b
     * @return the cost of an alignment between a and b
     */
    public static int p(byte a, byte b, int match, int mismatch)
    {
        if(a == b) return match;
        else return mismatch;
    }

    /**
     * Generates a list of arcs associated to a pair of sequences.
     * @param that an other sequence.
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return the arcs associated to a pair of sequences.
     */
    public List<Arc> arcGenerator(Sequence that, int match, int mismatch, int gap)
    {
        Sequence compF = this.complement();
        Sequence compG = that.complement();

        List<SequenceAlignment> FG = semiGlobalAlignment(this, that, match, mismatch, gap);
        List<SequenceAlignment> CompFG = semiGlobalAlignment(compF, that, match, mismatch, gap);
        List<SequenceAlignment> FCompG = semiGlobalAlignment(this, compG, match, mismatch, gap);
        List<SequenceAlignment> CompFCompG = semiGlobalAlignment(compF, compG, match, mismatch, gap);

        List<Arc> arcs = new ArrayList<>();

        assert(FG.size() == 2);
        assert(CompFG.size() == 2);
        assert (FCompG.size() == 2);
        assert (CompFCompG.size() == 2);

        arcs.add(new Arc(
                this,
                false,
                that,
                false,
                FG.get(0).score,
                match,
                mismatch,
                gap, FG.get(0).inside()));

        arcs.add(new Arc(
                that,
                false,
                this,
                false,
                FG.get(1).score,
                match,
                mismatch,
                gap,
                FG.get(1).inside()));

        arcs.add(new Arc(this,
                true,
                that,
                false,
                CompFG.get(0).score,
                match,
                mismatch,
                gap,
                CompFG.get(0).inside()));

        arcs.add(new Arc(
                that,
                false,
                this,
                true,
                CompFG.get(1).score,
                match,
                mismatch,
                gap,
                CompFG.get(1).inside()));

        arcs.add(new Arc(
                this,
                false,
                that,
                true,
                FCompG.get(0).score,
                match,
                mismatch,
                gap,
                FCompG.get(0).inside()));

        arcs.add(new Arc(
                that,
                true,
                this,
                false,
                FCompG.get(1).score,
                match,
                mismatch,
                gap,
                FCompG.get(1).inside()));

        arcs.add(new Arc(
                this,
                true,
                that,
                true,
                CompFCompG.get(0).score,
                match,
                mismatch,
                gap,
                CompFCompG.get(0).inside()));

        arcs.add(new Arc(
                that,
                true,
                this,
                true,
                CompFCompG.get(1).score,
                match,
                mismatch,
                gap,
                CompFCompG.get(1).inside()));

        return arcs;
    }

    /**
     * Given two sequences s1 and s2, cost of match, mismatch and gap,
     * computes the matrix of similarities  between s1 and s2
     * @param s1 a sequence
     * @param s2 an other sequence
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return the matrix of simiraties between s1 and s2
     */

    public static int[][] computeSimMat(Sequence s1, Sequence s2, int match, int mismatch, int gap)
    {
        final int m = s1.getSize();
        final int n = s2.getSize();

        int[][] a = new int[m+1][n+1];

        // Initialisation

        for(int i = 0; i<=m ; i++)
            a[i][0] = 0;

        for(int j = 0; j<=n ; j++)
            a[0][j] = 0;

        // Filling

        for(int i=1 ; i<= m ; i++)
        {
            for(int j=1 ; j<= n ; j++)
            {
                final int x = a[i-1][j] + gap;
                final int y = a[i-1][j-1] + p(s1.content[i-1], s2.content[j-1], match, mismatch);
                final int z = a[i][j-1] + gap;

                a[i][j] = Math.max(Math.max(x, y), z);
            }
        }
        return a;
    }

    /**
     * Computes the semi-global aligment of this sequence and an other one, considering
     * specific costs for a match, a mismatch, and a gap.
     * @param s1 a sequence
     * @param s2 an other sequence
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return the result of the alignment of this sequence with the other one.
     */
    public static List<SequenceAlignment> semiGlobalAlignment(Sequence s1, Sequence s2, int match, int mismatch, int gap)
    {

        int a[][] = computeSimMat(s1,s2, match, mismatch, gap);

        SequenceAlignment sBefore = backtrack(a, s1, s2, true, match, mismatch, gap);
        SequenceAlignment tBefore = backtrack(a, s1, s2, false, match, mismatch, gap);

        List<SequenceAlignment> ret = new ArrayList<>();

        ret.add(sBefore);
        ret.add(tBefore);

        return ret;
    }



    /**
     * Creates a sequence alignment by bactracking the similarity matrix from its bottom.
     * @param a the similarity matrix.
     * @param s1 a sequence.
     * @param s2 an other sequence.
     *
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return If the two considered sequence are not included in each other, the alignment of these sequences.
     * Otherwise, returns nothing.
     */
    private static SequenceAlignment backtrack(int[][] a,
                                               Sequence s1,
                                               Sequence s2,
                                               boolean bottom,
                                               int match,
                                               int mismatch,
                                               int gap)
    {
        int m = a.length - 1;
        int n = a[0].length - 1;

        int iPos = m;
        int iMax = Integer.MIN_VALUE;

        int jPos = n;
        int jMax = Integer.MIN_VALUE;

        if(bottom)
        {
            for(int j = n; j >= 0 ; j--)
            {
                if(a[m][j] > jMax)
                {
                    jMax = a[m][j];
                    jPos = j;
                }
            }
        }
        else
        {
            for(int i = m; i >= 0 ; i--)
            {
                if(a[i][n] > iMax)
                {
                    iMax = a[i][n];
                    iPos = i;
                }
            }
        }

        int x = jPos;
        int y = iPos;

        List<Byte> aligned_s = new ArrayList<Byte>();
        List<Byte> aligned_t = new ArrayList<Byte>();


        if(x == n)  // Right border
        {
            for(int i=m ; i>y ; i--)
            {
                aligned_t.add((byte) Nucleotide.GAP);
                aligned_s.add(s1.content[i-1]);
            }
        }
        else
        {
            for(int j=n ; j>x ; j--)
            {
                aligned_s.add((byte) Nucleotide.GAP);
                aligned_t.add(s2.content[j-1]);
            }
        }

        while(x > 0 && y > 0)
        {
            final int b = a[y][x - 1] + gap;
            final int c = a[y - 1][x - 1] + p(s1.content[y - 1], s2.content[x - 1], match, mismatch);

            if (c == a[y][x])
            {
                aligned_s.add(s1.content[y - 1]);
                aligned_t.add(s2.content[x - 1]);

                x -= 1;
                y -= 1;
            } else
            {
                if (b == a[y][x])
                {
                    aligned_s.add((byte) Nucleotide.GAP);
                    aligned_t.add((byte) s2.content[x - 1]);

                    x -= 1;
                } else
                {
                    aligned_s.add((byte) s1.content[y - 1]);
                    aligned_t.add((byte) Nucleotide.GAP);

                    y -= 1;
                }
            }
        }

        if (x > 0)
        {
            for(int j=x ; j>0 ; j--)
            {
                aligned_s.add((byte) Nucleotide.GAP);
                aligned_t.add(s2.content[j-1]);
            }
        } else
        {
            for(int i=y ; i>0 ; i--)
            {
                aligned_t.add((byte) Nucleotide.GAP);
                aligned_s.add(s1.content[i-1]);
            }
        }

        Collections.reverse(aligned_s);
        Collections.reverse(aligned_t);

        Byte[] sArray = new Byte[aligned_s.size()];
        aligned_s.toArray(sArray);
        Sequence alignedS = new Sequence(sArray);

        Byte[] tArray = new Byte[aligned_t.size()];
        aligned_t.toArray(tArray);
        Sequence alignedT = new Sequence(tArray);


        if(bottom)
            return new SequenceAlignment(alignedS, alignedT, jMax);
        else
             return new SequenceAlignment(alignedT, alignedS, iMax);
    }

    /**
    * Get the indice of the last nucleotide in the sequence. For example, if
    * the sequence is
    * --acg--
    * the function will return 4.
    * Useful for the alignment.
    * @return the indice of the last nucleotide in the sequence
    */
    public int getPosLastNucleotide()
    {
        int s = this.getSize() - 1;
        int current = s;
        for(int i = s;i >= 0;i--)
        {
            if (this.content[i] == (byte) Nucleotide.GAP)
                current--;
            else
                return (current);
        }
        return (current);
    }

    /**
     * Get the indice of the first nucleotide in the sequence. For example, if
     * the sequence is
     * --acg--
     * the function will return 2.
     * Useful for the alignment.
     * @return the indice of the first nucleotide in the sequence
     */

    public int getPosFirstNucleotide()
    {
        int current = 0;
        for(int i = 0;i < this.getSize();i++)
        {
            if (this.content[i] == (byte) Nucleotide.GAP)
                current++;
            else
                return (current);
        }
        return (current);
    }

    public int nbGapEnd()
    {
        return (this.getSize() - 1 - this.getPosLastNucleotide());
    }

    public int nbGapBegin()
    {
        return (this.getPosFirstNucleotide());
    }

    /**
     * @return true if this sequence is bounded by gaps; false otherwise.
     */
    public boolean isGapBounded()
    {
        return this.content[0] == Nucleotide.GAP && this.content[this.content.length-1] == Nucleotide.GAP;
    }

    /**
     *
     * @return true if there is no gap at the two extremities of the sequence; false otherwise
     */
    public boolean noExternalGap(){return this.content[0] != Nucleotide.GAP &&
            this.content[this.content.length-1] != Nucleotide.GAP;}

    public Sequence rebuildAddingGaps(int[] gaps)
    {
        StringBuilder s_final = new StringBuilder();
        int gaps_i = 0;
        int i = 0;
        while (i < this.getSize())
        {
            if (gaps_i < gaps.length && gaps[gaps_i] == i)
            {
                s_final.append(Nucleotide.base2letter((byte) Nucleotide.GAP));
                gaps_i++;
            }
            else
            {
                s_final.append(this.getLetter(i));
                i++;
            }
        }
        return (new Sequence(s_final.toString()));
    }
}
