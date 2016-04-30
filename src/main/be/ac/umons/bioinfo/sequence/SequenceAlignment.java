package be.ac.umons.bioinfo.sequence;

/**
 * Created by aline on 31/03/16.
 */
/**
 * Represents the alignment of two DNA sequences and the score associated to this alignment.
 */
public class SequenceAlignment
{
    public final Sequence s1, s2;
    public final Sequence initial_s1, initial_s2;
    public final int score;
    public final int longestCommonSubstringLength;

    /**
     * FIXME
     * Deprecated ? Use instead the next one ?
     * @param s1 One of the aligned sequences. May contain gaps.
     * @param s2 The other aligned sequence. May contain gaps.
     * @param score The score associated to this alignment.
     * @param longestCommonSubstringLength Length of the longuest common substring
     */
    public SequenceAlignment(Sequence s1, Sequence s2, int score, int longestCommonSubstringLength)
    {
        this.s1 = s1;
        this.s2 = s2;
        this.initial_s1 = s1; // FIXME: Wrong. Initiate at this value but the constructor will be removed.
        this.initial_s2 = s2; // FIXME: Wrong. Initiate at this value but the constructor will be removed.
        this.score = score;
        this.longestCommonSubstringLength = longestCommonSubstringLength;
    }

    public SequenceAlignment(Sequence s1, Sequence s2, Sequence initial_s1, Sequence initial_s2, int score, int longestCommonSubstringLength)
    {
        this.s1 = s1;
        this.s2 = s2;
        this.initial_s1 = initial_s1;
        this.initial_s2 = initial_s2;
        this.score = score;
        this.longestCommonSubstringLength = longestCommonSubstringLength;
    }

    @Override
    public String toString()
    {
        return "SequenceAlignment(" + s1 + " , " + s2 + " , " + score + ")";
    }

    @Override
    public boolean equals(Object other)
    {
        if(! (other instanceof SequenceAlignment)) return false;

        SequenceAlignment that = (SequenceAlignment)other;

        return (this.s1.equals(that.s1) && this.s2.equals(that.s2) && (this.score == that.score));
    }

    @Override
    public int hashCode()
    {
        return ((this.s1.hashCode() + this.s2.hashCode()) % Integer.MAX_VALUE) + this.score;
    }

    /**
     * @return true if one of the aligned sequence is a subsequence of the other aligned sequence,
     * false otherwise.
     */
    public boolean inside()
    {
        return this.s1.isGapBounded() || this.s2.isGapBounded();
    }
}
