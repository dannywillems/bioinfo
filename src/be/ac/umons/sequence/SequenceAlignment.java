package be.ac.umons.sequence;

/**
 * Created by aline on 31/03/16.
 */
/**
 * Represents the alignment of two DNA sequences.
 */
public class SequenceAlignment
{
    public final Sequence s1, s2;
    public final int score;

    /**
     * @param s1 One of the aligned sequences. May contain gaps.
     * @param s2 The other aligned sequence. May contain gaps.
     * @param score The score associated to this alignment.
     */
    public SequenceAlignment(Sequence s1, Sequence s2, int score)
    {
        this.s1 = s1;
        this.s2 = s2;
        this.score = score;
    }
}
