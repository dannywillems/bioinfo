package be.ac.umons.bioinfo.sequence;

/**
 * Created by aline on 31/03/16.
 */
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Represents a DNA sequence
 */
public class Sequence
{
    public static final int GAP = 4;
    public static final int C = 0;
    public static final int G = 1;
    public static final int T = 2;
    public static final int A = 3;
    private byte[] content;



    /**
     * @param seq A DNA sequence to encode
     */
    public Sequence(String seq)
    {
        this.content = letter2Base(seq);
    }

    private Sequence(byte[] seq)
    {
        this.content = seq;
    }

    private Sequence(Byte[] seq)
    {
        byte[] content = new byte[seq.length];

        for(int i=0; i<seq.length ; i++)
            content[i] = seq[i];

        this.content = content;
    }

    public static byte letter2Base(char c)
    {

        switch(java.lang.Character.toLowerCase(c))
        {
            case 'c' : return C;
            case 'g' : return G;
            case 't' : return T;
            case 'a' : return A;
            case '-' : return GAP;
            default : return 5;
        }
    }

    public static byte[] letter2Base(String text)
    {
        byte[] ret = new byte[text.length()];

        for(int i = 0; i<text.length() ; i++)
            ret[i] = letter2Base(text.charAt(i));

        return ret;
    }

    public static char base2letter(byte c)
    {
        switch(c)
        {
            case C: return 'c';
            case G : return 'g';
            case T : return 't';
            case A : return 'a';
            case GAP: return '-';
            default : return '?';
        }
    }

    public static String base2letter(byte[] bases)
    {
        StringBuilder builder = new StringBuilder();

        for(byte c : bases)
            builder.append(base2letter(c));

        return builder.toString();
    }

    @Override
    public String toString()
    {
        StringBuilder builder = new StringBuilder();

        for(byte c : this.content)
            builder.append(base2letter(c));

        return builder.toString();
    }

    /**
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
            comp[i] = complement(this.content[this.content.length - 1 - i]);

        return new Sequence(comp);
    }

    private static final byte complement(byte base)
    {
        switch(base)
        {
            case C: return G;
            case G : return C;
            case T : return A;
            case A : return T;
            default : return -1;
        }
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

        int score = (maxLength - minLength) * gap;

        for(int i = 0; i<minLength ; i++)
        {
            if(this.content[i] == GAP || that.content[i] == GAP) score += gap;
            else if(this.content[i] == that.content[i]) score += match;
            else score += mismatch;
        }

        return score;
    }

    @Override
    public boolean equals(Object that)
    {
        if(that instanceof Sequence)
        {
            Sequence s = (Sequence) that;
            if(this.getSize() != s.getSize()) return false;

            for(int i = C; i<this.getSize() ; i++)
                if(this.content[i] != s.content[i]) return false;

            return true;
        }
        else return false;
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
     * Computes the semi-global aligment of this sequence and an other one, considering
     * specific costs for a match, a mismatch, and a gap.
     * @param that an other sequence
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @param h the cost of a first gap
     * @return the similarity score of the sequences alignment
     */
    public int globalAlignment(Sequence that, int match, int mismatch, int gap, int h)
    {
        final int m = this.getSize();
        final int n = that.getSize();

        int[][] a = new int[m +1][n +1];
        int[][] b = new int[m +1][n +1];
        int[][] c = new int[m +1][n +1];


        // Initialisation

        a[0][0] = 0;

        for(int i = 1; i<= m; i++)
            a[i][0] = Integer.MIN_VALUE;

        for(int j = 1; j<= n; j++)
            a[0][j] = Integer.MIN_VALUE;

        for(int i = 1; i<= m; i++)
            b[i][0] = Integer.MIN_VALUE;

        for(int j = 1; j<= n; j++)
            b[0][j] = -(h+(gap*j));

        for(int i=1 ; i<= m ; i++)
            c[i][0] = -(h+(gap*i));

        for(int j=1 ; j<= n ; j++)
            c[0][j] = Integer.MIN_VALUE;


        // Filling

        for(int i=1 ; i<=m ; i++)
        {
            for(int j=i ; j<=n ; j++)
            {
                final int xa = a[i-1][j-1];
                final int ya = b[i-1][j-1];
                final int za = c[i-1][j-1];

                a[i][j] = p(this.content[i-1], that.content[j-1], match, mismatch) + Math.max(Math.max(xa,ya),za);

                final int xb = - (h+gap) + a[i][j-1];
                final int yb = - gap + b[i][j-1];
                final int zb = - (h + gap) + c[i][j-1];

                b[i][j] = Math.max(Math.max(xb, yb), zb);

                final int xc = - (h + gap) + a[i-1][j];
                final int yc = - (h + gap) + b[i-1][j];
                final int zc = - gap + c[i-1][j];

                c[i][j] = Math.max(Math.max(xc, yc), zc);
            }
        }

        int similarity = Math.max(Math.max(a[m][n], b[m][n]), c[m][n]);

        return similarity;
    }

    /**
     * Computes the semi-global aligment of this sequence and an other one, considering
     * specific costs for a match, a mismatch, and a gap.
     * @param that an other sequence
     * @param match the cost of a match
     * @param mismatch the cost of a mismatch
     * @param gap the cost of a gap
     * @return the result of the alignment of this sequence with the other one.
     */
    public SequenceAlignment semiGlobalAlignment(Sequence that, int match, int mismatch, int gap)
    {
        final int m = this.getSize();
        final int n = that.getSize();

        int[][] a = new int[this.getSize()+1][that.getSize()+1];

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
                final int y = a[i-1][j-1] + p(this.content[i-1], that.content[j-1], match, mismatch);
                final int z = a[i][j-1] + gap;

                a[i][j] = Math.max(Math.max(x, y), z);
            }
        }

        int iMax = 0;
        int jMax = 0;

        int iPos = Integer.MIN_VALUE;
        int jPos = Integer.MIN_VALUE;

        for(int i = 0; i<=m ; i++)
        {
            if(a[i][n] > jMax)
            {
                jMax = a[i][n];
                jPos = i;
            }
        }

        for(int j = 0; j<=n ; j++)
        {
            if(a[m][j] > iMax)
            {
                iMax = a[m][j];
                iPos = j;
            }
        }



        int similarity = 0;
        int x, y;

        if(iMax > jMax)
        {
            similarity = iMax;
            x = iPos;
            y = m;
        }
        else
        {
            similarity = jMax;
            x = n;
            y = jPos;
        }

        for(int i = 0; i <= m ; i++)
        {
            for(int j = 0; j <= n ; j++)
            {
                System.out.print(a[i][j] + " \t");
            }

            System.out.printf("\n");
        }

        System.out.println("Max in " + x + " " + y + " with score " + similarity);

        List<Byte> aligned_s = new ArrayList<Byte>();
        List<Byte> aligned_t = new ArrayList<Byte>();

        if(x == n)  // Right border
        {
            for(int i=m ; i>y ; i--)
            {
                aligned_t.add((byte) GAP);
                aligned_s.add(this.content[i-1]);
            }
        }
        else
        {
            for(int j=n ; j>x ; j--)
            {
                aligned_s.add((byte) GAP);
                aligned_t.add(that.content[j-1]);
            }
        }

        while(x > 0 && y > 0)
        {
            final int b = a[y][x - 1] + gap;
            final int c = a[y - 1][x - 1] + p(this.content[y - 1], that.content[x - 1], match, mismatch);
            // final int d = a[y-1][x] + gap;

            if (c == a[y][x])
            {
                aligned_s.add((byte) this.content[y - 1]);
                aligned_t.add((byte) that.content[x - 1]);

                x -= 1;
                y -= 1;
            } else
            {
                if (b == a[y][x])
                {
                    aligned_s.add((byte) GAP);
                    aligned_t.add((byte) that.content[x - 1]);

                    x -= 1;
                } else
                {
                    aligned_s.add((byte) this.content[y - 1]);
                    aligned_t.add((byte) GAP);

                    y -= 1;
                }
            }

            System.out.println(aligned_s);
            System.out.println(aligned_t);
            System.out.println("*******");
        }

        if (x > 0)
        {
            for(int j=x ; j>0 ; j--)
            {
                aligned_s.add((byte) GAP);
                aligned_t.add(that.content[j-1]);
            }
        } else
        {
            for(int i=y ; i>0 ; i--)
            {
                aligned_t.add((byte) GAP);
                aligned_s.add(this.content[i-1]);
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

        return new SequenceAlignment(alignedS, alignedT, similarity);
    }
}
