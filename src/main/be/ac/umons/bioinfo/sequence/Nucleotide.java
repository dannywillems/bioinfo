package be.ac.umons.bioinfo.sequence;

public class Nucleotide
{
    public static final byte GAP = 4;
    public static final byte C = 0;
    public static final byte G = 1;
    public static final byte T = 2;
    public static final byte A = 3;

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

    public static byte complement(byte base)
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
}
