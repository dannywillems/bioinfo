package be.ac.umons.sequence;

/**
 * Created by aline on 31/03/16.
 */
import java.io.IOException;

public class Main
{

    public static void main(String[] args) throws IOException
    {
        Sequence s = new Sequence("cagcacttggattctcgg");
        Sequence t = new Sequence("cagcgtgg");

        s.semiGlobalAlignment(t, 1, -1, -2);
    }
}
