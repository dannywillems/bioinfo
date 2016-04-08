package be.ac.umons.bioinfo;

/**
 * Created by aline on 31/03/16.
 */
import java.io.IOException;
import be.ac.umons.bioinfo.sequence.*;

public class Main
{
    public static void main(String[] args) throws IOException
    {
        Sequence s = new Sequence("attagaccatgcggc");
        Sequence t = new Sequence("atcggcattcagt");

        Sequence.semiGlobalAlignment(s, t, 1, -1, -2);
    }
}
