package be.ac.umons.bioinfo;

/**
 * Created by aline on 31/03/16.
 */
import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import be.ac.umons.bioinfo.consensus.*;
import be.ac.umons.bioinfo.greedy.*;
import be.ac.umons.bioinfo.sequence.*;
import be.ac.umons.bioinfo.fasta.FastaReader;
import be.ac.umons.bioinfo.fasta.FastaWriter;

public class Main
{
    public static void main(String[] args) throws IOException
    {
        compute(jar(args));
    }

    public static String[] jar(String[] args)
    {
        String[] f = new String[3];

        if (args.length == 5)
        {
            f[0] = args[0];
            if (args[1].equals("--out"))
            {
                f[1] = args[2];
                if (args[3].equals("--out-ic"))
                {
                    f[2]= args[4];
                    return (f);
                }
            }
        }
        showHelp();
        return (null);
    }

    public static void compute(String[] f) throws IOException
    {
        try
        {
            List<Sequence> list = FastaReader.readFromFile(new File(f[0]));
            List<SequenceAlignment> result = Greedy.greedy(list, 1, -1, -2);

            ConsensusAbstract c = new ConsensusAbstract(result);
            c.computeAlignment();
            SequenceAbstract consensus_final = c.build(false);
            FastaWriter.write(" Groupe 9 Longueur " + consensus_final.getSize(), consensus_final, new File(f[1]), 80);
            FastaWriter.write(" Groupe 9 Longueur " + consensus_final.getSize(), consensus_final.complement(), new File(f[2]), 80);
        }
        catch (IOException e)
        {
            System.out.println(e);
        }
    }

    public static void showHelp()
    {
        System.out.println("java -jar FragmentAssembler.jar <fichier.fasta> --out <sortie.fasta> --out-ic <sortie-ic.fasta>");
    }
}
