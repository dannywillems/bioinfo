package be.ac.umons.bioinfo;

/**
 * Created by aline on 31/03/16.
 */
import java.io.IOException;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import be.ac.umons.bioinfo.sequence.*;

public class Main
{
    public static void main(String[] args) throws IOException
    {
        //test();
        //cible(1, true, false, false, true);
        //cible(2, true, true, false, true);
        //cible(4, false, false, true, false);
        //cible(5, true, false, false, true);
        compute(jar(args));
    }

    public static void test()
    {
        Greedy greed = new Greedy();

        /*
        Sequence f = new Sequence("catagtc");
        Sequence g = new Sequence("taactat");
        Sequence h = new Sequence("agactatcc");

        List<Sequence> list = new ArrayList<Sequence>();

        list.add(f);
        list.add(g);
        list.add(h);
        */

        /*
        Sequence f = new Sequence("cacgt");
        Sequence g = new Sequence("acgt");
        Sequence h = new Sequence("actacg");
        Sequence i = new Sequence("gtact");
        Sequence j = new Sequence("actga");
        Sequence k = new Sequence("ctga");
        */

        Sequence f = new Sequence("actttacg");
        Sequence i = new Sequence("ggaatctgcgagtta");
        Sequence j = new Sequence("tact");
        Sequence h = new Sequence("ttgcg");
        Sequence g = new Sequence("ttgcacgat");
        Sequence k = new Sequence("gaccgat");

        List<Sequence> list = new ArrayList<Sequence>();

        list.add(j);
        list.add(f);
        list.add(k);
        list.add(h);
        list.add(g);
        list.add(i);

        /*
        Sequence f = new Sequence("catagtc");
        Sequence g = new Sequence("taactat");
        Sequence _g = g.complement();
        Sequence h = new Sequence("agactatcc");
        Sequence _h = h.complement();

        List<Sequence> list = new ArrayList<>();
        list.add(f);
        list.add(_h);
        list.add(_g);
        */

        List<SequenceAlignment> result = greed.greedy(list, 1, -1, -2);

        /*
        System.out.println("Je veux aligner :");
        System.out.println(f);
        System.out.println(h);
        System.out.println(g);
        System.out.println("------->");
        */

        ConsensusAbstract c = new ConsensusAbstract(result);
        c.computeAlignment();

        System.out.println("Les séquences alignées:\n");
        ArrayList<SequenceAbstract> alignment = c.getAlignment();
        for(int counter = 0;counter < alignment.size();counter++)
            System.out.println(alignment.get(counter));
        System.out.println("#########################################\n");

        System.out.println("Le consensus final:\n");
        SequenceAbstract consensus_final = c.build(false);
        System.out.println(consensus_final.toString());
        System.out.println("#########################################\n");
    }

    public static void cible(int num, boolean show_consensus, boolean show_alignment, boolean show_greedy, boolean save)
    {
        try
        {
            List<Sequence> list = FastaReader.readFromFile(new File("../../res/collections/cible" + num + "/collection" + num + ".fasta"));

            System.out.print("Greedy algorithm... ");
            Greedy g = new Greedy();
            List<SequenceAlignment> result = g.greedy(list, 1, -1, -2);
            System.out.println("Done");

            ConsensusAbstract c = new ConsensusAbstract(result);
            c.updateOffset();

            if (show_greedy)
            {
                for(int i = 0;i < c.getHamiltonianPath().size();i++)
                {
                    System.out.println(c.getHamiltonianPath().get(i).s);
                    System.out.println(c.getHamiltonianPath().get(i).t);
                    System.out.println("###");
                }
                System.out.println("#########################################\n");
            }

            System.out.print("Alignement... ");
            c.computeAlignment();
            System.out.println("Done");

            if (show_alignment)
            {
                System.out.println("Les séquences alignées:\n");
                ArrayList<SequenceAbstract> alignment = c.getAlignment();
                for(int counter = 0;counter < alignment.size();counter++)
                    System.out.println(alignment.get(counter));
                System.out.println("#########################################\n");
            }

            System.out.print("Construction du consensus... ");
            SequenceAbstract consensus_final = c.build(false);
            Sequence n_s = new Sequence(consensus_final.toString());
            n_s = n_s.complement();
            System.out.println("Done");

            if (show_consensus)
            {
                System.out.println("Le consensus final:\n");
                System.out.println(n_s);
                System.out.println("#########################################\n");
            }

            if (save)
            {
                String file_name = "Cible_calcul" + num + ".fasta";
                String directory = "../../res/results/cible" + num + "/";
                System.out.print("Ecriture du consensus dans le fichier " + file_name + "... ");
                FastaWriter.write("Cible" + num, n_s, new File(directory + file_name), 80);
                System.out.println("Done");
            }
        }
        catch (IOException e)
        {
            System.out.println(e.toString());
        }
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
            Greedy g = new Greedy();
            List<SequenceAlignment> result = g.greedy(list, 1, -1, -2);

            ConsensusAbstract c = new ConsensusAbstract(result);
            c.computeAlignment();
            SequenceAbstract consensus_final = c.build(false);
            FastaWriter.write("Normal", consensus_final, new File(f[1]), 80);
            FastaWriter.write("Complementary inverse", consensus_final.complement(), new File(f[2]), 80);
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
