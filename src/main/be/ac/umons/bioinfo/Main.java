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
        test();
        //cible(1, false, false, true);
        //cible(2, false, false, true);
        //cible(3, false, false, true); // A exécuter en dernier, car long

        //test2();


    }


    public static void test2()
    {
        Sequence s1 = new Sequence("act-gtg-");
        Sequence s2 = new Sequence("--ta-tac");
        Sequence _s2 = new Sequence("t--at-ac---");
        Sequence s3 = new Sequence("-ac-tc-cgta");

        SequenceAlignment sA1 = new SequenceAlignment(s1, s2, 42, 2);
        SequenceAlignment sA2 = new SequenceAlignment(_s2,s3, 42, 2);

        List<SequenceAlignment>list = new ArrayList<>();
        list.add(sA1);
        list.add(sA2);

        ConsensusAline consensus = new ConsensusAline(list);
        Sequence result = consensus.buildSequence();

        System.out.println(result);

    }
    public static void test()
    {
        Greedy greed = new Greedy();

        /*
        Sequence s = new Sequence("attagaccatgcggc");
        Sequence t = new Sequence("atcggcattcagt");
        Sequence u = new Sequence("cgtaccgtttacgttt");
        Sequence v = new Sequence("gtacctt");

        List<Sequence> list = new ArrayList<Sequence>();
        list.add(s);
        list.add(t);
        list.add(u);
        list.add(v);*/
        /*

        Sequence f = new Sequence("catagtc");
        Sequence g = new Sequence("taactat");
        Sequence h = new Sequence("agactatcc");

        List<Sequence> list = new ArrayList<Sequence>();

        list.add(f);
        list.add(g);
        list.add(h);*/


        /*
        Sequence f = new Sequence("cacgt");
        Sequence g = new Sequence("acgt");
        Sequence h = new Sequence("actacg");
        Sequence i = new Sequence("gtact");
        Sequence j = new Sequence("actga");
        Sequence k = new Sequence("ctga");

        List<Sequence> list = new ArrayList<Sequence>();

        list.add(f);
        list.add(g);
        list.add(h);
        list.add(i);
        list.add(j);
        list.add(k);*/

        Sequence f = new Sequence("actttacg");
        Sequence g = new Sequence("ttgcacgat");
        Sequence h = new Sequence("ttgcg");
        Sequence i = new Sequence("ggaatctgcgagtta");
        Sequence j = new Sequence("tact");
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
        list.add(_g);*/


        List<SequenceAlignment> result = greed.greedy(list, 1, -1, -2);
        System.out.println(result.size());
        for(SequenceAlignment seq : result)
        {
            System.out.println(seq.s1);
            System.out.println(seq.s2);
            System.out.println("score "+ seq.score);


            System.out.println("****************");
        }

        ConsensusAline consensus = new ConsensusAline(result);
        Sequence cible = consensus.buildSequence();
        System.out.println("tadaaaam");
        System.out.println(cible);



        /*
        System.out.println("Je veux aligner :");
        System.out.println(f);
        System.out.println(h);
        System.out.println(g);
        System.out.println("------->");
        */

        //Consensus c = new Consensus(result);
        //c.showWithOffset();
        //c.showWithEndGapAndOffset();

        /*
        c.computeAlignment();
        Sequence consensus_final = c.build();

        System.out.println("Les séquences alignées:\n");
        ArrayList<Sequence> alignment = c.getAlignment();
        for(int counter = 0;counter < alignment.size();counter++)
            System.out.println(alignment.get(counter));
        System.out.println("#########################################\n");

        System.out.println("Le consensus final:\n");
        System.out.println(consensus_final);

        for(SequenceAlignment s : result)
        {
            System.out.println(s.s1);
            System.out.println(s.s2);
            System.out.println("#####");
        }
        */
    }

    public static void cible(int num, boolean show_consensus, boolean show_alignment, boolean save)
    {
        try
        {
            List<Sequence> list = FastaReader.readFromFile(new File("../../res/collections/Collection" + num + "-Simplifiee.FASTA"));

            System.out.print("Greedy algorithm... ");
            Greedy g = new Greedy();
            List<SequenceAlignment> result = g.greedy(list, 1, -1, -2);
            System.out.println("Done");

            Consensus c = new Consensus(result);
            System.out.print("Alignement... ");
            c.computeAlignment();
            System.out.println("Done");
            //c.showWithEndGapAndOffset();

            System.out.print("Construction du consensus... ");
            Sequence consensus_final = c.build(true);
            System.out.println("Done");

            if (show_alignment)
            {
                System.out.println("Les séquences alignées:\n");
                ArrayList<Sequence> alignment = c.getAlignment();
                for(int counter = 0;counter < alignment.size();counter++)
                    System.out.println(alignment.get(counter));
                System.out.println("#########################################\n");
            }

            if (show_consensus)
            {
                System.out.println("Le consensus final:\n");
                System.out.println(consensus_final);
                System.out.println("#########################################\n");
            }

            if (save)
            {
                String file_name = "Cible_calcul" + num + ".fasta";
                String directory = "../../res/results/";
                System.out.print("Ecriture du consensus dans le fichier " + file_name + "... ");
                FastaWriter.write("Cible1", consensus_final, new File(directory + file_name), 80);
                System.out.println("Done");
            }
        }
        catch (IOException e)
        {
            System.out.println(e.toString());
        }
    }
}
