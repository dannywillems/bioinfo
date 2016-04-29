package be.ac.umons.bioinfo;

/**
 * Created by aline on 31/03/16.
 */
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import be.ac.umons.bioinfo.sequence.*;

public class Main
{
    public static void main(String[] args) throws IOException
    {

        Greedy greed = new Greedy();

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
        Sequence k = new Sequence("ctga");*/




        Sequence f = new Sequence("actttacg");
        Sequence g = new Sequence("ttgcacgat");
        Sequence h = new Sequence("ttgcg");
        Sequence i = new Sequence("ggaatctgcgagtta");
        Sequence j = new Sequence("tact");
        Sequence k = new Sequence("gaccgat");

        List<Sequence> list = new ArrayList<Sequence>();

        list.add(j);
        list.add(k);
        list.add(f);
        list.add(i);
        list.add(g);
        list.add(h);


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

        /*
        System.out.println("Je veux aligner :");
        System.out.println(f);
        System.out.println(h);
        System.out.println(g);
        System.out.println("------->");
*/

        for(SequenceAlignment sa : result)
        {
            System.out.println(sa.s1);
            System.out.println(sa.s2);
            System.out.println(sa.score);

            System.out.println("****************");
        }


    }
}
