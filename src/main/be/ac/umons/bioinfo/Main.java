package be.ac.umons.bioinfo;

/**
 * Created by aline on 31/03/16.
 */
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import be.ac.umons.bioinfo.sequence.*;

public class Main
{
    public static void main(String[] args) throws IOException
    {
        Greedy greed = new Greedy();

        Sequence f = new Sequence("catagtc");
        Sequence g = new Sequence("taactat");
        Sequence h = new Sequence("agactatcc");

        List<Sequence> list = new ArrayList<>();
        list.add(f);
        list.add(g);
        list.add(h);

        List<SequenceAlignment> result = greed.computePath(list, 1, -1, -2);
        Iterator<SequenceAlignment> iterator = result.iterator();

        while(iterator.hasNext())
        {
            SequenceAlignment alignment = iterator.next();
            System.out.println(alignment.s1);
            System.out.println(alignment.s2);
            System.out.println("----------------");
        }


    }
}
