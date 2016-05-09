package be.ac.umons.bioinfo.consensus;

import java.util.HashMap;
import java.util.Map;
import be.ac.umons.bioinfo.sequence.*;

/**
 * Created by aline on 6/05/16.
 */
public class MyCounter
{
    private Map<Byte, Integer> counts;

    public MyCounter()
    {
        this.counts = new HashMap<>();
    }

    public void vote(byte nucleotide)
    {
        int old = this.counts.getOrDefault(nucleotide, 0);
        this.counts.put(nucleotide, old+1);
    }

    public byte max()
    {
        int max = Integer.MIN_VALUE;
        byte nucleotide = 0;

        for(Map.Entry<Byte, Integer> entry : counts.entrySet())
        {
            if(entry.getValue() > max && entry.getKey() != Nucleotide.GAP)
            {
                max = entry.getValue();
                nucleotide = entry.getKey();
            }
        }

        return nucleotide;
    }
}
