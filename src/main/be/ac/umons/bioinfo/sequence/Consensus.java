package be.ac.umons.bioinfo.sequence;

import java.util.List;
import java.util.ArrayList;
import java.lang.Integer;
import java.lang.Character;
import java.util.HashMap;
import java.util.Iterator;

public class Consensus
{
    private List<SequenceAlignment> hamiltonian_path;
    private ArrayList<Sequence> alignment;

    public Consensus(List<SequenceAlignment> hamiltonian_path)
    {
        this.hamiltonian_path = hamiltonian_path;
        this.alignment = new ArrayList<Sequence>();
    }

    /**
     * Computer the offset for a given hamiltonian path. The offset is defined
     * as the number of gaps to insert at the beginning of the sequence.
     */
    public int[][] computeOffset()
    {
        int offset[][] = new int[2][this.hamiltonian_path.size()];
        // Will be useful to recenter the offset to 0 at the end
        int min = 0;
        SequenceAlignment sa = this.hamiltonian_path.get(0);

        int pos_s1 = sa.s1.getPosFirstNucleotide();
        int pos_s2 = sa.s2.getPosFirstNucleotide(); 
        if (pos_s1 >= pos_s2)
        {
            offset[0][0] = pos_s1;
            offset[1][0] = 0;
        }
        else
        {
            offset[0][0] = 0;
            offset[1][0] = pos_s2 - pos_s1;
        }

        for(int i = 1;i < this.hamiltonian_path.size();i++)
        {
            sa = this.hamiltonian_path.get(i);
            pos_s1 = sa.s1.getPosFirstNucleotide();
            pos_s2 = sa.s2.getPosFirstNucleotide();

            // Always >= 0
            offset[0][i] = offset[1][i - 1];
            if (pos_s1 != 0)
            {
                offset[1][i] = pos_s2 - pos_s1 + offset[0][i];
                if (offset[1][i] < min)
                    min = offset[1][i];
            }
            else
                offset[1][i] = offset[0][i] + pos_s2;
        }

        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            offset[0][i] -= min;
            offset[1][i] -= min;
        }

        return (offset);
    }

    public void showWithOffset()
    {
        int[][] gaps = this.computeOffset();
        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            for(int j = 0;j < gaps[0][i] - sa.s1.getPosFirstNucleotide();j++)
                System.out.print("-");
            System.out.println(sa.s1);

            for(int j = 0;j < gaps[1][i] - sa.s2.getPosFirstNucleotide();j++)
                System.out.print("-");
            System.out.println(sa.s2);
        }
    }

    /**
     * Do the alignment and save it in alignment attribute.
     */
    public void computeAlignment()
    {
        int[][] offset = this.computeOffset();
    }

    /**
     * Build the consensus based on the alignment done in the computeAlignment method.
     * Use the alignment attribute so you need to use this.computeAlignment before!
     * @return the consensus as a Sequence object
     */
    public Sequence build(boolean remove_if_max_gap)
    {
        StringBuilder consensus = new StringBuilder();
        int size_consensus = this.alignment.get(0).getSize();
        for (int i = 0;i < size_consensus;i++)
        {
            HashMap<Character, Integer> occurences = new HashMap<Character, Integer>();
            for (int j = 0;j < alignment.size();j++)
            {
                Character c = new Character(alignment.get(j).getLetter(i));
                Integer c_occurence = occurences.get(c);
                if (c_occurence == null)
                    occurences.put(c, 0);
                else
                    occurences.put(c, c_occurence + 1);
            }
            char base = this.getBase(occurences, remove_if_max_gap);
            // if remove_if_max_gap is true, we can have a gap as a maximum but
            // we don't add it
            if (remove_if_max_gap)
            {
                if (base != Sequence.base2letter((byte) Sequence.GAP))
                    consensus.append(base);
            }
            // else, we know the max is not a gap because getBase manage this
            // case. So we can add without regarding the type.
            else
                consensus.append(base);
        }
        return (new Sequence(consensus.toString()));
    }

    /**
     * Return the more occured character at a position.
     * If equality, it takes the first. If a gap is more occured than others, it
     * takes in an undefined ordre (due to HashMap, not a feature). FIXME
     * It supposes that there is another character than the gap.
     */
    public char getBase(HashMap<Character, Integer> o, boolean remove_if_max_gap)
    {
        int max = Integer.MIN_VALUE;
        char gap = Sequence.base2letter((byte) Sequence.GAP);
        char c_max = gap;

        Iterator<Character> iterator_o = o.keySet().iterator();
        while (iterator_o.hasNext())
        {
            Character c = iterator_o.next();
            if (o.get(c).intValue() > max)
            {
                // if remove_if_max_gap is false, we don't have to take the gap.
                if (remove_if_max_gap || (c.charValue() != gap))
                {
                    max = o.get(c).intValue();
                    c_max = c.charValue();
                }
            }
        }

        //assert c_max != gap : "There must never been a gap in the final consensus";
        return (c_max);
    }

    public void setAlignment(ArrayList<Sequence> a)
    {
        this.alignment = a;
    }

    public ArrayList<Sequence> getAlignment()
    {
        return (this.alignment);
    }
}
