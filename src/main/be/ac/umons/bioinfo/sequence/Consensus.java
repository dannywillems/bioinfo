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
     * Do the alignment and save it in alignment attribute.
     */
    public void computeAlignment()
    {
        int[] gaps = new int[2 * this.hamiltonian_path.size()];

        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            if (i == 0)
                this.alignment.add(this.hamiltonian_path.get(i).s1);
            this.alignment.add(this.hamiltonian_path.get(i).s2);
        }

        //this.alignment = this.addGap(gaps);
    }

    /**
     * Add gaps to the initial alignments.
     * @param gaps int array containing at the position 2i the number of gaps to
     * add at the beginning of the ith sequence and at position 2i + 1 the
     * numbers of gaps to add at the end of the ith sequence.
     * @return the alignment with gaps
     */
    private ArrayList<Sequence> addGap(int[] gaps)
    {
        ArrayList<Sequence> final_alignment = new ArrayList<Sequence>();
        for(int i = 0;i < this.alignment.size();i++)
        {
            Sequence current = this.alignment.get(i);
            StringBuilder s = new StringBuilder();
            for(int j = 0;j < gaps[2 * i];j++)
                s.append(Sequence.base2letter((byte) Sequence.GAP));
            for(int j = 0;j < current.getSize();j++)
                s.append(current.getLetter(j));
            for(int j = 0;j < gaps[2 * i + 1];j++)
                s.append(Sequence.base2letter((byte) Sequence.GAP));
            final_alignment.add(new Sequence(s.toString()));
        }
        return (final_alignment);
    }

    /**
     * Build the consensus based on the alignment done in the computeAlignment method.
     * Use the alignment attribute.
     * @return the consensus as a Sequence object
     */
    public Sequence build()
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
            consensus.append(this.getBase(occurences));
        }
        return (new Sequence(consensus.toString()));
    }

    /**
     * Return the more occured character at a position.
     * If equality, it takes the first. If a gap is more occured than others, it
     * takes in an undefined ordre (due to HashMap, not a feature). FIXME
     * It supposes that there is another character than the gap.
     */
    public char getBase(HashMap<Character, Integer> o)
    {
        int max = Integer.MIN_VALUE;
        char gap = Sequence.base2letter((byte) Sequence.GAP);
        char c_max = gap;

        Iterator<Character> iterator_o = o.keySet().iterator();
        while (iterator_o.hasNext())
        {
            Character c = iterator_o.next();
            if (o.get(c).intValue() > max && c.charValue() != gap)
            {
                max = o.get(c).intValue();
                c_max = c.charValue();
            }
        }

        assert c_max != gap : "There must never been a gap in the final consensus";
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
