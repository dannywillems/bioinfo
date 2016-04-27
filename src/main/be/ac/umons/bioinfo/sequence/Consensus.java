package be.ac.umons.bioinfo.sequence;

import java.util.List;
import java.lang.Integer;
import java.lang.Character;
import java.util.HashMap;
import java.util.Iterator;

public class Consensus
{
    private List<SequenceAlignment> hamiltonian_path;
    private List<Sequence> alignment;

    public Consensus(List<SequenceAlignment> hamiltonian_path)
    {
        this.hamiltonian_path = hamiltonian_path;
    }

    /**
     * Do the alignment and save it in alignment attribute.
     */
    public void alignment()
    {
    }

    /**
     * Build the consensus based on the alignment done in the alignment method.
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
                if (c_occurence.equals(new Integer(0)))
                    occurences.put(c, new Integer(0));
                else
                    occurences.put(c, new Integer(c_occurence.intValue() + Integer.valueOf(1)));
            }
            consensus.append(this.getBase(occurences));
        }
        return (new Sequence(consensus.toString()));
    }

    /**
     * Return the more occured character at a position.
     * If equality, it takes the first. If a gap is more occured than others, it
     * takes the first next more occured character.
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
        return (c_max);
    }

    public void setAlignment(List<Sequence> a)
    {
        this.alignment = a;
    }
}
