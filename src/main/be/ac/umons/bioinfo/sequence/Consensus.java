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
    private int[][] offset;

    public Consensus(List<SequenceAlignment> hamiltonian_path)
    {
        this.hamiltonian_path = hamiltonian_path;
        this.alignment = new ArrayList<Sequence>();
        // Check for tests which send a null to initialize a consensus
        if (this.hamiltonian_path != null)
            this.offset = new int[2][this.hamiltonian_path.size()];
    }

    /**
     * Computer the offset for a given hamiltonian path. The offset is defined
     * as the number of gaps to insert at the beginning of the sequence.
     */
    public void computeOffset()
    {
        // Will be useful to recenter the offset to 0 at the end
        int min = 0;
        SequenceAlignment sa = this.hamiltonian_path.get(0);

        int pos_s1 = sa.s1.getPosFirstNucleotide();
        int pos_s2 = sa.s2.getPosFirstNucleotide();
        if (pos_s1 >= pos_s2)
        {
            this.offset[0][0] = pos_s1;
            this.offset[1][0] = 0;
        }
        else
        {
            this.offset[0][0] = 0;
            this.offset[1][0] = pos_s2 - pos_s1;
        }

        for(int i = 1;i < this.hamiltonian_path.size();i++)
        {
            sa = this.hamiltonian_path.get(i);
            pos_s1 = sa.s1.getPosFirstNucleotide();
            pos_s2 = sa.s2.getPosFirstNucleotide();

            // Always >= 0
            this.offset[0][i] = this.offset[1][i - 1];
            if (pos_s1 != 0)
            {
                this.offset[1][i] = pos_s2 - pos_s1 + this.offset[0][i];
                if (this.offset[1][i] < min)
                    min = this.offset[1][i];
            }
            else
                this.offset[1][i] = this.offset[0][i] + pos_s2;
        }

        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            this.offset[0][i] -= min;
            this.offset[1][i] -= min;
        }
    }

    public void showWithEndGapAndOffset()
    {
        this.addOffset();
        this.fillEndWithGaps();
        this.showHamiltonianPath();
    }

    public void showWithOffset()
    {
        this.addOffset();
        this.showHamiltonianPath();

        /*
        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            for(int j = 0;j < this.offset[0][i] - sa.s1.getPosFirstNucleotide();j++)
                System.out.print("-");
            System.out.println(sa.s1);

            for(int j = 0;j < this.offset[1][i] - sa.s2.getPosFirstNucleotide();j++)
                System.out.print("-");
            System.out.println(sa.s2);
        }
        */
    }

    public void propageGapsDownFrom(int beg, int pos)
    {
        SequenceAlignment sa_beg = this.hamiltonian_path.get(beg);
        for(int i = beg + 1;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            sa.s1.addGapAtPos(pos);
            sa.s2.addGapAtPos(pos);
        }
    }

    public void propageGapsUpFrom(int beg, int pos)
    {
        SequenceAlignment sa_beg = this.hamiltonian_path.get(beg);
        for(int i = beg - 1;i >= 0;i--)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            sa.s1.addGapAtPos(pos);
            sa.s2.addGapAtPos(pos);
        }
    }

    /**
     * Propage the gaps from begin to last sequence. IMPROVEME!
     */
    public void propageGaps(int begin, int[][] gaps)
    {
        int move_down = 0;
        int move_up = 0;

        int i = 0;
        int j = 0;
        if (gaps[0][i] < gaps[1][j])
        {
            this.propageGapsDownFrom(begin, gaps[0][i] + move_down);
            move_up++;
        }
        else
        {
            this.propageGapsUpFrom(begin, gaps[1][i] + move_up);
            move_down++;
        }
    }

    /**
     * Fill end with gaps. IMPROVEME!!
     */
    public void fillEndWithGaps()
    {
        int max = 0;
        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            int size = sa.s1.getSize();
            if (size > max)
                max = size;
            size = sa.s2.getSize();
            if (size > max)
                max = size;
        }

        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            byte[] content = new byte[max];
            StringBuilder s = new StringBuilder();

            int size = sa.s1.getSize();
            for(int j = 0;j < size;j++)
                s.append(sa.s1.getLetter(j));
            for(int j = 0;j < max - size;j++)
                s.append(Sequence.base2letter((byte) Sequence.GAP));
            sa.s1.setContent(s.toString());

            size = sa.s2.getSize();
            s = new StringBuilder();
            for(int j = 0;j < size;j++)
                s.append(sa.s2.getLetter(j));
            for(int j = 0;j < max - size;j++)
                s.append(Sequence.base2letter((byte) Sequence.GAP));
            sa.s2.setContent(s.toString());
        }
    }

    public void addOffset()
    {
        this.computeOffset();
        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            int pos_s1 = sa.s1.getPosFirstNucleotide();
            int pos_s2 = sa.s2.getPosFirstNucleotide();
            StringBuilder s = new StringBuilder();

            for(int j = 0;j < offset[0][i] - pos_s1;j++)
                s.append(Sequence.base2letter((byte) Sequence.GAP));
            for(int j = 0;j < sa.s1.getSize();j++)
                s.append(sa.s1.getLetter(j));
            sa.s1.content = Sequence.letter2Base(s.toString());

            s = new StringBuilder();
            for(int j = 0;j < offset[1][i] - pos_s2;j++)
                s.append(Sequence.base2letter((byte) Sequence.GAP));
            for(int j = 0;j < sa.s2.getSize();j++)
                s.append(sa.s2.getLetter(j));
            sa.s2.content = Sequence.letter2Base(s.toString());
        }
    }

    /**
     * Do the alignment and save it in alignment attribute.
     */
    public void computeAlignment()
    {
        // Add offset and end gaps to work easier. IMPROVEME!
        this.addOffset();
        this.fillEndWithGaps();

        for(int i = 0;i < this.hamiltonian_path.size() - 1;i++)
        {
            System.out.println(i);
            SequenceAlignment sa_up = this.hamiltonian_path.get(i);
            SequenceAlignment sa_down = this.hamiltonian_path.get(i + 1);
            SequencePairSame pair = new SequencePairSame(sa_up.initial_s2, sa_up.s2, sa_down.s1);

            /*
            System.out.println(pair.initial.toString());
            System.out.println(pair.s.toString());
            System.out.println(pair.t.toString());
            */
            int[][] gaps = pair.findGaps();
            //pair.rebuildWithGaps(); // Done in the propagation
            this.propageGaps(i, gaps);
        }

        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            if (i == 0)
                this.alignment.add(this.hamiltonian_path.get(i).s1);
            this.alignment.add(this.hamiltonian_path.get(i).s2);
        }
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

    public int[][] getOffset()
    {
        return (this.offset);
    }

    public void setAlignment(ArrayList<Sequence> a)
    {
        this.alignment = a;
    }

    public ArrayList<Sequence> getAlignment()
    {
        return (this.alignment);
    }

    public List<SequenceAlignment> getHamiltonianPath()
    {
        return (this.hamiltonian_path);
    }

    public void showHamiltonianPath()
    {
        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            SequenceAlignment sa = this.hamiltonian_path.get(i);
            System.out.println(sa.s1.toString());
            System.out.println(sa.s2.toString());
        }
    }

    public boolean sameHamiltonianPath(Consensus other)
    {
        List<SequenceAlignment> hp_o = other.getHamiltonianPath();
        if (hp_o.size() == this.hamiltonian_path.size())
        {
            for(int i = 0;i < this.hamiltonian_path.size();i++)
            {
                SequenceAlignment sa = this.hamiltonian_path.get(i);
                SequenceAlignment sa_o = hp_o.get(i);
                if (!sa.s1.equals(sa_o.s1) || !sa.s2.equals(sa_o.s2))
                    return (false);
            }
            return (true);
        }
        return (false);
    }
}
