package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.Debug;
import be.ac.umons.bioinfo.fasta.FastaWriter;
import be.ac.umons.bioinfo.sequence.*;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.SortedSet;
import java.util.Iterator;

import java.io.IOException;
import java.io.File;

public class ConsensusAbstract
{
    private ArrayList<SequenceAlignmentAbstract> hamiltonian_path;
    private ArrayList<SequenceAbstract> alignment;

    /* ---------------------------------------------------------------------- */
    // CONSTRUCTORS
    /**
     * Create a ConsensusAbstract object from a hamilonian path of
     * SequenceAlignment object.
     * Complexity: ??? Be careful by using ArrayList.add !
     */
    public ConsensusAbstract(List<SequenceAlignment> hamiltonian_path)
    {
        this.hamiltonian_path = new ArrayList<SequenceAlignmentAbstract>();
        for(int i = 0;i < hamiltonian_path.size();i++)
            this.hamiltonian_path.add(new SequenceAlignmentAbstract(hamiltonian_path.get(i)));

        this.alignment = new ArrayList<SequenceAbstract>();
        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            if (i == 0)
                this.alignment.add(sa.s);
            this.alignment.add(sa.t);
        }
    }

    /**
     * Create a ConsensusAbstract object from an ArrayList of
     * SequenceAlignmentAbstract.
     * Complexity: O(n) where n is the size of the hamiltonian path due to
     * initialisation of the alignment attribute.
     */
    public ConsensusAbstract(ArrayList<SequenceAlignmentAbstract> hp)
    {
        this.hamiltonian_path = hp;

        // Initialize alignment to be able to use ArrayList.set which is in O(1)
        // contrary to ArrayList.add which is O(n). So it takes much time to
        // initialize a Consensus but the updateAlignment method is in O(n)
        // where n is the hamiltonian path size, not in O(n^2) (due to insertion
        // in O(n)).
        this.alignment = new ArrayList<SequenceAbstract>();
        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            if (i == 0)
                this.alignment.add(sa.s);
            this.alignment.add(sa.t);
        }
    }

    /**
     * Create a ConsensusAbstract object from an ArrayList of SequenceAbstract
     * which will be set as the alignment.
     * It can only be used to build the consensus, not compute it because we
     * don't have an hamiltonian path.
     */
    public ConsensusAbstract()
    {
        this.alignment = new ArrayList<SequenceAbstract>();
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /**
     * Update the offset of the sequences in the hamiltonian path.
     * Complexity: O(n) where n is the size of the hamiltonian path. Using
     * ArrayList.get, SequenceAbstract.getOffset, SequenceAbstract.setOffset is
     * in O(1).
     */
    public void updateOffset()
    {
        // Will be useful to recenter the offset to 0 at the end
        int min = 0;
        SequenceAlignmentAbstract sa = this.hamiltonian_path.get(0);
        SequenceAlignmentAbstract sa_previous = this.hamiltonian_path.get(0);

        /* ------------------------------------------------------------------ */
        /* basis case */
        // The first line is used as basis. A recurrence rule will be used.
        int pos_s = sa.s.getOffset();
        int pos_t = sa.t.getOffset();
        if (pos_s >= pos_t)
            // if s begins with gaps, pos_t doesn't begin with gaps
            sa.t.setOffset(0);
        else
            // Else, is t which begins with gaps and not s.
            sa.t.setOffset(pos_t);
        /* ------------------------------------------------------------------ */

        /* ------------------------------------------------------------------ */
        /* Recurrence rule */
        for(int i = 1;i < this.hamiltonian_path.size();i++)
        {
            sa_previous = sa;
            sa = this.hamiltonian_path.get(i);
            pos_s = sa.s.getOffset();
            pos_t = sa.t.getOffset();

            // Always >= 0
            sa.s.setOffset(sa_previous.t.getOffset());
            if (pos_s != 0)
            {
                sa.t.setOffset(pos_t - pos_s + sa.s.getOffset());
                min = Math.min(sa.t.getOffset(), min);
            }
            else
                sa.t.setOffset(sa.s.getOffset() + pos_t);
        }

        for(int i = 0;i < this.hamiltonian_path.size();i++)
        {
            sa = this.getHamiltonianPath().get(i);
            sa.s.setOffset(sa.s.getOffset() - min);
            sa.t.setOffset(sa.t.getOffset() - min);
        }
        /* ------------------------------------------------------------------ */
    }

    /**
     * Propage up [nb] gaps from [begin] to the beginning of the hamiltonian
     * path at the position [pos].
     * Complexity: ???
     */
    public void propageGapsUpPosFrom(int begin, int indice, int nb)
    {
        int pos = 0;
        for (int i = begin - 1;i >= 0;i--)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            //System.out.println("On doit propager à l'indice " + indice + " dans " + sa.t.toString());
            //pos = sa.t.addGapsAfterIndiceAndReturnPosition(nb, indice);
            pos = sa.t.addGapsAfterIndiceEndAndReturnPosition(nb, indice);
            //System.out.println("On doit propager à la position " + pos + " dans " + sa.s.toString());
            indice = sa.s.addGapsAndReturnIndice(nb, pos);
        }
    }

    /**
     * Propage up gaps from [begin] to the beginning of the hamiltonian path.
     * Complexity: ???
     */
    public void propageGapsUpFrom(int begin)
    {
        SequenceAlignmentAbstract sa_up = this.getHamiltonianPath().get(begin - 1);
        SequenceAlignmentAbstract sa_down = this.getHamiltonianPath().get(begin);

        SequenceAbstract s = sa_up.t;
        SequenceAbstract t = sa_down.s;

        for (int i = 0;i < t.nb_gaps.length - 1;i++)
        {
            if (t.nb_gaps[i] > s.nb_gaps[i])
            {
                this.propageGapsUpPosFrom(begin, i, t.nb_gaps[i] - s.nb_gaps[i]);
                this.checkEqualitySameGapsNumberDuringPropagation(begin - 1, i);
            }
        }
        this.checkPropageGapDownFrom(begin - 1);
    }

    /**
     * Propage down [nb] gaps from [begin] to the end of the hamiltonian
     * path at the position [pos].
     * Complexity: ???
     */
    public void propageGapsDownPosFrom(int begin, int indice, int nb)
    {
        int pos = 0;
        for (int i = begin + 1;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            pos = sa.s.addGapsAfterIndiceEndAndReturnPosition(nb, indice);
            indice = sa.t.addGapsAndReturnIndice(nb, pos);
        }
    }

    /**
     * Propage down gaps from [begin] to the end of the hamiltonian path.
     * Complexity: ???
     */
    public void propageGapsDownFrom(int begin)
    {
        SequenceAlignmentAbstract sa_up = this.getHamiltonianPath().get(begin);
        SequenceAlignmentAbstract sa_down = this.getHamiltonianPath().get(begin + 1);

        SequenceAbstract s = sa_up.t;
        SequenceAbstract t = sa_down.s;

        for (int i = 0;i < s.nb_gaps.length - 1;i++)
        {
            if (s.nb_gaps[i] > t.nb_gaps[i])
            {
                this.propageGapsDownPosFrom(begin, i, s.nb_gaps[i] - t.nb_gaps[i]);
                this.checkEqualitySameGapsNumberDuringPropagation(begin, i);
            }
        }
        this.checkPropageGapDownFrom(begin);
    }

    /**
     * Compute the alignment by
     * - update the offset
     * - propage down gaps
     * - propage up gaps
     * - update the alignment
     * - add ending gaps to the hamiltonian path.
     * Complexity: ???
     */
    public void computeAlignment()
    {
        try
        {
            this.updateOffset();

            //System.out.println("On propage vers le bas...");
            for (int i = 0;i < this.getHamiltonianPath().size() - 1;i++)
            {
                this.checkOffsetSame("Offset avant propagation bas: ", i);
                this.propageGapsDownFrom(i);
                this.checkOffsetSame("Offset après propagation bas: ", i);
            }
            FastaWriter.writeAlignment(this.getHamiltonianPath(), new File("../../hp_after_down_" + Debug.count_writer + ".fasta"), 80);

            //System.out.println("Résultat après propagation vers le bas:");
            //this.showHamiltonianPath();
            //System.out.println("########################################");

            //System.out.println("On propage vers le haut...");
            for (int i = this.getHamiltonianPath().size() - 1;i >= 1;i--)
            {
                this.checkOffsetSame("Offset avant propagation haut: ", i - 1);
                this.propageGapsUpFrom(i);
                this.checkOffsetSame("Offset avant propagation haut: ", i - 1);
            }
            FastaWriter.writeAlignment(this.getHamiltonianPath(), new File("../../hp_after_up_" + Debug.count_writer + ".fasta"), 80);
            Debug.count_writer++;

            //System.out.println("Résultat après propagation vers le haut:");
            //this.showHamiltonianPath();
            //System.out.println("########################################");
            /* To propage up and down at the same time. NOT WORKING. FIXME -->
             * ArrayOutOfBounds.
            for (int i = 0;i < this.getHamiltonianPath().size();i++)
            {
                if (i == 0)
                    this.propageGapsDownFrom(i);
                else if (i == this.getHamiltonianPath().size())
                    this.propageGapsUpFrom(i);
                else
                {
                    this.propageGapsDownFrom(i);
                    this.propageGapsUpFrom(i);
                }
            }
            */

            this.addEndGaps();
            this.updateAlignment();

            this.checkEqualitySameGapsNumber();
        }
        catch (IOException e)
        {
            System.out.println(e);
        }
    }

    /**
     * Update the alignment attribute based on the hamiltonian path. It only
     * consists of ranging the hamiltonian path, taking the first sequence for
     * the first alignment and the second for all alignments.
     * Complexity: O(n) where n is the size of the hamiltonian path due to
     * ArrayList.set method used to update the alignment attribute and to
     * ArrayList.get to get the i-th element of the hamiltonian path.
     * ArrayList.size is also in O(1).
     */
    public void updateAlignment()
    {
        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            if (i == 0)
                this.alignment.set(i, sa.s);
            this.alignment.set(i + 1, sa.t);
        }
    }

    /**
     * Add ending gaps to the sequences in the hamiltonian path.
     * Complexity: O(n + m) where n is the number of sequences in the hamiltonian
     * path and where m is the compressed size (see getSize in SequenceAbstract
     * of the bigger sequence. The m is due to the computation of the size at
     * (*). NB: get on ArrayList is O(1).
     */
    public void addEndGaps()
    {
        int max = 0;
        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            // IMPROVEME: compute the size is not in O(1) but in O(m) where m
            // is the compressed size. Can be improve by saving the size in the
            // SequenceAbstract object. (*)
            max = Math.max(sa.s.getSize(), Math.max(sa.t.getSize(), max));
        }

        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            // IMPROVEME: compute the size is not in O(1) but in O(m) where m
            // is the compressed size. Can be improve by saving the size in the
            // SequenceAbstract object. (*)
            sa.s.nb_gaps[sa.s.nb_gaps.length - 1] += max - sa.s.getSize();
            sa.t.nb_gaps[sa.t.nb_gaps.length - 1] += max - sa.t.getSize();
        }
    }

    /**
     * Add ending gaps to the sequences in the alignment.
     * Complexity: O(n + m) where n is the number of sequences in the alignment
     * and where m is the compressed size (see getSize in SequenceAbstract
     * of the bigger sequence. The m is due to the computation of the size at
     * (*). NB: get on ArrayList is O(1).
     * @deprecated:
     *  Use addEndGaps following by updateAlignment ==> Same
     *  complexity and update at the same time all sequences
     */
    @Deprecated public void addEndGapsAlignment()
    {
        int max = 0;
        for(int i = 0;i < this.getAlignment().size();i++)
            // IMPROVEME: compute the size is not in O(1) but in O(m) where m
            // is the compressed size. Can be improve by saving the size in the
            // SequenceAbstract object. (*)
            max = Math.max(this.getAlignment().get(i).getSize(), max);

        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAbstract s = this.getAlignment().get(i);
            // IMPROVEME: compute the size is not in O(1) but in O(m) where m
            // is the compressed size. Can be improve by saving the size in the
            // SequenceAbstract object. (*)
            s.nb_gaps[s.nb_gaps.length - 1] += max - s.getSize();
        }
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    // GETTERS
    public ArrayList<SequenceAbstract> getAlignment()
    {
        return (this.alignment);
    }

    public ArrayList<SequenceAlignmentAbstract> getHamiltonianPath()
    {
        return (this.hamiltonian_path);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    // SETTERS
    public void setAlignment(ArrayList<SequenceAbstract> a)
    {
        this.alignment = a;
    }

    public void setHamiltonianPath(ArrayList<SequenceAlignmentAbstract> hp)
    {
        this.hamiltonian_path = hp;
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    // SHOW METHODS
    public void showHamiltonianPath()
    {
        for(int i = 0;i < this.hamiltonian_path.size();i++)
            this.hamiltonian_path.get(i).showAlignment();
    }

    public void showAlignment()
    {
        for(int i = 0;i < this.hamiltonian_path.size();i++)
            System.out.println(this.alignment.get(i).toString());
    }

    public void showWithOffset()
    {
        this.updateOffset();
        this.showHamiltonianPath();
    }

    public void showWithEndGapsAndOffset()
    {
        this.updateOffset();
        this.addEndGaps();
        this.showHamiltonianPath();
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /**
     * Check if two ConsensusAbstract object have the same hamiltonian paths.
     * Complexity: O(n + m) where n is the size of the hamiltonian path and m
     * the compressed size of the longer sequence in the hamiltonian path. The m
     * is due to the equals function.
     * @return true if the ConsensusAbstract objects have the same hamiltonian
     * paths else false.
     */
    public boolean sameHamiltonianPath(ConsensusAbstract other)
    {
        boolean equal = this.getHamiltonianPath().size() == other.getHamiltonianPath().size();
        int i = 0;
        while (equal && i < this.getHamiltonianPath().size())
            equal &=  this.getHamiltonianPath().get(i).s.equals(other.getHamiltonianPath().get(i).s)
                  &&  this.getHamiltonianPath().get(i).t.equals(other.getHamiltonianPath().get(i++).t);
        return (equal);
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
    /** Build the consensus.
     * FIXME: Convert all sequences in the aligment in strings to be able to use
     * the old algorithm which works on strings and char. We can think about an
     * Iterator.
     */
    public SequenceAbstract build(boolean remove_if_max_gap)
    {
        StringBuilder consensus = new StringBuilder();

        ArrayList<String> s_list = new ArrayList<String>();
        for (int i = 0;i < this.alignment.size();i++)
            s_list.add(this.alignment.get(i).toString());

        int size_consensus = this.alignment.get(0).getSize();
        for (int i = 0;i < size_consensus;i++)
        {
            HashMap<Character, Integer> occurences = new HashMap<Character, Integer>();
            // We browse all sequences in the alignment at pos i.
            for (int j = 0;j < alignment.size();j++)
            {
                Character c = new Character(s_list.get(j).charAt(i));
                Integer c_occurence = occurences.get(c);
                if (c_occurence == null)
                    occurences.put(c, 1);
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
        return (new SequenceAbstract(consensus.toString()));
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

        SortedSet<Character> keys = new TreeSet<Character>(o.keySet());
        char c_max = keys.first().charValue();

        for (Character c : keys)
        {
            if (o.get(c).intValue() > max)
            {
                if (!remove_if_max_gap)
                {
                    if (c.charValue() != gap)
                    {
                        max = o.get(c).intValue();
                        c_max = c.charValue();
                    }
                }
                else
                {
                    max = o.get(c).intValue();
                    c_max = c.charValue();
                }
            }
        }

        //assert c_max != gap : "There must never been a gap in the final consensus";
        if (Debug.GAP_IN_CONSENSUS)
        {
            if (c_max == gap && !remove_if_max_gap)
            {
                System.out.println("GAP renvoyée");
                System.out.println("Il y a " + o.keySet().size() + " nucléotides dans cette colonne.");
            }
        }

        return (c_max);
    }

    /* ---------------------------------------------------------------------- */
    // Check functions
    public void checkOffsetSame(String str, int i)
    {
        if (Debug.OFFSET)
        {
            System.out.print(str);
            if (this.hamiltonian_path.get(i).t.getOffset() == this.hamiltonian_path.get(i + 1).s.getOffset())
                System.out.println("Offset: OK");
            else
                System.out.println("Offset: ERROR!!!!!");
        }
    }

    public void checkEqualitySameGapsNumber()
    {
        if (Debug.EQUALITY_SAME_GAPS_NUMBER)
        {
            for(int i = 0;i < this.hamiltonian_path.size() - 1;i++)
            {
                if (this.hamiltonian_path.get(i).t.hasSameGapsNumber(this.hamiltonian_path.get(i + 1).s))
                    System.out.println("Same gaps number: OK");
                else
                {
                    System.out.println(this.hamiltonian_path.get(i).t);
                    System.out.println(this.hamiltonian_path.get(i + 1).s);
                    System.out.println("Same gaps number : ERROR!!!!!");
                }
            }
        }
    }

    public void checkEqualitySameGapsNumberDuringPropagation(int begin, int pos)
    {
        if (Debug.EQUALITY_SAME_GAPS_NUMBER_DURING_PROPAGATION)
        {
            if (this.hamiltonian_path.get(begin).t.nb_gaps[pos] == this.hamiltonian_path.get(begin + 1).s.nb_gaps[pos])
                System.out.println("Same gaps number during propagation: OK");
            else
            {
                System.out.println(this.hamiltonian_path.get(begin).t);
                System.out.println(this.hamiltonian_path.get(begin + 1).s);
                System.out.println("Same gaps number during propagation: ERROR!!!!!");
            }
        }
    }

    public void checkPropageGapDownFrom(int begin)
    {
        if (Debug.PROPAGATION_GAP_DOWN_FROM)
        {
            for(int i = 0;i < this.getHamiltonianPath().get(begin).t.nb_gaps.length - 1;i++)
            {
                if (this.getHamiltonianPath().get(begin).t.nb_gaps[i] > this.getHamiltonianPath().get(begin + 1).s.nb_gaps[i])
                {
                    System.out.println("ERROR!!!! Must not have more gaps in the sequence up");
                    System.out.println("Found at: " + begin);
                }
            }
        }
    }
    /* ---------------------------------------------------------------------- */
}
