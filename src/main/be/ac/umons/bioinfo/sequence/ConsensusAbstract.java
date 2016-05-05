package be.ac.umons.bioinfo.sequence;

import java.util.List;
import java.util.ArrayList;

public class ConsensusAbstract
{
    private ArrayList<SequenceAlignmentAbstract> hamiltonian_path;
    private ArrayList<SequenceAbstract> alignment;

    /* ---------------------------------------------------------------------- */
    // CONSTRUCTORS
    public ConsensusAbstract(List<SequenceAlignment> hamiltonian_path)
    {
        this.hamiltonian_path = new ArrayList<SequenceAlignmentAbstract>();
        for(int i = 0;i < hamiltonian_path.size();i++)
            this.hamiltonian_path.add(new SequenceAlignmentAbstract(hamiltonian_path.get(i)));
        this.alignment = new ArrayList<SequenceAbstract>();
    }

    public ConsensusAbstract(ArrayList<SequenceAlignmentAbstract> hp)
    {
        this.hamiltonian_path = hp;
    }
    /* ---------------------------------------------------------------------- */

    /* ---------------------------------------------------------------------- */
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

    public void propageGapsUpPosFrom(int begin, int pos, int nb)
    {
        for (int i = begin - 1;i >= 0;i--)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            sa.s.addGaps(nb, pos);
            sa.t.addGaps(nb, pos);
        }
    }

    public void propageGapsUpFrom(int begin)
    {
        SequenceAlignmentAbstract sa_down = this.getHamiltonianPath().get(begin);
        SequenceAlignmentAbstract sa_up = this.getHamiltonianPath().get(begin - 1);

        SequenceAbstract s = sa_up.t;
        SequenceAbstract t = sa_down.s;

        int real_pos = t.getOffset();
        for (int i = 0;i < t.nb_gaps.length;i++)
        {
            if (t.nb_gaps[i] > s.nb_gaps[i])
                this.propageGapsUpPosFrom(begin, real_pos + 1, t.nb_gaps[i] - s.nb_gaps[i]);
            real_pos += t.nb_gaps[i] + 1;
        }
    }

    public void propageGapsDownPosFrom(int begin, int pos, int nb)
    {
        for (int i = begin + 1;i < this.getHamiltonianPath().size();i++)
        {
            System.out.println("Begin = " + i);
            System.out.println("Must add " + nb + " gaps at pos " + pos);
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            sa.s.addGaps(nb, pos);
            sa.t.addGaps(nb, pos);
        }
    }

    public void propageGapsDownFrom(int begin)
    {
        SequenceAlignmentAbstract sa_up = this.getHamiltonianPath().get(begin);
        SequenceAlignmentAbstract sa_down = this.getHamiltonianPath().get(begin + 1);

        SequenceAbstract s = sa_up.t;
        SequenceAbstract t = sa_down.s;

        int real_pos = s.getOffset();
        for (int i = 0;i < s.nb_gaps.length;i++)
        {
            if (s.nb_gaps[i] > t.nb_gaps[i])
                this.propageGapsDownPosFrom(begin, real_pos + 1, s.nb_gaps[i] - t.nb_gaps[i]);
            real_pos += s.nb_gaps[i] + 1;
        }
    }

    public void computeAlignment()
    {
        this.updateOffset();

        for (int i = 0;i < this.getHamiltonianPath().size() - 1;i++)
            propageGapsDownFrom(i);
        for (int i = this.getHamiltonianPath().size() - 1;i >= 1;i--)
            propageGapsUpFrom(i);
    }

    public void addEndGaps()
    {
        int max = 0;
        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            max = Math.max(sa.s.getSize(), Math.max(sa.t.getSize(), max));
        }

        for(int i = 0;i < this.getHamiltonianPath().size();i++)
        {
            SequenceAlignmentAbstract sa = this.getHamiltonianPath().get(i);
            sa.s.nb_gaps[sa.s.nb_gaps.length - 1] += max - sa.s.getSize();
            sa.t.nb_gaps[sa.t.nb_gaps.length - 1] += max - sa.t.getSize();
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
}
