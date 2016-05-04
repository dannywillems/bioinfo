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
    public void computeAlignment()
    {

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
}
