package be.ac.umons.bioinfo.sequence;

import java.util.List;
import java.util.ArrayList;

public class ConsensusAbstract
{
    public ArrayList<SequenceAlignmentAbstract> hamiltonian_path;
    public ArrayList<SequenceAbstract> alignment;

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
}
