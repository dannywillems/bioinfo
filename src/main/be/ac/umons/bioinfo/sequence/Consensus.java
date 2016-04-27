package be.ac.umons.bioinfo.sequence;

import java.util.List;

public class Consensus
{
    private List<SequenceAlignment> hamiltonian_path;

    public Consensus(List<SequenceAlignment> hamiltonian_path)
    {
        this.hamiltonian_path = hamiltonian_path;
    }

    public List<SequenceAlignment> addGap()
    {
        return (this.hamiltonian_path);
    }
}
