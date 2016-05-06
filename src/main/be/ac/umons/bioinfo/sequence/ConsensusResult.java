package be.ac.umons.bioinfo.sequence;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by aline on 5/05/16.
 */
public class ConsensusResult {

    public HashMap<Integer,Counter> map;
    public int posMin;
    public int posMax;

    public ConsensusResult(HashMap<Integer,Counter> map, int min, int max)
    {
        this.map = map;
        this.posMin = min;
        this.posMax = max;
    }
}
