package be.ac.umons.bioinfo.consensus;

import java.util.*;
import be.ac.umons.bioinfo.sequence.*;

/**
 * Created by aline on 5/05/16.
 */
public class ConsensusAline {

    public List<SequenceAlignment> path;

    public ConsensusAline(List<SequenceAlignment> path){
        this.path = path;
    }


    /**
     *
     * @param votes the votes associated to each position of the consensus sequences. Must not be empty.
     * @return The consensus sequence resulting from the votes.
     */
    public static Sequence produceConsensusSequence(Map<Integer, MyCounter> votes)
    {
        int min = votes.keySet().stream().mapToInt(Integer::intValue).min().getAsInt();
        int max = votes.keySet().stream().mapToInt(Integer::intValue).max().getAsInt();

        byte[] result = new byte[max-min+1];

        for(int pos=min ; pos <= max ; pos++)
            result[pos-min] = votes.get(pos).max();

        return new Sequence(result);
    }

    public HashMap<Integer, Counter> computeAlignementWithoutGapPropagation(List<SequenceAlignment> path) {

        HashMap<Integer, Counter> count = new HashMap<Integer, Counter>();
        int posMin = Integer.MAX_VALUE; // the minimum position in the map count
        int posMax = Integer.MIN_VALUE; // the maximum position in the map count
        // these two variables allows to know if a counter already exists for a specific offset
        int offset = 0;

        for (SequenceAlignment alignment : path) {


            Sequence s1 = alignment.s1;
            Sequence s2 = alignment.s2;
            int dec = 0;
            if (s1.content[0] == Nucleotide.GAP) {
                while (s1.content[dec] == Nucleotide.GAP) {
                    offset -= 1;
                    dec += 1;
                }
            }

            int pos = offset;

            for (int i = 0; i <= s1.content.length - 1; i++)
            {

                if (!count.containsKey(pos))
                {
                    Counter c = new Counter();
                    c.vote(s1.content[i]);
                    count.put(pos, c);

                    posMin = Math.min(posMin, pos);
                    posMax = Math.max(posMax, pos);

                } else {
                    count.get(pos).vote(s2.content[i]);
                    pos += 1;

                }

            }
            if (s1.content[0] != Nucleotide.GAP)
            {
                int j = 0;
                while(s2.content[j] == Nucleotide.GAP)
                {
                    offset += 1;
                    j += 1;
                }
            }
        }
        return count;

    }


    public HashMap<Integer,Counter> computeConsensus(){

        //stock the counter for all column of the alignement


        HashMap<Integer, Counter> count = new HashMap<Integer, Counter>();
        int posMin = Integer.MAX_VALUE; // the minimum position in the map count
        int posMax = Integer.MIN_VALUE; // the maximum position in the map count
        // these two variables allows to know if a counter already exists for a specific offset
        int offset = 0;
        List<SequenceAlignment> path = this.path;

        for(SequenceAlignment alignment : path)
        {

            Sequence s1 = alignment.s1;
            Sequence s2 = alignment.s2;

            if(s1.content[0] == Nucleotide.GAP )
            {
                //the offset has to be moved to the left
                int i = 0;

                while(s1.content[i] == Nucleotide.GAP)
                {
                    offset -= 1;
                    i += 1;
                }
            }

            int pos = offset;


            for(int i = 0 ; i < s1.getSize(); i++)
            {
                if( ! count.containsKey(pos))
                {
                    //There is no counter yet for this offset
                    Counter c = new Counter();
                    c.vote(s1.content[i]);
                    count.put(pos,c);

                    posMin = Math.min(posMin, pos);
                    posMax = Math.max(posMax, pos);
                }


                if(s1.content[i] == Nucleotide.GAP && count.get(pos).getLast() != Nucleotide.GAP)
                {

                    //There is a gap we have to propagate upward
                    //All counter after this position are shift to the right

                    Counter c = new Counter();
                    c.vote(s2.content[i]);
                    count.put(posMax+1,count.get(posMax));
                    for(int j = posMax; j > pos ; j--)
                        count.put(j+1,count.get(j));

                    count.put(pos,c);
                    pos +=1;
                    posMax +=1;
                }
                else
                {
                    if (s1.content[i] != Nucleotide.GAP && count.get(pos).getLast() == Nucleotide.GAP)
                    {

                        //There is a gap in the down alignement at the precedent step we have to propagate


                        count.get(pos).vote(Nucleotide.GAP);
                        count.get(pos).setMoveRight(1); // stock into the counter the shift to do in this position
                        //due to a down propagation

                        if(count.containsKey(pos+1))
                            count.get(pos + 1).vote(s2.content[i]);
                        else
                        {
                            Counter c = new Counter();
                            c.vote(s2.content[i]);
                            count.put(pos +1,c);
                            posMax += 1;

                        }

                        pos += 2;

                    }
                    else
                    {

                        //There exists maybe a propagation induces by a gap propopagation (downward)

                        int moveRight = count.get(pos).getMoveRight();
                        count.get(pos + moveRight).vote(s2.content[i]);
                        pos += 1 + moveRight;
                    }
                }

            }
            if (s1.content[0] !=Nucleotide.GAP)
            {
                int i = 0;
                while(s2.content[i] == Nucleotide.GAP)
                {
                    offset += 1;
                    i += 1;
                }
            }

        }
        return count;


    }

    public Sequence buildSequence()
    {
        HashMap<Integer,Counter> votes = this.computeConsensus();

        int min = votes.keySet().stream().mapToInt(Integer::intValue).min().getAsInt();
        int max = votes.keySet().stream().mapToInt(Integer::intValue).max().getAsInt();


        byte[] result = new byte[max - min +1];

        int i = min;
        int index_tab = 0;
        while(i <= max)
        {

            if( votes.get(i).getMaximum() != Nucleotide.GAP)
                result[index_tab] = votes.get(i).getMaximum();
            i += 1;
            index_tab += 1;
        }

        return new Sequence(result);

    }
}
