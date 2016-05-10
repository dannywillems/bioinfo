package be.ac.umons.bioinfo.sequence;

import java.util.*;

/**
 * Created by aline on 5/05/16.
 */
public class ConsensusAline {

    public List<SequenceAlignment> path;

    public ConsensusAline(List<SequenceAlignment> path){
        this.path = path;
    }

    public static Sequence consensus(List<SequenceAlignment> path)
    {
        return produceConsensusSequence(buildConsensus(path));
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

    public static Map<Integer, MyCounter> buildConsensus(List<SequenceAlignment> path)
    {
        Map<Integer, MyCounter> votes = new HashMap<Integer, MyCounter>();

        return votes;
    }


    public HashMap<Integer, Counter> computeAlignementDebug(List<SequenceAlignment> path) {

        HashMap<Integer, Counter> count = new HashMap<Integer, Counter>();
        int posMin = Integer.MAX_VALUE; // the minimum position in the map count
        int posMax = Integer.MIN_VALUE; // the maximum position in the map count
        // these two variables allows to know if a counter already exists for a specific offset
        int offset = 0;

        for (SequenceAlignment alignment : path) {


            Sequence s1 = alignment.s1;
            Sequence s2 = alignment.s2;
            int dec = 0;
            if (s1.content[0] == Sequence.GAP) {
                while (s1.content[dec] == Sequence.GAP) {
                    offset -= 1;
                    dec += 1;
                }
            }

            //On est mis à la position la plus à gauche possible, on peut commencer à compter
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
            if (s1.content[0] !=Sequence.GAP)
            {
                int j = 0;
                while(s2.content[j] ==Sequence.GAP)
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

            if(s1.content[0] == Sequence.GAP )
            {
                //the offset has to be moved to the left
                int i = 0;

                while(s1.content[i] == Sequence.GAP)
                {
                    offset -= 1;
                    i += 1;
                }
            }

            int pos = offset;
            //System.out.print("offset "+offset);
            //System.out.println(" pos "+pos);


            for(int i = 0 ; i < s1.getSize(); i++)
            {
                //System.out.println("je suis en train de traiter la position "+pos);
                //if( !(posMin <= pos && pos <= posMax))
                if(!count.containsKey(pos)) //count.containsKey(pos) ?
                {
                    //System.out.println("Pas encore de vote pour cette position");
                    //There is no counter yet for this offset
                    Counter c = new Counter();
                    c.vote(s1.content[i]);
                    //System.out.println("vote a cette etape");
                    //System.out.println(c);
                    count.put(pos,c);

                    posMin = Math.min(posMin, pos);
                    posMax = Math.max(posMax, pos);
                }


                if(s1.content[i] == Sequence.GAP && count.get(pos).getLast() != Sequence.GAP)
                {
                    //System.out.println("propagation de gap vers le haut");

                    // Alors nous sommes dans le cas où nous voulons propager un gap vers le haut
                    //Tous les compteurs enregistrés à une position supérieure doivent donc etre decalés vers la droite

                    Counter c = new Counter();
                    c.vote(s2.content[i]);
                    // il faut décaler toutes les clefs plus grandes que pos d'une unité
                    count.put(posMax+1,count.get(posMax));
                    for(int j = posMax; j > pos ; j--)
                        count.put(j+1,count.get(j));

                    count.put(pos,c);
                    pos +=1;
                    posMax +=1;
                }
                else
                {
                    if (s1.content[i] != Sequence.GAP && count.get(pos).getLast() ==Sequence.GAP)
                    {
                        //System.out.println("propagation de gap vers le bas");


                        //Alors cela signifie qu'il y avait un gap dans la sequence inferieure de l etape precedente
                        // à propager vers le bas.
                        //Il faut donc une nouvelle fois décaler les compteurs vers la droite
                        //Avant cela on va voter puis seulement decaller
                        //Sinon on peut decaller et puis voter pour deux positions plus loin que la position en cours
                        // de traitement

                        count.get(pos).vote(Sequence.GAP);// voter pour de vrai pour unSequence.GAP ou se contenter de dire que
                        // à l'étape d'avant on avait considéré un Sequence.GAP ?
                        count.get(pos).setMoveRight(1);

                        assert(count.containsKey(pos+1)==true);
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
                        //System.out.println("pas de propagation");

                        //mais peut etre qu'on a une propagation induite avant

                        int moveRight = count.get(pos).getMoveRight();
                        if(count.containsKey(pos+moveRight))
                            count.get(pos + moveRight).vote(s2.content[i]);
                        else
                        {
                            Counter c = new Counter();
                            c.vote(s2.content[i]);
                            count.put(pos+moveRight,c);
                        }

                        //System.out.println("actualisiation du vote");
                        //System.out.println(count.get(pos + moveRight));
                        pos += 1 + moveRight;
                    }
                }

            }
            if (s1.content[0] !=Sequence.GAP)
            {
                int i = 0;
                while(s2.content[i] ==Sequence.GAP)
                {
                    offset += 1;
                    i += 1;
                }
            }
            //System.out.println("offset en fin de traitement "+offset);

        }
        return count;


    }

    public Sequence buildSequence()
    {
        HashMap<Integer,Counter> map = this.computeConsensus();

        int min = map.keySet().stream().mapToInt(Integer::intValue).min().getAsInt();
        int max = map.keySet().stream().mapToInt(Integer::intValue).max().getAsInt();

        byte[] result = new byte[max - min +1];

        int i = min;
        int index_tab = 0;
        while(i <= max)
        {
            //System.out.println("position "+ i);
            //System.out.println("maximum "+map.get(i).getMaximum());
            if(map.get(i).getMaximum() != Sequence.GAP)
                result[index_tab] = map.get(i).getMaximum();
            i += 1;
            index_tab += 1;
        }

        return new Sequence(result);

    }
}
