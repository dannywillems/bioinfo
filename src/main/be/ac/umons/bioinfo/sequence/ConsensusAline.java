package be.ac.umons.bioinfo.sequence;

import be.ac.umons.bioinfo.Pair;

import java.util.*;

/**
 * Created by aline on 5/05/16.
 */
public class ConsensusAline {

    public List<SequenceAlignment> path;
    public static final int GAP = 4;
    public static final int C = 0;
    public static final int G = 1;
    public static final int T = 2;
    public static final int A = 3;

    public ConsensusAline(List<SequenceAlignment> path){
        this.path = path;
    }

    public ConsensusResult computeConsensus(){

        //stock the counter for all column of the alignement


        HashMap<Integer, Counter> count = new HashMap<Integer, Counter>();
        int posMin = Integer.MAX_VALUE; // the minimum position in the map count
        int posMax = Integer.MIN_VALUE; // the maximum position in the map count
        // these two variables permit to know if a counter already exists for a specific offset
        int offset = 0;
        List<SequenceAlignment> path = this.path;

        for(SequenceAlignment alignment : path)
        {

            Sequence s1 = alignment.s1;
            Sequence s2 = alignment.s2;

            if(s1.content[0] == GAP )
            {
                //the offset has to be moved to the left
                int i = 0;

                while(s1.content[i] == GAP)
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
                if( !(posMin <= pos && pos <= posMax)) //count.containsKey(offset) ?
                {
                    //System.out.println("Pas encore de vote pour cette position");
                    //There is no counter yet for this offset
                    Counter c = new Counter();
                    c.refresh(s1.content[i]);
                    //System.out.println("vote a cette etape");
                    //System.out.println(c);
                    count.put(pos,c);

                    if (pos < posMin)
                        posMin = pos;
                    else
                        posMax = pos;

                }

                if(s1.content[i] == GAP && count.get(pos).getLast() != GAP)
                {
                    //System.out.println("propagation de gap vers le haut");

                    // Alors nous sommes dans le cas où nous voulons propager un gap vers le haut
                    //Tous les compteurs enregistrés à une position supérieure doivent donc etre decalles vers la droite

                    Counter c = new Counter();
                    c.refresh(s2.content[i]);
                    // il faut decaller toutes les clefs plus grandes que pos d'une unité
                    count.put(posMax+1,count.get(posMax));
                    for(int j = posMax; j > pos ; j--)
                    {
                        count.replace(j,count.get(j-1));
                    }
                    count.put(pos,c);
                    pos +=1;
                    posMax +=1;
                }

                else {
                    if (s1.content[i] != GAP && count.get(pos).getLast() == GAP) {
                        //System.out.println("propagation de gap vers le bas");


                        //Alors cela signifie qu'il y avait un gap dans la sequence inferieure de l etape precedente
                        // à propager vers le bas.
                        //Il faut donc une nouvelle fois décaller les compteurs vers la droite
                        //Avant cela on va voter puis seulement decaller
                        //Sinon on peut decaller et puis voter pour deux positions plus loin que la position en cours
                        // de traitement

                        count.get(pos).refresh(GAP);// voter pour de vrai pour un GAP ou se contenter de dire que
                        // à l'étape d'avant on avait considéré un GAP ?
                        count.get(pos).setMoveRight(1);

                        if(count.containsKey(pos+1))
                            count.get(pos + 1).refresh(s2.content[i]);
                        else
                        {
                            Counter c = new Counter();
                            c.refresh(s2.content[i]);
                            count.put(pos +1,c);
                        }

                        pos += 2;

                    }
                    else
                    {
                        //System.out.println("pas de propagation");

                        //mais peut etre qu'on a une propagation induite avant

                        int moveRight = count.get(pos).getMoveRight();
                        count.get(pos + moveRight).refresh(s2.content[i]);

                        //System.out.println("actualisiation du vote");
                        //System.out.println(count.get(pos + moveRight));
                        pos += 1 + moveRight;
                    }
                }

            }
            if (s1.content[0] != GAP)
            {
                int i = 0;
                while(s2.content[i] == GAP)
                {
                    offset += 1;
                    i += 1;
                }
            }
            //System.out.println("offset en fin de traitement "+offset);

        }
        return new ConsensusResult(count,posMin,posMax);


    }

    public Sequence buildSequence()
    {
        ConsensusResult consensus = this.computeConsensus();
        HashMap<Integer,Counter> map = consensus.map;

        byte[] result = new byte[consensus.posMax - consensus.posMin +1];

        int i = consensus.posMin;
        int index_tab = 0;
        while(i <= consensus.posMax)
        {
            //System.out.println("position "+ i);
            //System.out.println("maximum "+map.get(i).getMaximum());
            result[index_tab] = map.get(i).getMaximum();
            i += 1;
            index_tab += 1;
        }

        return new Sequence(result);

    }
}
