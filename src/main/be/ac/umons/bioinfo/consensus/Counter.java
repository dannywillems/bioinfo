package be.ac.umons.bioinfo.consensus;

import be.ac.umons.bioinfo.sequence.*;
/**
 * Created by aline on 5/05/16.
 */
public class Counter {

    public static final int GAP = 4;
    public static final int C = 0;
    public static final int G = 1;
    public static final int T = 2;
    public static final int A = 3;

    private int last;
    private int moveRight;
    private int nbrGAP;
    private int nbrC;
    private int nbrG;
    private int nbrT;
    private int nbrA;

    public Counter(){
        this.nbrC = 0;
        this.nbrG = 0;
        this.nbrT = 0;
        this.nbrA = 0;
        this.moveRight = 0;
    }

    public int getLast(){
        return this.last;
    }

    public void setMoveRight(int addMoveRight){this.moveRight += addMoveRight;}
    public int getMoveRight(){return this.moveRight;}

    public void vote(int nucleotide){
        switch(nucleotide)
        {
            case GAP: this.nbrGAP += 1; break;
            case C : this.nbrC += 1; break;
            case G : this.nbrG += 1; break;
            case T : this.nbrT += 1; break;
            case A : this.nbrA +=1; break;
        }
        this.last = nucleotide;
    }

    public byte getMaximum()
    {
        int indexMax = -1;
        int max_value = -1;
        //System.out.println("nbr Gap, "+ nbrGAP);
        //System.out.println("nbr G, "+ nbrG);
        //System.out.println("nbr C, "+ nbrC);
        //System.out.println("nbr T, "+ nbrT);
        //System.out.println("nbr A, "+ nbrA);

        int[] tab = new int[]{this.nbrGAP, this.nbrC,this.nbrG, this.nbrT, this.nbrA};

        //If the max is a gap, we consider the second max

        for(int i = 1; i <= 4 ;i++)
        {
            if(tab[i] >= max_value)
            {
                max_value = tab[i];
                indexMax = i;
            }
        }

        if (indexMax == 0)
            return GAP;
        if(indexMax == 1)
            return C;
        if(indexMax == 2)
            return G;
        if(indexMax == 3)
            return T;
        else
            return A;


    }

    public String toString ()
    {
        return "(last: "+last+", "+
                "nbrGAP: "+nbrGAP+", "+
                "nbrG: "+nbrG+", "+
                "nbrC: "+nbrC+","+
                "nbrA: "+nbrA+","+
                "nbrT: "+nbrT+")";

    }


}
