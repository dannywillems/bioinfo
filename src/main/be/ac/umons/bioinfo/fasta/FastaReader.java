package be.ac.umons.bioinfo.fasta;

/**
 * Created by aline on 31/03/16.
 */

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import be.ac.umons.bioinfo.sequence.Sequence;

/**
 * This class allows reading of a FASTA file.
 */
public class FastaReader
{
    /**
     * Reads a text file and generates a sequence of DNA sequences based on this file,
     * supposing the file use FASTA as data structure.
     * @param f a FASTA file
     * @return a sequence of DNA sequences read from a FASTA file.
     */
    public static List<Sequence> readFromFile(File f) throws IOException
    {
        BufferedReader br = new BufferedReader(new FileReader(f));

        StringBuilder builder = new StringBuilder();
        List<Sequence> ret = new ArrayList<Sequence>();
        String line;

        while((line = br.readLine()) != null)
        {
            if(line.startsWith(">"))
            {
                String text = builder.toString();
                if(!text.isEmpty())
                    ret.add(new Sequence(text));

                builder = new StringBuilder();
            }
            else
                builder.append(line.trim());
        }

        return ret;
    }
}
