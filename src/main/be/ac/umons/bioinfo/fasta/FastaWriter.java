package be.ac.umons.bioinfo.fasta;

/**
 * Created by aline on 31/03/16.
 */
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import java.util.Map;
import java.util.ArrayList;

import be.ac.umons.bioinfo.sequence.SequenceAbstract;
import be.ac.umons.bioinfo.sequence.SequenceAlignment;
import be.ac.umons.bioinfo.sequence.SequenceAlignmentAbstract;
import be.ac.umons.bioinfo.sequence.Sequence;

public class FastaWriter
{
    /**
     * Writes a sequence of DNA alignment preceded by a free text.
     * @param alignment The DNA alignment to write in a FASTA file, with their associated texts.
     * @param f The file in which the DNA alignment must be writen, in the FASTA format.
     * @param chunkSize Max size of a chunk. Typically 80.
     */
    public static void write(Map<String, Sequence> alignment, File f, int chunkSize) throws IOException
    {
        FileWriter fw = new FileWriter(f);
        BufferedWriter bw = new BufferedWriter(fw);

        for(Map.Entry<String, Sequence> sequence : alignment.entrySet())
        {
            bw.write(">" + sequence.getKey() + "\n");
            int offset = 0;
            String repr = sequence.getValue().toString();

            while(offset < repr.length() - 1)
            {
                int end = Math.min(offset + chunkSize, repr.length());

                CharSequence sub = repr.substring(offset, end);
                bw.write(sub.toString() + "\n");

                offset = end;
            }
        }

        bw.close();
    }

    public static void write(String sequence_name, Sequence sequence, File f, int chunkSize) throws IOException
    {
        FileWriter fw = new FileWriter(f);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write(">" + sequence_name + "\n");
        int offset = 0;
        String repr = sequence.toString();

        while(offset < repr.length() - 1)
        {
            int end = Math.min(offset + chunkSize, repr.length());

            CharSequence sub = repr.substring(offset, end);
            bw.write(sub.toString() + "\n");

            offset = end;
        }

        bw.close();
    }

    public static void write(String sequence_name, SequenceAbstract sequence, File f, int chunkSize) throws IOException
    {
        FileWriter fw = new FileWriter(f);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write(">" + sequence_name + "\n");
        int offset = 0;
        String repr = sequence.toString();

        while(offset < repr.length() - 1)
        {
            int end = Math.min(offset + chunkSize, repr.length());

            CharSequence sub = repr.substring(offset, end);
            bw.write(sub.toString() + "\n");

            offset = end;
        }

        bw.close();
    }

    public static void write(ArrayList<SequenceAbstract> alignment, File f, int chunkSize) throws IOException
    {
        FileWriter fw = new FileWriter(f);
        BufferedWriter bw = new BufferedWriter(fw);

        for (int i = 0;i < alignment.size();i++)
        {
            bw.write(">Sequence " + i + "\n");
            int offset = 0;
            String repr = alignment.get(i).toString();

            while(offset < repr.length() - 1)
            {
                int end = Math.min(offset + chunkSize, repr.length());

                CharSequence sub = repr.substring(offset, end);
                bw.write(sub.toString() + "\n");

                offset = end;
            }
        }
        bw.close();
}

    public static void writeAlignment(ArrayList<SequenceAlignmentAbstract> alignment, File f, int chunkSize) throws IOException
    {
        FileWriter fw = new FileWriter(f);
        BufferedWriter bw = new BufferedWriter(fw);

        for (int i = 0;i < alignment.size();i++)
        {
            bw.write(">Alignement " + i + ": up\n");
            int offset = 0;
            String repr = alignment.get(i).s.toString();

            while(offset < repr.length() - 1)
            {
                int end = Math.min(offset + chunkSize, repr.length());

                CharSequence sub = repr.substring(offset, end);
                bw.write(sub.toString() + "\n");

                offset = end;
            }

            bw.write(">Alignement " + i + ": down\n");
            offset = 0;
            repr = alignment.get(i).t.toString();

            while(offset < repr.length() - 1)
            {
                int end = Math.min(offset + chunkSize, repr.length());

                CharSequence sub = repr.substring(offset, end);
                bw.write(sub.toString() + "\n");

                offset = end;
            }

        }
        bw.close();
    }
}
