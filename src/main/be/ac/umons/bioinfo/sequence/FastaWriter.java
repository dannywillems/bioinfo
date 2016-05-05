package be.ac.umons.bioinfo.sequence;

/**
 * Created by aline on 31/03/16.
 */
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

public class FastaWriter
{
    /**
     * Writes a sequence of DNA sequences preceded by a free text.
     * @param sequences The DNA sequences to write in a FASTA file, with their associated texts.
     * @param f The file in which the DNA sequences must be writen, in the FASTA format.
     * @param chunkSize Max size of a chunk. Typically 80.
     */
    public static void write(Map<String, Sequence> sequences, File f, int chunkSize) throws IOException
    {
        FileWriter fw = new FileWriter(f);
        BufferedWriter bw = new BufferedWriter(fw);

        for(Map.Entry<String, Sequence> sequence : sequences.entrySet())
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
}
