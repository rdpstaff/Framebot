package edu.msu.cme.rdp.framebot.cli;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;

public class Translation {

    private static ProteinUtils proteinUtils = ProteinUtils.getInstance();

    public static void parse(int translTable, String nuclFile, String protFile, int frame) throws IOException {

        SequenceReader reader = new SequenceReader(new File(nuclFile));
        BufferedWriter writer = new BufferedWriter(new FileWriter(new File(protFile)));
        Sequence seq = null;
        while ((seq = (Sequence) (reader.readNextSequence())) != null) {
            String s = translate(seq.getSeqString().replaceAll("-", ""), true, translTable, frame);

            writer.write(">" + seq.getSeqName() + "\n" + s + "\n");
        }
        reader.close();
        writer.close();
    }

    public static String translate(String nucl, boolean dontAllowInitiators, int translTable, int frame) {
        String s = proteinUtils.translateToProtein(nucl.substring(frame - 1), dontAllowInitiators, translTable).toUpperCase();
        return s;
    }

    public static void main(String[] args) throws IOException {

        String usage = "Translation translTable nuclfile outProtfile frame \nframe start from 1-3";
        if (args.length != 4) {
            throw new IllegalArgumentException(usage);
        }
        Translation.parse(Integer.parseInt(args[0]), args[1], args[2], Integer.parseInt(args[3]));

    }
}
