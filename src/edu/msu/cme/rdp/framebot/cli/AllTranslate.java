/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.framebot.cli;

import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author fishjord
 */
public class AllTranslate {

    public static void main(String[] args) throws IOException {
        if(args.length != 2) {
            System.err.println("USAGE: AllTranslate <in_nucl_file> <out_file>");
            System.exit(1);
        }

        SeqReader reader = new SequenceReader(new File(args[0]));
        FastaWriter writer = new FastaWriter(new PrintStream(args[1]));

        Sequence seq;
        ProteinUtils protUtils = ProteinUtils.getInstance();

        while((seq = reader.readNextSequence()) != null) {
            for(Sequence protSeq : protUtils.allTranslate(seq)) {
                writer.writeSeq(protSeq);
            }
        }

        writer.close();
    }

}
