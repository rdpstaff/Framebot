/*
 * Copyright (C) 2012 Michigan State University <rdpstaff at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
