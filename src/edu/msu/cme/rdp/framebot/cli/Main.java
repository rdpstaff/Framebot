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

import edu.msu.cme.rdp.framebot.index.FramebotIndex;
import edu.msu.cme.rdp.framebot.stat.GetFrameBotStatMain;
import edu.msu.cme.rdp.framebot.stat.RdmSelectSampleMapping;
import edu.msu.cme.rdp.framebot.stat.TaxonAbundance;
import java.util.Arrays;

/**
 *
 * @author fishjord
 */
public class Main {

    public static void main(String [] args) throws Exception {
        String usage = "USAGE: Main <subcommand> <subcommand args ...>" +
                "\n\tframebot      - run framebot" +
                "\n\tindex         - build an index" +
                "\n\tstat          - convert framebot output files to different output formats" +
                "\n\ttaxonAbund    - taxonomic abundance group by phylum or match" +
                "\n\trdmselect     - randomly selects a subset of sequence IDs from the sample Mapping file" +
                "\n\ttranslate     - translate nucleotide sequences to protein (no frameshift correction) " +
                "\n\tall-translate - translate nucleotide sequences to protein, all reading frames (no frameshift correction)";
        if(args.length == 0 ) {
            System.err.println(usage);
            return;
        }

        String cmd = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);

        if(cmd.equals("framebot")) {
            FramebotMain.main(newArgs);
        } else if(cmd.equals("index")) {
            FramebotIndex.main(newArgs);
        } else if(cmd.equals("translate")) {
            Translation.main(newArgs);
        } else if(cmd.equals("all-translate")) {
            AllTranslate.main(newArgs);
        } else if(cmd.equals("stat")) {
            GetFrameBotStatMain.main(newArgs);
        } else if(cmd.equals("rdmselect")) {
            RdmSelectSampleMapping.main(newArgs);
        } else if(cmd.equals("taxonAbund")) {
            TaxonAbundance.main(newArgs);
        }else {
            System.err.println("ERROR: " + "wrong subcommand");
            System.err.println(usage);
            return;
        }
    }

}
