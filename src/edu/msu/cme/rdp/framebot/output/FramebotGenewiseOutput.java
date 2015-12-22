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
package edu.msu.cme.rdp.framebot.output;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.framebot.core.FramebotCore;
import edu.msu.cme.rdp.framebot.core.FramebotResult;
import edu.msu.cme.rdp.framebot.core.FramebotScoringMatrix;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import org.apache.commons.lang.StringUtils;
import java.util.Locale;

/**
 *
 * @author fishjord
 */
public class FramebotGenewiseOutput {

    private static final int printWidth = 60;
    private static final DecimalFormat format = (DecimalFormat)DecimalFormat.getInstance(Locale.ENGLISH);
    static {
        format.applyPattern("0.000");
    }

    private PrintStream out;

    public FramebotGenewiseOutput(File out) throws IOException {
        this(new BufferedOutputStream(new FileOutputStream(out)));
    }

    public FramebotGenewiseOutput(OutputStream out) throws IOException {
        this.out = new PrintStream(out);
    }

    public void printResult(FramebotResult result, boolean reversed) {

        char[] targetAligned = result.getAlignedTarget().getSeqString().toCharArray();
        char[] queryAligned = result.getAlignedQuery().getSeqString().toCharArray();
        char[] correctedNucl = result.getCorrectedNucl().getSeqString().toCharArray();

        StringBuilder match = new StringBuilder();
        StringBuilder frame1 = new StringBuilder();
        StringBuilder frame2 = new StringBuilder();
        StringBuilder frame3 = new StringBuilder();
        StringBuilder frames = new StringBuilder();


        int lastFrame = result.getFrames().get(0);
        int nuclIndex = 0;
        for (int protAlignIndex = 0; protAlignIndex < targetAligned.length; protAlignIndex++) {
            char targetAA = targetAligned[protAlignIndex];
            char queryAA = queryAligned[protAlignIndex];
            int currFrame = result.getFrames().get(protAlignIndex);

            //We wouldn't want a match bar if both are gaps, but since this is a pairwise alignment
            //it shouldn't be the case that both AAs are gaps...that'd be weird right?
            if (targetAA == queryAA) {
                match.append("|");
            } else {
                match.append(" ");
            }

            if ( result.getNumberOfOrigBasesList().get(protAlignIndex) == FramebotScoringMatrix.NUM_FRAMES) {

                frame1.append(correctedNucl[nuclIndex]);
                frame2.append(correctedNucl[nuclIndex + 1]);
                frame3.append(correctedNucl[nuclIndex + 2]);
               
                nuclIndex += 3;
            }else {
                if (  result.getNumberOfOrigBasesList().get(protAlignIndex) == 0){
                    frame1.append(' ');
                }else {
                 // there is a frameshift, we need to print the number of original nucleotides in that AA spot
                    frame1.append(result.getNumberOfOrigBasesList().get(protAlignIndex));
                }
                frame2.append(' ');
                frame3.append(' ');
                if (queryAA != '-') {
                    nuclIndex += 3;
                }
            }

            frames.append(result.getFrames().get(protAlignIndex));
            lastFrame = currFrame;
        }
        out.println(">\tTarget\tQuery\tNuclLen\tAlignLen\t%Identity\tScore\tFrameshifts\tReversed");
        out.println("STATS\t" +
                result.getAlignedTarget().getSeqName() + "\t" +
                result.getAlignedQuery().getSeqName() + "\t" +
                result.getCorrectedNucl().getSeqString().length() + "\t" +
                result.getAlignedQuery().getSeqString().length() + "\t" +
                format.format(result.getPercentIdent()*100) + "\t" +
                result.getFrameScore().getMaxScore() + "\t" +
                result.getFrameshifts() + "\t" +
                reversed);

        for(int index = 0;index < targetAligned.length;index += printWidth) {
            int width = printWidth;
            if(index + width > targetAligned.length) {
                width = targetAligned.length - index;
            }
            String targetStringPart = new String(targetAligned, index, width);
            String queryStringPart = new String(queryAligned, index, width);
            String matchPart = match.substring(index, index + width);
            String framesPart = frames.substring(index, index + width);
            String frame1Part = frame1.substring(index, index + width);
            String frame2Part = frame2.substring(index, index + width);
            String frame3Part = frame3.substring(index, index + width);

            int targetIndex = result.getFrameScore().getTargetStart() + index;
            int queryIndex = result.getFrameScore().getQueryStart() + index;

            out.println("Target" + StringUtils.leftPad(targetIndex + 1 + "", 5) + " " + targetStringPart + " " + StringUtils.leftPad((targetIndex + width) + "", 5));
            out.println("            " + matchPart);
            out.println("Query " + StringUtils.leftPad(queryIndex + 1 + "", 5) + " " + queryStringPart + " " +  StringUtils.leftPad((queryIndex + width) + "", 5 ));
            out.println("Frame       " + framesPart);
            out.println("            " + frame1Part);
            out.println("            " + frame2Part);
            out.println("            " + frame3Part);
            out.println();
        }

    }

    public void close() {
        out.close();
    }

    public static void main(String[] args) throws Exception {
        Sequence seedSeq = new Sequence("52080945", "", "mqvqekrilvinpgststkigvfhddrsifeksirhdeaelqqyqtiidqysfrkqailetlheqginiskldavcarggllrpieggtyevndamivdlkngyagqhasnlggiiareiadglnipafivdpvvvdemapiakisgtpaierrsifhalnqkavarkaawqfgkryedmkmiithmgggitigvhcrgrvidvnnglhgegplsperagtipagdlidmcfsgeytkdelmkmlvgggglagylgttdavkvekmikegdqkaaliyeamayqiakeigaasavlkgevdviiltgglaygksfissirqyidwisdvvvfpgenelqalaegafrvlngeeeakqypnqrreshgn");
        Sequence queryNucl = new Sequence("GF040U108JKLLA", "", "aaagatcggcgtttttcatgatgaccgttcgattttcgaaaaatcaatccgtcatgacgaggctgagctacagcaatatcagaccattattgatcaatattcgttcagaaaacaggcgatactcgaaaccctgcatgaacagggaatcaatatttctaaattggatgccgtttgcgccaggggagggctgcttcggccgattgaaggcggcacttacgaagtcaatgatgcgatgattgtcgatttgaaaaacggctatgcggggcagcatgcatcaaatctcgggggcatcatcgccagggagattgccgacgggttaaatattcccgctttttatcgtcgaccccgttgttgtggatgaatggctcctatcgcaaaatttccggcaccccggctattgaaaggcgcagcatttttcac");

        FramebotGenewiseOutput out = new FramebotGenewiseOutput(System.out);
        
        FramebotResult result = FramebotCore.processSequence(queryNucl, seedSeq, true, 11, AlignmentMode.glocal, ScoringMatrix.getDefaultProteinMatrix());
        out.printResult(result, false);

        seedSeq = new Sequence("52080945", "", "mqvqekrilvinpgststkigvfhddrsifeksirhdeaelqqyqtiidqysfrkqailetlheqginiskldavcarggllrpieggtyevndamivdlkngyagqhasnlggiiareiadglnipafivdpvvvdemapiakisgtpaierrsifhalnqkavarkaawqfgkryedmkmiithmgggitigvhcrgrvidvnnglhgegplsperagtipagdlidmcfsgeytkdelmkmlvgggglagylgttdavkvekmikegdqkaaliyeamayqiakeigaasavlkgevdviiltgglaygksfissirqyidwisdvvvfpgenelqalaegafrvlngeeeakqypnqrreshgn");
        queryNucl = new Sequence("GF040U108JMSLD", "", "aaagatcggcgtttttcatgatgaccgttcgattttcgaaaaatcaatccgtcatgacgaggctgagctacagcaatatcagaccattattgatcaatattcgttcagaaaacaggcgatactcgaaaccctgcatgaacagggaatcaatatttctaaattggatgccgtttgcgccaggggagggctgcttcggccgattgaaggcggcacttacgaagtcaatgatgcgatgattgtcgatttgaaaaacggctatgcggggcagcatgcatcaaatctcgggggcatcatcgccagggagattgccgacgggttaaatattcccgcttttatcgtcgaccccgttgttgtggatgaaatggctcctatcgcaaaatttccggcaccccggctattgaaggcgcagcatttttcac");
        result = FramebotCore.processSequence(queryNucl, seedSeq, true, 11, AlignmentMode.glocal, ScoringMatrix.getDefaultProteinMatrix());
        out.printResult(result, false);
    }
}
