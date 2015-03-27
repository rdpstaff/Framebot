/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.output;

import edu.msu.cme.rdp.framebot.core.FramebotResult;
import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author fishjord
 */
public class OutputCoordinator {

    private FramebotGenewiseOutput framebotOut;
    private FramebotGenewiseOutput failedFramebotOutput;
    private FastaWriter correctedProtOut;
    private FastaWriter correctedNuclOut;
    private FastaWriter correctedNuclQualOut;
    private File qualOutFile;
    private FastaWriter nuclFailedOut;
    private int lengthCutoff;
    private double identityCutoff;

    public OutputCoordinator(int lengthCutoff, double identityCutoff, File framebotOutputFile, File failedFramebotOutputFile, File correctedProtOutputFile, File correctedNuclOutputFile, File failedNuclOutputFile) throws IOException {
        this(lengthCutoff, identityCutoff, framebotOutputFile, failedFramebotOutputFile, correctedProtOutputFile, correctedNuclOutputFile, null, failedNuclOutputFile);
    }

    public double getIdentityCutoff(){
        return identityCutoff;
    }
    
    public OutputCoordinator(int lengthCutoff, double identityCutoff, File framebotOutputFile, File failedFramebotOutputFile, File correctedProtOutputFile, File correctedNuclOutputFile, File correctedNuclQualOut, File failedNuclOutputFile) throws IOException {
        this.lengthCutoff = lengthCutoff;
        this.identityCutoff = identityCutoff;

        framebotOut = new FramebotGenewiseOutput(framebotOutputFile);
        failedFramebotOutput = new FramebotGenewiseOutput(failedFramebotOutputFile);
        correctedProtOut = new FastaWriter(correctedProtOutputFile);
        correctedNuclOut = new FastaWriter(correctedNuclOutputFile);
        nuclFailedOut = new FastaWriter(failedNuclOutputFile);

        if(correctedNuclQualOut != null) {
            this.correctedNuclQualOut = new FastaWriter(correctedNuclQualOut);
            qualOutFile = correctedNuclQualOut;
        }
    }

    public void printResult(FramebotResult result, Sequence origNucl, boolean reversed) {

        Sequence unalignedProt = new Sequence(result.getAlignedQuery().getSeqName(), "", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()));

        if(result.getPercentIdent() < identityCutoff || unalignedProt.getSeqString().length() < lengthCutoff) {
            failedFramebotOutput.printResult(result, reversed);
            nuclFailedOut.writeSeq(origNucl);
        } else {
            framebotOut.printResult(result, reversed);
            correctedProtOut.writeSeq(unalignedProt);
            Sequence correctedNuclSeq = result.getCorrectedNucl();
            correctedNuclOut.writeSeq(correctedNuclSeq);

            if((correctedNuclSeq instanceof QSequence) && correctedNuclQualOut != null) {
                correctedNuclQualOut.writeRawQual(correctedNuclSeq.getSeqName(), ((QSequence)correctedNuclSeq).getQuality());
            }
        }
    }

    public void close() {
        framebotOut.close();
        failedFramebotOutput.close();
        correctedProtOut.close();
        correctedNuclOut.close();
        nuclFailedOut.close();

        if(qualOutFile != null && qualOutFile.length() == 0) {
            qualOutFile.delete();
        }
    }
}
