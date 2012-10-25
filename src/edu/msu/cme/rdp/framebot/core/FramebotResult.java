package edu.msu.cme.rdp.framebot.core;


import edu.msu.cme.rdp.alignment.pairwise.rna.DistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.OverlapCheckFailedException;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import java.util.Collections;
import java.util.List;

public class FramebotResult {
    private static final DistanceModel dm = new IdentityDistanceModel(true);

    private Sequence alignedQuery;   // the aligned of the nucl protein
    private Sequence alignedTarget;  // the aligned of the subject protein
    private Sequence correctedNucl; // the corrected query nucl

    private int frameshifts;
    private FrameScore frameScore;
    private List<Integer> frames;  // the reading frame for each AA, relative to the original nucl seq
    private List<Integer> numberOfOrigBasesList;  // the number of nucleotide bases from the original nucl seq for each corresponding AA
    private double percentIdent;

    public FramebotResult(Sequence alignedQuery, Sequence alignedTarget, Sequence correctedNucl, int frameshifts, List<Integer> frames, List<Integer> numberOfOrigBasesList, FrameScore frameScore) {
        this.alignedQuery = alignedQuery;
        this.alignedTarget = alignedTarget;
        this.correctedNucl = correctedNucl;
        this.frameshifts = frameshifts;
        this.frames = Collections.unmodifiableList(frames);
        this.numberOfOrigBasesList = Collections.unmodifiableList(numberOfOrigBasesList);
        this.frameScore = frameScore;

        try {
            percentIdent = 1 - dm.getDistance(alignedQuery.getSeqString().getBytes(), alignedTarget.getSeqString().getBytes(), 1);
        } catch(OverlapCheckFailedException ignore) {}
    }

    public Sequence getAlignedQuery() {
        return alignedQuery;
    }

    public Sequence getAlignedTarget() {
        return alignedTarget;
    }

    public Sequence getCorrectedNucl() {
        return correctedNucl;
    }

    public List<Integer> getFrames() {
        return frames;
    }

    public List<Integer> getNumberOfOrigBasesList() {
        return numberOfOrigBasesList;
    }

    public int getFrameshifts() {
        return frameshifts;
    }

    public FrameScore getFrameScore() {
        return frameScore;
    }

    public double getPercentIdent() {
        return percentIdent;
    }
}
