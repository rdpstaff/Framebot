package edu.msu.cme.rdp.framebot.core;

public class FrameScore {

    private int queryStart = 0;
    private int targetStart = 0;
    private int queryEnd;
    private int targetEnd;
    private int bestFrame;
    private int maxScore;

    public FrameScore(int bestFrame, int maxScore, int maxQueryIndex , int maxSubjectIndex){
        this.bestFrame = bestFrame;
        this.maxScore = maxScore;
        this.targetEnd = maxSubjectIndex;
        this.queryEnd = maxQueryIndex;
    }

    public int getBestFrame() {
        return bestFrame;
    }

    public int getMaxScore() {
        return maxScore;
    }

    public int getQueryStart() {
        return queryStart;
    }

    public int getTargetStart() {
        return targetStart;
    }

    public int getQueryEnd() {
        return queryEnd;
    }

    public int getTargetEnd() {
        return targetEnd;
    }

    void setQueryStart(int queryStart) {
        this.queryStart = queryStart;
    }

    void setTargetStart(int targetStart) {
        this.targetStart = targetStart;
    }
}
