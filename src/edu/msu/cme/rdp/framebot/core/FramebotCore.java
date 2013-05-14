/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.core;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.readseq.readers.Sequence;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author wangqion
 */
public class FramebotCore {

    private static final char GAP = '-';
    private static final char UNKNOWN_AA = 'X';
    private static final char UNKNOWN_NUCL = 'n';
    private static final byte UNKNOWN_NUCL_QSCORE = 40;
   
            
    private static final ProteinUtils proteinUtils = ProteinUtils.getInstance();
    private static final int MATCH = 0;
    private static final int DOWN = 1;
    private static final int ACROSS = 2;
    private static final int NEW_START = 3;

    private static class FramebotSeqs {

        char[] queryNuclSeqString;
        char[][] queryProtFrames;
        char[] targetProtSeqString;
        Sequence targetProtSeq;
        Sequence queryNuclSeq;
    }

    public static FramebotResult processSequence(Sequence nucl, Sequence prot, boolean dontAllowInitiators, int translTable, AlignmentMode mode, ScoringMatrix simMatrix) {
               
        // initialize the matrix
        FramebotSeqs seqs = new FramebotSeqs();

        seqs.queryNuclSeq = nucl;
        seqs.targetProtSeq = prot;
        seqs.queryNuclSeqString = nucl.getSeqString().toCharArray();
        seqs.queryProtFrames = new char[FramebotScoringMatrix.NUM_FRAMES][];
        seqs.targetProtSeqString = (" " + prot.getSeqString().toUpperCase()).toCharArray();


        for (int i = 0; i < FramebotScoringMatrix.NUM_FRAMES; i++) {
            seqs.queryProtFrames[i] = (" " + proteinUtils.translateToProtein(seqs.queryNuclSeq.getSeqString().substring(i), dontAllowInitiators, translTable).toUpperCase()).toCharArray();
        }

        int queryLength = seqs.queryProtFrames[0].length;
        int subjectLength = seqs.targetProtSeqString.length;

        //This method will also check the alignment mode to make sure it is supported
        FramebotScoringMatrix scoringMatrices = new FramebotScoringMatrix(queryLength, subjectLength, mode, simMatrix);

        computeMatrix(queryLength, subjectLength, seqs.queryProtFrames, seqs.targetProtSeqString, scoringMatrices, mode, simMatrix);

        return traceback(seqs, scoringMatrices, mode);
    }

    private static void computeMatrix(int queryLength, int subjectLength, char[][] queryProtFrames, char[] targetProtSeqString, FramebotScoringMatrix scoringMatrices, AlignmentMode mode, ScoringMatrix simMatrix) {

        for (int i = 1; i < queryLength; ++i) {
            for (int j = 1; j < subjectLength; ++j) {
                for (int frameIndex = 0; frameIndex < FramebotScoringMatrix.NUM_FRAMES; frameIndex++) {
                    if ((queryProtFrames[frameIndex].length - 1) < i) {
                        continue;
                    }

                    int matchScore = simMatrix.score(queryProtFrames[frameIndex][i], targetProtSeqString[j]);
                    scoringMatrices.p[frameIndex][i][j] = Math.max(scoringMatrices.d[frameIndex][i - 1][j] + simMatrix.getIndelPenalty(), scoringMatrices.p[frameIndex][i - 1][j] + simMatrix.getGapExtend() );
                    scoringMatrices.q[frameIndex][i][j] = Math.max(scoringMatrices.d[frameIndex][i][j - 1] + simMatrix.getIndelPenalty(), scoringMatrices.q[frameIndex][i][j - 1] + simMatrix.getGapExtend());
                    int temp_match = scoringMatrices.d[frameIndex][i - 1][j - 1] + matchScore;
                    if (mode == AlignmentMode.glocal || mode == AlignmentMode.global) {
                        scoringMatrices.d[frameIndex][i][j] = Math.max(temp_match, Math.max(scoringMatrices.p[frameIndex][i][j], scoringMatrices.q[frameIndex][i][j]));
                    } else {
                        scoringMatrices.d[frameIndex][i][j] = Math.max(temp_match, Math.max(0, Math.max(scoringMatrices.p[frameIndex][i][j], scoringMatrices.q[frameIndex][i][j])));
                    }
                    scoringMatrices.match_frame[frameIndex][i][j] = frameIndex;
                    // calculate the path
                    if (mode == AlignmentMode.local && scoringMatrices.d[frameIndex][i][j] == 0) {
                        scoringMatrices.path[frameIndex][i][j] = NEW_START;
                    } else if (scoringMatrices.d[frameIndex][i][j] == temp_match) {
                        scoringMatrices.path[frameIndex][i][j] = MATCH;
                    } else if (scoringMatrices.d[frameIndex][i][j] == scoringMatrices.p[frameIndex][i][j]) {
                        scoringMatrices.path[frameIndex][i][j] = DOWN;
                    } else {
                        scoringMatrices.path[frameIndex][i][j] = ACROSS;
                    }

                    for (int index = 1; index < FramebotScoringMatrix.NUM_FRAMES; index++) {
                        int nextFrame = (frameIndex + index) % FramebotScoringMatrix.NUM_FRAMES;

                        /** if we assume Frame index starting from 1. To avoid reuse the nucleotides from overlapping frames, we need to find the correct offset
                         * 1 -> 2, 1 -> 3, and 2 -> 3, offset should be 2
                         * 2 ->1, 3 -> 1 and 3 -> 2, offset should be 1 */
                        int prevCellOffSet = (frameIndex < nextFrame) ? 2 : 1;
                        if ((i - prevCellOffSet) < 0) {
                            continue;
                        }
                        int t1 = Math.max(scoringMatrices.d[nextFrame][i - prevCellOffSet][j] + simMatrix.getIndelPenalty(), scoringMatrices.p[nextFrame][i - prevCellOffSet][j] + simMatrix.getGapExtend());
                        int t2 = Math.max(scoringMatrices.d[nextFrame][i - prevCellOffSet + 1][j - 1] + simMatrix.getIndelPenalty(), scoringMatrices.q[nextFrame][i - prevCellOffSet + 1][j - 1] + simMatrix.getGapExtend());
                        int temp_ins_newframe_match = scoringMatrices.d[nextFrame][i - prevCellOffSet][j - 1] + matchScore;

                        int t3 = 0;
                        if (mode == AlignmentMode.glocal || mode == AlignmentMode.global) {
                            t3 = Math.max(temp_ins_newframe_match, Math.max(t1, t2));
                        } else {
                            t3 = Math.max(temp_ins_newframe_match, Math.max(0, Math.max(t1, t2)));
                        }

                        if ((t3 + simMatrix.getFrameshiftPenalty() ) > scoringMatrices.d[frameIndex][i][j]) {
                            // update the scoringMatrices.path
                            if (mode == AlignmentMode.local && t3 == 0) {
                                scoringMatrices.path[frameIndex][i][j] = NEW_START;
                            } else if (t3 == temp_ins_newframe_match) {
                                scoringMatrices.path[frameIndex][i][j] = MATCH;
                            } else if (t3 == t1) {
                                scoringMatrices.path[frameIndex][i][j] = DOWN;
                            } else {
                                scoringMatrices.path[frameIndex][i][j] = ACROSS;
                            }
                            scoringMatrices.p[frameIndex][i][j] = t1 + simMatrix.getFrameshiftPenalty();
                            scoringMatrices.q[frameIndex][i][j] = t2 + simMatrix.getFrameshiftPenalty();
                            scoringMatrices.d[frameIndex][i][j] = t3 + simMatrix.getFrameshiftPenalty();
                            scoringMatrices.match_frame[frameIndex][i][j] = nextFrame;
                        }
                    } 
                }
            }
        }
    }

    private static FrameScore getBestFrameScore(FramebotSeqs seqs, FramebotScoringMatrix scoringMatrices, AlignmentMode mode) {
        int bestFrame = 0;
        int maxScore = Integer.MIN_VALUE;
        int maxQueryIndex = 0;
        int maxSubjectIndex = 0;

        switch (mode) {
            case local: {
                for (int frame = 0; frame < FramebotScoringMatrix.NUM_FRAMES; frame++) {
                    for (int queryIndex = 0; queryIndex < seqs.queryProtFrames[frame].length; queryIndex++) {
                        for (int subjectIndex = 0; subjectIndex < seqs.targetProtSeqString.length; subjectIndex++) {
                            if (scoringMatrices.d[frame][queryIndex][subjectIndex] >= maxScore) {
                                maxScore = scoringMatrices.d[frame][queryIndex][subjectIndex];
                                maxQueryIndex = queryIndex;
                                maxSubjectIndex = subjectIndex;
                                bestFrame = frame;
                            }
                        }
                    }
                }
                break;
            }
            case global: {
                for (int frame = 0; frame < FramebotScoringMatrix.NUM_FRAMES; frame++) {
                    int queryIndex = seqs.queryProtFrames[frame].length - 1;
                    int subjectIndex = seqs.targetProtSeqString.length - 1;
                    if (scoringMatrices.d[frame][queryIndex][subjectIndex] >= maxScore) {
                        maxScore = scoringMatrices.d[frame][queryIndex][subjectIndex];
                        maxQueryIndex = queryIndex;
                        maxSubjectIndex = subjectIndex;
                        bestFrame = frame;
                    }
                }
                break;
            }
            case glocal: {
                // local on the subject, global on query
                for (int frame = 0; frame < FramebotScoringMatrix.NUM_FRAMES; frame++) {
                    int queryIndex = seqs.queryProtFrames[frame].length - 1;
                    for (int subjectIndex = 0; subjectIndex < seqs.targetProtSeqString.length; subjectIndex++) {
                        if (scoringMatrices.d[frame][queryIndex][subjectIndex] >= maxScore) {
                            maxScore = scoringMatrices.d[frame][queryIndex][subjectIndex];
                            maxQueryIndex = queryIndex;
                            maxSubjectIndex = subjectIndex;
                            bestFrame = frame;
                        }
                    }

                }
                break;
            }
            default:
                throw new IllegalArgumentException("Unsupported alignment mode " + mode);
        }

        //frame starting with index 1 for output consistent
        FrameScore bestFrameScore = new FrameScore((bestFrame + 1), maxScore, maxQueryIndex, maxSubjectIndex);
        return bestFrameScore;
    }

    private static void printMatrix(FramebotSeqs seqs, FramebotScoringMatrix scoringMatrices){
        for (int frame = 0; frame < FramebotScoringMatrix.NUM_FRAMES; frame++) {
            System.err.println("FRAME\t" + (frame+1) + "\t" + seqs.queryNuclSeq.getSeqName() +"\t" + seqs.targetProtSeq.getSeqName()+ "\n"); 
            for (int queryIndex = 0; queryIndex < seqs.queryProtFrames[frame].length; queryIndex++) {
                for (int subjectIndex = 0; subjectIndex < seqs.targetProtSeqString.length; subjectIndex++) {
                    System.err.print(scoringMatrices.d[frame][queryIndex][subjectIndex] + " ");
                }
                System.err.println();
            }
        }
    }
    
    
    private static FramebotResult traceback(FramebotSeqs seqs, FramebotScoringMatrix scoringMatrices, AlignmentMode mode) {
        StringBuilder targetAlignment = new StringBuilder();   // the alignment of the subject protein
        StringBuilder queryAlignment = new StringBuilder();    // the alignment of the query protein
        List<Integer> frames = new ArrayList();
        StringBuilder correctedNucl = new StringBuilder();  // the corrected query nucl
        List<Integer> numberOfOrigBasesList = new ArrayList(); // the number of nucleotide bases from the original nucl seq for each corresponding AA
        FrameScore frameScore = getBestFrameScore(seqs, scoringMatrices, mode);
        
        byte[] origQual = null;
        List<Byte> qualitySeq = new ArrayList();
        if (seqs.queryNuclSeq instanceof QSequence) {
            origQual = ((QSequence) seqs.queryNuclSeq).getQuality();
        }

        int numOfFrameshifts = 0;

        int baseFrame = frameScore.getBestFrame() - 1;
        int i = frameScore.getQueryEnd();
        int j = frameScore.getTargetEnd();

        while (i >= 0 && j >= 0) {
            if (j==0 && mode != AlignmentMode.global) {
                break;
            }
            if (mode == AlignmentMode.global && i == 0 & j == 0) {
                break;
            }
            // take care of the global alignment 
            if ( i == 0){
                if (mode == AlignmentMode.global){
                    queryAlignment.append(GAP);
                    targetAlignment.append(seqs.targetProtSeqString[j]);
                    numberOfOrigBasesList.add(0);
                    frames.add(baseFrame + 1);
                    j--;
                    continue;
                }else {
                    break;
                }                
            }
                       
           
            frames.add(baseFrame + 1); // print frame starting with index 1
            int nextFrame = scoringMatrices.match_frame[baseFrame][i][j];
            int match_path = scoringMatrices.path[baseFrame][i][j];
            // for global alignment, the path should be DOWN and assuming not frameshift when target reaches the beginning but we did not initialize that way
            if ( j == 0 && mode == AlignmentMode.global){
                match_path = DOWN;
                nextFrame = baseFrame;
            }
            if (nextFrame != baseFrame) {
                numOfFrameshifts++;
            }

            // we need to find the correct frame offset
            int prevCellOffSet = (baseFrame < nextFrame) ? 2 : 1;
            if (mode == AlignmentMode.local && match_path == NEW_START) {
                break;
            }

            int tempBeginCharIndex = (i - 1) * FramebotScoringMatrix.NUM_FRAMES + baseFrame;
            int tempEndCharIndex = (i - 1) * FramebotScoringMatrix.NUM_FRAMES + baseFrame + 2;

            if (match_path == MATCH) {
                queryAlignment.append(seqs.queryProtFrames[baseFrame][i]);
                targetAlignment.append(seqs.targetProtSeqString[j]);

                //get the aligned nucleotides
                if (nextFrame == baseFrame) {
                    for (int f = FramebotScoringMatrix.NUM_FRAMES - 1; f >= 0; f--) {
                        correctedNucl.append(seqs.queryNuclSeqString[(i - 1) * FramebotScoringMatrix.NUM_FRAMES + baseFrame + f]);
                        if (origQual != null) {
                            qualitySeq.add(origQual[(i - 1) * FramebotScoringMatrix.NUM_FRAMES + baseFrame + f]);
                        }
                    }
                    numberOfOrigBasesList.add(3);
                    i--;
                } else {
                    i = i - prevCellOffSet;
                }

                j--;
            } else if (match_path == DOWN) {
                queryAlignment.append(seqs.queryProtFrames[baseFrame][i]);
                targetAlignment.append(GAP);
                //get the aligned nucleotides, it's an insertion
                if (nextFrame == baseFrame) {
                    for (int f = FramebotScoringMatrix.NUM_FRAMES - 1; f >= 0; f--) {
                        correctedNucl.append(seqs.queryNuclSeqString[(i - 1) * FramebotScoringMatrix.NUM_FRAMES + baseFrame + f]);
                        if (origQual != null) {
                            qualitySeq.add(origQual[(i - 1) * FramebotScoringMatrix.NUM_FRAMES + baseFrame + f]);
                        }
                    }
                    numberOfOrigBasesList.add(3);
                    i--;
                } else {
                    i = i - prevCellOffSet;
                }
            } else if (match_path == ACROSS) {
                targetAlignment.append(seqs.targetProtSeqString[j]);

                j--;
                //get the aligned nucleotides, it's deletion in this case
                if (nextFrame != baseFrame) {
                    i = i - prevCellOffSet + 1;
                    // if there is 2 nucleotides,insert an unknown aa
                    tempBeginCharIndex = i * FramebotScoringMatrix.NUM_FRAMES + nextFrame;
                    if ((tempEndCharIndex - tempBeginCharIndex + 1) == 2) {
                        queryAlignment.append(UNKNOWN_AA);
                    } else {
                        queryAlignment.append(GAP);
                    }

                } else {
                    queryAlignment.append(GAP);
                    numberOfOrigBasesList.add(0);
                }
            }

            // take care when frame changed
            if (nextFrame != baseFrame) {
                int align_nucl_tempBeginCharIndex = i * FramebotScoringMatrix.NUM_FRAMES + nextFrame;
                numberOfOrigBasesList.add((tempEndCharIndex - align_nucl_tempBeginCharIndex + 1));

                // we need to get the correct dna bases if it's not deletion
                if (match_path != ACROSS) {

                    for (int charIndex = tempEndCharIndex; charIndex >= tempBeginCharIndex; charIndex--) {
                        correctedNucl.append(seqs.queryNuclSeqString[charIndex]);
                        if (origQual != null) {
                            qualitySeq.add(origQual[charIndex]);
                        }
                    }
                } else {
                    // if it's deletion and has 2 nucl, insert NNN
                    if ((tempEndCharIndex - tempBeginCharIndex + 1) == 2) {
                        for (int f = FramebotScoringMatrix.NUM_FRAMES - 1; f >= 0; f--) {
                            correctedNucl.append(UNKNOWN_NUCL);
                            if (origQual != null) {
                                qualitySeq.add(UNKNOWN_NUCL_QSCORE);
                            }
                        }

                    }
                    // else leave it as blank
                }

            }

            baseFrame = nextFrame;

        }       

        //We have to be careful to reverse everything, since tracing back starts at the end of the alignment and builds toward the beginning
        Sequence alignedQuery = new Sequence(seqs.queryNuclSeq.getSeqName(), "", StringUtils.reverse(queryAlignment.toString()));
        Sequence alignedTarget = new Sequence(seqs.targetProtSeq.getSeqName(), "", StringUtils.reverse(targetAlignment.toString()));
        Sequence correctedNuclSeq = null;
        
        if(origQual == null) {
            correctedNuclSeq = new Sequence(seqs.queryNuclSeq.getSeqName(), "", StringUtils.reverse(correctedNucl.toString()));
        } else {
            Collections.reverse(qualitySeq);
            byte[] newQual = new byte[qualitySeq.size()];
            for(int index = 0;index < qualitySeq.size();index++) {
                newQual[index] = qualitySeq.get(index);
            }

            correctedNuclSeq = new QSequence(seqs.queryNuclSeq.getSeqName(), "", StringUtils.reverse(correctedNucl.toString()), newQual);
        }

        Collections.reverse(frames);
        Collections.reverse(numberOfOrigBasesList);
        frameScore.setQueryStart(i);
        frameScore.setTargetStart(j);
        return new FramebotResult(alignedQuery, alignedTarget, correctedNuclSeq, numOfFrameshifts, frames, numberOfOrigBasesList, frameScore);

    }
}
