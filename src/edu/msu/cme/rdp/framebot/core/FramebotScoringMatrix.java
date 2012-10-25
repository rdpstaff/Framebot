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
package edu.msu.cme.rdp.framebot.core;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;

/**
 *
 * @author fishjord
 */
public class FramebotScoringMatrix {

    public static final int NUM_FRAMES = 3;
    int[][][] d;
    /** the alignment score matrix */
    int[][][] p;
    /** the insertion score matrix */
    int[][][] q;
    /** the deletion score matrix */
    int[][][] match_frame;
    /** the frame of current match */
    int[][][] path;
    /** the path taken to get the score */

    /** the path taken to get the score */
    public FramebotScoringMatrix(int queryLength, int subjectLength, AlignmentMode mode, ScoringMatrix simMatrix) {
        d = new int[NUM_FRAMES][queryLength][subjectLength];
        p = new int[NUM_FRAMES][queryLength][subjectLength];
        q = new int[NUM_FRAMES][queryLength][subjectLength];
        match_frame = new int[NUM_FRAMES][queryLength][subjectLength];
        path = new int[NUM_FRAMES][queryLength][subjectLength];
        int initialVal =  simMatrix.getGapOpen();
        int extendVal = simMatrix.getGapExtend();
        switch (mode) {
            case local: {
                for (int i = 0; i < NUM_FRAMES; i++) {
                    for (int j = 0; j < subjectLength; j++) {
                        d[i][0][j] = 0;
                        p[i][0][j] = 0;
                        q[i][0][j] = 0;
                    }
                    for (int j = 0; j < queryLength; j++) {
                        d[i][j][0] = 0;
                        p[i][j][0] = 0;
                        q[i][j][0] = 0;
                    }
                }
                break;
            }

            case global:{
                for (int i = 0; i < NUM_FRAMES; i++) {
                    for (int j = 0; j < subjectLength; j++) {
                        d[i][0][j] = initialVal + extendVal*j;
                        p[i][0][j] = initialVal + extendVal*j;
                        q[i][0][j] = initialVal + extendVal*j;
                    }
                    for (int j = 0; j < queryLength; j++) {
                        d[i][j][0] = initialVal + extendVal*j;
                        p[i][j][0] = initialVal + extendVal*j;
                        q[i][j][0] = initialVal + extendVal*j;
                    }
                    d[i][0][0] = 0;
                    p[i][0][0] = 0;
                    q[i][0][0] = 0;  
                }     
                 break;
            }
            case glocal:{
                // local on the subject, global on the query
                 for (int i = 0; i < NUM_FRAMES; i++) {
                    for (int j = 0; j < subjectLength; j++) {
                        d[i][0][j] = initialVal;
                        p[i][0][j] = initialVal;
                        q[i][0][j] = initialVal;
                    }
                    for (int j = 0; j < queryLength; j++) {
                        d[i][j][0] = initialVal + extendVal*j;
                        p[i][j][0] = initialVal + extendVal*j;
                        q[i][j][0] = initialVal + extendVal*j;
                    }
                    
                }
                break;
            }
            default:
                throw new IllegalArgumentException("Unsupported alignment mode " + mode);
        }
       
    }
}
