/*
 * Copyright (C) 2014 wangqion
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

package edu.msu.cme.rdp.framebot.stat;

import edu.msu.cme.rdp.readseq.stat.StdevCal;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author wangqion
 */
public class CompareFrameBotResult {
    
    /**
     * 
     * @param args
     * @throws IOException 
     */
    public static void main(String[] args) throws IOException{
        String usage = "framebot_result1.txt framebot_result2.txt\n" 
                + "This program compares two framebot results of the same query sequences";
        
        HashMap<String, FrameBotStat> result2_statMap = new HashMap<String, FrameBotStat>();
        File file2 = new File(args[1]);
        FrameBotStatIterator iterator = new FrameBotStatIterator(file2, true);
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            result2_statMap.put(stat.getQueryID(), stat);
        }
        iterator.close();
        
        int total = 0;
        int found = 0;
        int same_target = 0;
        int same_fshift = 0;
        int same_identity = 0;
        int same_score = 0;
        ArrayList<Double> higher_identity_result1 = new ArrayList<Double>();
        ArrayList<Double> higher_identity_result2 = new ArrayList<Double>();
        
        ArrayList<Double> higher_fshifts_result1 = new ArrayList<Double>();
        ArrayList<Double> higher_fshifts_result2 = new ArrayList<Double>();
        
        ArrayList<Double> higher_score_result1 = new ArrayList<Double>();
        ArrayList<Double> higher_score_result2 = new ArrayList<Double>();
        
        File file1 = new File(args[0]);
        iterator = new FrameBotStatIterator(file1, true);
        while (iterator.hasNext()){
            total ++;
            FrameBotStat stat_1 = iterator.next();
            FrameBotStat stat_2 = result2_statMap.get(stat_1.getQueryID());
            
            if ( stat_2 == null){
                System.err.println("query " + stat_1.getQueryID() + " not found in " + args[1]);
            }else {
                found ++;
                if ( stat_1.getSubjectID().equals(stat_2.getSubjectID())){
                    same_target ++;
                }
                if ( stat_1.getFrameshifts() == stat_2.getFrameshifts()){
                    same_fshift ++;
                }else {
                    if ( stat_1.getFrameshifts() > stat_2.getFrameshifts()){
                        higher_fshifts_result1.add( new Double(stat_1.getFrameshifts() - stat_2.getFrameshifts()) );
                    }else {
                        higher_fshifts_result2.add( new Double(stat_2.getFrameshifts() - stat_1.getFrameshifts()) );
                    }
                }   
                
                if ( stat_1.getScore()== stat_2.getScore()){
                    same_score ++;
                }else {
                    //System.err.println(stat_1.getQueryID() + "\t" + stat_1.getScore() + "\t" +  stat_2.getScore() + "\tframeshift\t" + stat_1.getFrameshifts());
                    if ( stat_1.getScore() > stat_2.getScore()){
                        higher_score_result1.add( new Double(stat_1.getScore() - stat_2.getScore()) );
                    }else {
                        higher_score_result2.add( new Double(stat_2.getScore() - stat_1.getScore()) );
                    }
                }   
                
                if ( stat_1.getIdentity() == stat_2.getIdentity()){
                    same_identity ++;
                }else {
                    if ( stat_1.getIdentity() > stat_2.getIdentity()){                        
                        higher_identity_result1.add(stat_1.getIdentity() - stat_2.getIdentity());
                    }else {
                        higher_identity_result2.add(stat_2.getIdentity() - stat_1.getIdentity());
                    }
                }
            }
        }
        iterator.close();
        
        // 
        System.out.println( "File1" + "\t" + "File2" + "\t" + "same_target" + "\t" + "same_score"+ "\t" + "same_identity" + "\t" + "same_frameshift" + "\t" + "samescore%" );
        System.out.println(file1.getName() + "\t" + file2.getName() + "\t" + total + "\t" + found + "\t" + same_target + "\t" + same_score + "\t" + same_identity + "\t" + same_fshift + "\t" + same_score/((float)found) );
        
        System.out.println("\nListName\tTotalCount\tMean\tStdev");
        StdevCal.Std std = StdevCal.calStd(higher_score_result1);
        System.out.println("higher_score_result1\t" + std.getTotalCount() + "\t" + String.format("%.3f", std.getMean()) + "\t" + String.format("%.3f", std.getStdev()) );
        
        std = StdevCal.calStd(higher_score_result2);
        System.out.println("higher_score_result2\t" + std.getTotalCount() + "\t" + String.format("%.3f", std.getMean()) + "\t" + String.format("%.3f", std.getStdev()) );
     
        
        std = StdevCal.calStd(higher_identity_result1);
        System.out.println("higher_identity_result1\t" + std.getTotalCount() + "\t" + String.format("%.3f", std.getMean()) + "\t" + String.format("%.3f", std.getStdev()) );
        
        std = StdevCal.calStd(higher_identity_result2);
        System.out.println("higher_identity_result2\t" + std.getTotalCount() + "\t" + String.format("%.3f", std.getMean()) + "\t" + String.format("%.3f", std.getStdev()) );
        
        std = StdevCal.calStd(higher_fshifts_result1);
        System.out.println("higher_fshifts_result1\t" + std.getTotalCount() + "\t" + String.format("%.3f", std.getMean()) + "\t" + String.format("%.3f", std.getStdev()) );
        
        std = StdevCal.calStd(higher_fshifts_result2);
        System.out.println("higher_fshifts_result2\t" + std.getTotalCount() + "\t" + String.format("%.3f", std.getMean()) + "\t" + String.format("%.3f", std.getStdev()) );
        
    }
    
}
