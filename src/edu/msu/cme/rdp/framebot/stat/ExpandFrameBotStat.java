/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author wangqion
 */
public class ExpandFrameBotStat {
    
    public static HashSet<String> readIds(String idFile) throws IOException{
        HashSet<String> idSet = new HashSet<String>();
        BufferedReader reader = new BufferedReader(new FileReader(new File(idFile)));
        String line = null;
        while ((line=reader.readLine()) != null){
            String[] values = line.split("\\t");
            idSet.add(values[0].trim());
        }
        reader.close();
        return idSet;
    }
    
    
    public  static void main(String[] args) throws Exception {
        String usage = "selectedIdFile framebot_result [idmapping_file>] \n" +
                " The output contains the STAT line of FrameBot alignment result of the selected sequences, and to each sequence in the mappping file";
        if ( args.length != 2 && args.length != 3) {
            throw new RuntimeException(usage);
        }
        
        HashSet<String> idSet = readIds(args[0]);
        
        FrameBotStatIterator iterator = new FrameBotStatIterator(new File(args[1]));
        HashMap<String, HashSet<String>> idMap = null;
        if ( args.length == 3){
            HashMap<String, HashMap> countMap = IdMappingReader.getIDMapping(args[2], null);
             idMap = countMap.get(GetFrameBotStatMain.DEFAULT_SAMPLE);
        }
        
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            
            if ( idSet.contains(stat.queryID)){
                if ( idMap == null){
                    System.out.println("STATS\t" + stat.getSubjectID() + "\t" + stat.queryID + "\t" + 
                            stat.getNuclLen() + "\t" + stat.getAlignLen() +"\t" + stat.getIdentity() + "\t" + 
                            stat.getScore() + "\t" + stat.getFrameshifts() +"\t" + stat.isReversed());
                }else if ( idMap.get(stat.queryID) != null){               
                    for ( String id: idMap.get(stat.queryID)){
                        System.out.println("STATS\t" + stat.getSubjectID() + "\t" + id + "\t" + 
                                stat.getNuclLen() + "\t" + stat.getAlignLen() +"\t" + stat.getIdentity() + "\t" + 
                                stat.getScore() + "\t" + stat.getFrameshifts() +"\t" + stat.isReversed());
                    }
                }
            }
        }
        
        iterator.close();
    }
}
