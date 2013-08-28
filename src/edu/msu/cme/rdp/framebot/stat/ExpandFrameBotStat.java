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
        String usage = "selectedIdFile framebot_result ignoreAlignment [idmapping_file>] \n" +
                " if ignoreAlignment is tru, the output contains the STAT line of FrameBot alignment result of the selected sequences,\n" +
                " if ignoreAlignment is false, the output contains the full FrameBot alignment result of the selected sequences, expanded to each sequence in the mappping file";
        if ( args.length != 3 && args.length != 4) {
            throw new RuntimeException(usage);
        }
        
        
        HashSet<String> idSet = readIds(args[0]);
        boolean ignoreAlignment = Boolean.parseBoolean(args[2]);
        FrameBotStatIterator iterator = new FrameBotStatIterator(new File(args[1]), ignoreAlignment);
        HashMap<String, HashSet<String>> idMap = null;
        if ( args.length == 4){
            HashMap<String, HashMap> countMap = IdMappingReader.getIDMapping(args[3], null);
             idMap = countMap.get(GetFrameBotStatMain.DEFAULT_SAMPLE);
        }
        
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            
            if ( idSet.contains(stat.queryID)){
                if ( idMap == null){
                    if ( ignoreAlignment ){
                        System.out.println(stat.statToString());
                    }else {
                        System.out.println(stat.alignment.toString());
                    }
                }else if ( idMap.get(stat.queryID) != null){               
                    for ( String id: idMap.get(stat.queryID)){
                        if ( ignoreAlignment ){
                            System.out.println("STATS\t" + stat.getSubjectID() + "\t" + id + "\t" + 
                                stat.getNuclLen() + "\t" + stat.getAlignLen() +"\t" + stat.getIdentity() + "\t" + 
                                stat.getScore() + "\t" + stat.getFrameshifts() +"\t" + stat.isReversed());
                        }else {
                            System.out.println(">" + stat.getAlignment()[0]);
                            System.out.println("STATS\t" + stat.getSubjectID() + "\t" + id + "\t" + 
                                stat.getNuclLen() + "\t" + stat.getAlignLen() +"\t" + stat.getIdentity() + "\t" + 
                                stat.getScore() + "\t" + stat.getFrameshifts() +"\t" + stat.isReversed());
                            for ( int i = 2; i < stat.getAlignment().length; i++){
                                System.out.println(stat.getAlignment()[i]);
                            }
                            System.out.println();
                        }
                    }
                }
            }
        }
        
        iterator.close();
    }
}
