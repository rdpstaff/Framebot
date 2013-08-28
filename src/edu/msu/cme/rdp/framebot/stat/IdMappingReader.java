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
import java.util.Set;

/**
 *
 * @author wangqion
 */
public class IdMappingReader {
    /**
     * get the abundance for each rep sequence for each sample from the idMapping file and sampleMapping file
     * @param idmapping
     * @return
     * @throws IOException
     */
     public static HashMap getIDCount(String idmapping, HashMap<String, HashSet<String>> sampleMap ) throws IOException {
         HashMap<String, HashMap<String, Integer>> countMap = new HashMap<String, HashMap<String, Integer>>();
        if ( idmapping == null) {
            if ( sampleMap == null) {
                return countMap;
            } else { // if idmapping is null, we should assume seq in sampleMap have one count
            
                for ( String sample: sampleMap.keySet()){
                    Set<String> sample_seqset = sampleMap.get(sample);
                    HashMap<String, Integer> tmp = countMap.get(sample);
                    if ( tmp == null){
                        tmp = new HashMap<String, Integer>();
                        countMap.put(sample, tmp);
                    }
                    for ( String id: sample_seqset){
                        tmp.put( id, 1);
                    }
                    
                }
            }
            return countMap;
        }
        
        BufferedReader reader = new BufferedReader(new FileReader(new File(idmapping)));
        String line = null;

        while ( (line=reader.readLine()) != null){
            if ( line.trim().equals("")) continue;
            String [] values = line.split("\\s+");
            String[] ids = values[1].split(",");
            
            if ( sampleMap == null){
                HashMap<String, Integer> tmp = countMap.get(GetFrameBotStatMain.DEFAULT_SAMPLE);
                if ( tmp == null){
                    tmp = new HashMap<String, Integer>();
                    countMap.put(GetFrameBotStatMain.DEFAULT_SAMPLE, tmp);
                }
                tmp.put( ids[0], ids.length);

            }else {
                HashSet<String> seqset = new HashSet<String>();
                for ( String id: ids){
                    seqset.add(id);
                }
                int prev_size = seqset.size();
                for ( String sample: sampleMap.keySet()){
                    Set sample_seqset = sampleMap.get(sample);
                    if ( sample_seqset != null)
                        seqset.removeAll(sample_seqset);

                    HashMap<String, Integer> tmp = countMap.get(sample);
                    if ( tmp == null){
                        tmp = new HashMap<String, Integer>();
                        countMap.put(sample, tmp);
                    }
                    if ( (prev_size - seqset.size() ) > 0 ){
                        tmp.put( ids[0], prev_size - seqset.size());
                    }
                    prev_size = seqset.size();
                }
            }
        }
        reader.close();
        return countMap;
    }
    
     
     /**
     * get the mapping sequence IDs of the rep sequences for each sample from the idMapping and sampleMapping file
     * @param idmapping
     * @return
     * @throws IOException
     */
     public static HashMap getIDMapping(String idmapping, HashMap<String, HashSet<String>> sampleMap ) throws IOException {
         HashMap<String, HashMap<String, Set>> countMap = new HashMap<String, HashMap<String, Set>>();
        if ( idmapping == null) return countMap;
                
        BufferedReader reader = new BufferedReader(new FileReader(new File(idmapping)));
        String line = null;
        
        while ( (line=reader.readLine()) != null){
            if ( line.trim().equals("")) continue;
            String [] values = line.split("\\s+");
            String[] ids = values[1].split(",");
            HashSet<String> seqset = new HashSet<String>();
            for ( String id: ids){
                seqset.add(id);
            }
            if ( sampleMap == null){
                HashMap<String, Set> tmp = countMap.get(GetFrameBotStatMain.DEFAULT_SAMPLE);
                if ( tmp == null){
                    tmp = new HashMap<String, Set>();
                    countMap.put(GetFrameBotStatMain.DEFAULT_SAMPLE, tmp);
                }               
                tmp.put( ids[0], seqset);
                countMap.put(GetFrameBotStatMain.DEFAULT_SAMPLE, tmp);
            }else {

                for ( String sample: sampleMap.keySet()){
                    HashSet sample_seqset = sampleMap.get(sample);
                    Set commonSet = new HashSet(seqset);
                    commonSet.retainAll(sample_seqset);
                    
                    HashMap<String, Set> tmp = countMap.get(sample);
                    if ( tmp == null){
                        tmp = new HashMap<String,Set>();
                        countMap.put(sample, tmp);
                    }
                    tmp.put(ids[0], commonSet);
                }
            }
        }
        reader.close();
        return countMap;
    }
}
