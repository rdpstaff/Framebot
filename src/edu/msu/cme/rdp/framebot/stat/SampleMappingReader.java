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
public class SampleMappingReader {
    public static HashMap getSampleMapping( String sampleMapping ) throws IOException {
        if ( sampleMapping == null) return null;
        HashMap<String, HashSet<String>> sampleMap = new HashMap<String, HashSet<String>>();
        
        BufferedReader sampleReader = new BufferedReader(new FileReader(new File(sampleMapping)));
        String line = null;
        while ( (line=sampleReader.readLine()) != null){
            if ( line.trim().equals("")) continue;
            String [] values = line.split("\\s+");
            if ( sampleMap.get(values[1]) != null){
                sampleMap.get(values[1]).add(values[0]);
            }else {
                HashSet<String> tmp = new HashSet<String>();
                tmp.add(values[0]);
                sampleMap.put(values[1], tmp);
            }
        }
        sampleReader.close();
        return sampleMap;
    }
    
}
