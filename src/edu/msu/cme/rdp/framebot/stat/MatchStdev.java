/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.framebot.stat;

import edu.msu.cme.rdp.readseq.stat.StdevCal;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author wangqion
 */
public class MatchStdev {

    public static HashMap<String,String> readMap(String file) throws IOException {
        HashMap<String, String> map = new HashMap<String, String> ();
        BufferedReader reader = new BufferedReader(new FileReader(new File(file)));
        String line = null; 
        while ( (line = reader.readLine()) != null){
            if ( line.trim().equals("")) break;
            String [] vals = line.split("\\t");
            map.put(vals[0].trim(), vals[1].trim());
        }
        reader.close();
        return map;
    }
    
    
    /**
     * Input file: FrameBot matrix output from GetFrameBotStatMain class
     * For each sample group, and match group, calculate the mean and stdev of the percent of the match 
     * @param args
     * @throws IOException 
     */
    public static void main(String[] args) throws IOException {
        String usage = "MatchStdev matrix.txt outfile sample_group match_group\n" +
                " sample_group, tab-delimited file contains the sample name and sample group designation" +
                " match_group, tab-delimited file contains the match name and the match group designation";
        
        
        HashMap<String, String> sampleGroupMap = readMap(args[2]);
        HashMap<String, String> matchGroupMap = readMap(args[3]);
        // key=sampleGroupName, value=map [value=matchGroup, value=map[sampleName, pct]]
        HashMap<String, HashMap<String, HashMap<String, Double>>> sampleCountMap = new HashMap<String, HashMap<String, HashMap<String, Double>>>();  
        // find the unique match group set that also present in the matrix file
        HashSet<String> uniqueMatchGroupSet = new HashSet();
        
        BufferedReader reader = new BufferedReader(new FileReader(new File(args[0])));
        // first line is header of match names, first column is the sample name
        String[] header = reader.readLine().split("\\t");
        String line = null;
        // sanity check, make sure all match names listed in the match group
        for ( int i = 1; i < header.length; i++){
            if ( !matchGroupMap.keySet().contains(header[i].trim())){
                throw new IllegalArgumentException("match name \"" + header[i] + "\" in the matrix file not found in the match_group file");
            }
            uniqueMatchGroupSet.add( matchGroupMap.get(header[i].trim()));
        }
        for ( String sampleGroup: sampleGroupMap.values()){
            if ( !sampleCountMap.containsKey(sampleGroup)) {
                sampleCountMap.put(sampleGroup, new HashMap<String, HashMap<String, Double>>());
                for ( String matchGroup: uniqueMatchGroupSet) {
                    sampleCountMap.get(sampleGroup).put(matchGroup, new HashMap<String, Double>());
                } 
            }
        }
        
        while ( (line=reader.readLine()) != null){
            String[] values = line.split("\\t");
            String sampleName = values[0];
            if ( values.length != header.length) {
                throw new IllegalArgumentException("the number of values in line " + line + " should have " + header.length + " columns");
            }
            if ( !sampleGroupMap.containsKey(values[0])) {
                throw new IllegalArgumentException("sample name \"" + values[0] + "\" in the matrix file not found in the sample_group file");
            }
            
            double total = 0.0;
            for ( int i = 1; i < values.length; i++){
                total += Double.parseDouble(values[i]);
            }
            HashMap<String, HashMap<String, Double>> matchCountMap = sampleCountMap.get(sampleGroupMap.get(sampleName));
            for ( int i = 1; i < values.length; i++){
                Double existVal = matchCountMap.get(matchGroupMap.get(header[i])).get(sampleName);
                if ( existVal != null){
                    matchCountMap.get(matchGroupMap.get(header[i])).put(sampleName, existVal + Double.parseDouble(values[i])/total);
                }else {
                    matchCountMap.get(matchGroupMap.get(header[i])).put(sampleName, Double.parseDouble(values[i])/total);
                }
            }
        }
        
        reader.close();
        
        // calculate mean and stdev
        HashMap<String, HashMap<String, StdevCal.Std>> stdevMap = new HashMap<String, HashMap<String, StdevCal.Std>>();  // key=sampleName, value=map
        for ( String sampleGroup: sampleCountMap.keySet()){
            HashMap<String, HashMap<String, Double>> matchCountMap = sampleCountMap.get(sampleGroup);
            stdevMap.put(sampleGroup, new HashMap<String, StdevCal.Std>());
            for ( String matchGroup: matchCountMap.keySet()){
                // need to collect all the pct to a list 
                ArrayList<Double> pctList = new ArrayList<Double>();
                HashMap<String, Double> h = matchCountMap.get(matchGroup);
                for ( String sampleName: h.keySet()){
                    pctList.add(h.get(sampleName));
                }
                StdevCal.Std result = StdevCal.calStd(pctList);
                stdevMap.get(sampleGroup).put(matchGroup, result);
            }
        }

        PrintWriter outWriter = new PrintWriter(new File(args[1]));
        // print the mean
        outWriter.print("mean");
        for ( String sampleGroup: stdevMap.keySet()){
             outWriter.print( "\t" + sampleGroup  );
        }
        
        outWriter.println("\t"); // this extra tab preventing the last char  not shown in excel, an excel bug
        for ( String matchGroup: uniqueMatchGroupSet){
            outWriter.print( matchGroup);
            for ( String sampleGroup: stdevMap.keySet()){
                HashMap<String, StdevCal.Std> resultMap = stdevMap.get(sampleGroup);
                outWriter.print( "\t"+ resultMap.get(matchGroup).getMean());                
            }
            outWriter.println("\t");
        }
        
        outWriter.println("\t");
        // print the stdev
        outWriter.print("stdev");
        for ( String sampleGroup: stdevMap.keySet()){
             outWriter.print( "\t" + sampleGroup  );
        }
        outWriter.println("\t"); // this extra tab preventing the last char  not shown in excel, an excel bug
        for ( String matchGroup: uniqueMatchGroupSet){
            outWriter.print( matchGroup);
            for ( String sampleGroup: stdevMap.keySet()){
                HashMap<String, StdevCal.Std> resultMap = stdevMap.get(sampleGroup);
                outWriter.print( "\t"+ resultMap.get(matchGroup).getStdev());                
            }
            outWriter.println("\t");
        }
        outWriter.close();
    }

}
