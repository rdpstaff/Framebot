/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

import java.io.*;
import java.util.*;

/**
 *
 * @author wangqion
 */
public class GetFrameBotMatchMatrix {
    
    public static class Match{
        String subjectName;
        Match(String m){
            subjectName = m;
        }
    }
    
    public static class MatchComparator implements Comparator {
        public int compare(Object a, Object b){
            return ((Match)a).subjectName.compareToIgnoreCase( ((Match)b).subjectName);
        }
    }
    
    /**
     * take one or a directory of FrameBot alignment output(s), optional idMapping and sampleMapping file
     * return the number of sequences to the nearest match. The format is similar to a data matrix used for R analysis
     * samples in rows, match in columns
     * @param framebotResult
     * @param idMapping
     * @param sampleMapping
     * @throws IOException
     */
    public static void getFramebotMatrix(File framebotResult, String idMapping, String sampleMapping, String outFile, double identity) throws IOException{
        HashMap<String, HashSet<String>> sampleMap = SampleMappingReader.getSampleMapping(sampleMapping);
        HashMap<String, HashMap> countMap = IdMappingReader.getIDCount(idMapping, sampleMap);
        BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(outFile)));
        
        HashMap<String, HashMap<Match, Float>> matrixMap = new HashMap<String, HashMap<Match, Float>>();
        if ( framebotResult.isDirectory()){
            for ( File child : framebotResult.listFiles()){
                HashMap<String, Match> matchResultMap = getMatchMap(child, identity);
                if ( sampleMap == null){
                    matrixMap.put(GetFrameBotStatMain.DEFAULT_SAMPLE, getOneFramebotMatrix(matchResultMap, countMap.get(GetFrameBotStatMain.DEFAULT_SAMPLE)));
                }else {
                    for ( String sample: sampleMap.keySet()){
                        matrixMap.put(sample, getOneFramebotMatrix(matchResultMap, countMap.get(sample)));
                    }
                }
            }
        }else {
            HashMap<String, Match> matchResultMap = getMatchMap(framebotResult, identity);
            if ( sampleMap == null){
                matrixMap.put(GetFrameBotStatMain.DEFAULT_SAMPLE,getOneFramebotMatrix(matchResultMap, countMap.get(GetFrameBotStatMain.DEFAULT_SAMPLE)));
            }else {               
                for ( String sample: sampleMap.keySet()){
                    matrixMap.put(sample, getOneFramebotMatrix(matchResultMap, countMap.get(sample)));
                } 
            }
        }

        TreeSet<Match> commonSubjectSet = new TreeSet<Match>( new MatchComparator());
        for ( String sample: matrixMap.keySet()){
            commonSubjectSet.addAll(matrixMap.get(sample).keySet());
        }
       
        outWriter.write("matchName");
        for ( Match subject: commonSubjectSet){
            outWriter.write("\t" + subject.subjectName);
        }
        outWriter.write("\n");
        
       
        for ( String sample: matrixMap.keySet()){
             outWriter.write(sample);
            for ( Match subject: commonSubjectSet){
                if ( matrixMap.get(sample).get(subject) != null){
                    outWriter.write("\t" + matrixMap.get(sample).get(subject).floatValue());
                }else {
                    outWriter.write("\t" + "0");
                }
            }
            outWriter.write("\n");
        }

        outWriter.close();

    }
    
    private static HashMap<Match, Float> getOneFramebotMatrix(HashMap<String, Match> matchResultMap, HashMap<String, Integer> countMap ) throws IOException{
        HashMap<Match, Float> subjectMap = new HashMap<Match, Float>();
       	Set<String> idSet = new HashSet<String>();
        idSet.addAll( matchResultMap.keySet()); 
        if ( countMap != null){
            idSet.retainAll(countMap.keySet());
        }

        for ( String id: idSet) { 
            int mappingCount = 1;
            if ( countMap != null) {
                mappingCount = countMap.get(id);
            }
            // subject
            Float existingSubject = subjectMap.get(matchResultMap.get(id));
            if ( existingSubject == null){
                subjectMap.put(matchResultMap.get(id), (float)mappingCount );
            }else {
                subjectMap.put(matchResultMap.get(id), (float)(mappingCount + existingSubject.intValue()) );
            }
        }  
        return subjectMap;
    }
   
    private static HashMap<String, Match> getMatchMap(File framebotResult, double identity ) throws IOException{
        HashMap<String, Match> matchMap = new HashMap<String, Match>();
        HashMap<String, Match> subjectMap = new HashMap<String, Match>();
        FrameBotStatIterator iterator = new FrameBotStatIterator(framebotResult, true);
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            if ( stat.getIdentity() < identity) continue;
            Match match = subjectMap.get(stat.subjectID);
            if ( match == null ){
                match = new Match(stat.subjectID);
                subjectMap.put(stat.subjectID, match);
            }
            matchMap.put(stat.queryID, match);
        }
        iterator.close();
        
        return matchMap;
    }

    /**
     * takes one or a directory of FrameBot alignment output(s), sampleMapping file, and an optional idMapping file,
     * returns the number of sequences to the nearest match. 
     * The format is similar to a data matrix used for R analysis: samples in rows, match in columns.
     * If the sampleMapping file contains a subset of sequences, it outputs the results for the subset.
     * The subset sampleMapping can be generated by RdmSelectSampleMapping
     * @param framebotResult
     * @param idMapping
     * @param sampleMapping
     * @throws IOException
     */
     public static void getSubsetMatrix(File framebotResult, String idMapping, String subsetSampleMapping, String outFile, double identity) throws IOException{
        HashMap<String, HashSet<String>> sampleMap = SampleMappingReader.getSampleMapping(subsetSampleMapping);
       
        BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(outFile)));
        HashMap<String, Match> matchResultMap = getMatchMap(framebotResult, identity);
        HashMap<String, HashMap<Match, Float>> matrixMap = new HashMap<String, HashMap<Match, Float>>();
       
        if ( idMapping == null){
            for ( String sample: sampleMap.keySet()){
                HashSet<String> sample_seqset = sampleMap.get(sample);
                HashMap<Match, Float> tempMap = matrixMap.get(sample);
                if ( tempMap == null){
                    tempMap = new HashMap<Match, Float>();
                    matrixMap.put(sample, tempMap);
                }
                for ( String id: sample_seqset){
                     Match match = matchResultMap.get(id);
                     if ( match == null) { // sequence may failed FrameBot
                         continue; 
                     }
                     Float count = tempMap.get(match);
                     if ( count == null){
                         tempMap.put(match, 1.0f);
                     }else {
                         tempMap.put(match, count.floatValue() + 1.0f);
                     }
                }
            }
        }else {
            BufferedReader reader = new BufferedReader(new FileReader(new File(idMapping)));
            String line = null;

            while ( (line=reader.readLine()) != null){
                if ( line.trim().equals("")) continue;
                String [] values = line.split("\\s+");
                String[] ids = values[1].split(",");
                Match match = matchResultMap.get(ids[0]);
                if ( match == null){  // sequence may failed FrameBot
                    continue;
                }

                HashSet<String> seqset = new HashSet<String>();
                for ( String id: ids){
                    seqset.add(id);
                }

                for ( String sample: sampleMap.keySet()){
                    HashSet sample_seqset = sampleMap.get(sample);
                    Set<String> commonSet = new HashSet<String>(seqset);
                    commonSet.retainAll(sample_seqset);
                    if ( commonSet.size() == 0) continue;

                    HashMap<Match, Float> tempMap = matrixMap.get(sample);
                    if ( tempMap == null){
                        tempMap = new HashMap<Match, Float>();
                        matrixMap.put(sample, tempMap);
                    }
                    for ( String id: commonSet){
                        Float count = tempMap.get(match);
                        if ( count == null){
                            tempMap.put(match, 1.0f);
                        }else {
                            tempMap.put(match, count.floatValue() + 1.0f);
                        }
                    }
                }
            }
        }
        
        TreeSet<Match> commonSubjectSet = new TreeSet<Match>( new MatchComparator());
        for ( String sample: matrixMap.keySet()){
            commonSubjectSet.addAll(matrixMap.get(sample).keySet());
        }
        
        outWriter.write("matchName");
        for ( Match subject: commonSubjectSet){
            outWriter.write("\t" + subject.subjectName);
        }
        outWriter.write("\n");
       
        for ( String sample: matrixMap.keySet()){
             outWriter.write(sample);
            for ( Match subject: commonSubjectSet){
                if ( matrixMap.get(sample).get(subject) != null){
                    outWriter.write("\t" + matrixMap.get(sample).get(subject).floatValue());
                }else {
                    outWriter.write("\t" + "0");
                }
            }
            outWriter.write("\n");
        }

        outWriter.close();
     }
}
