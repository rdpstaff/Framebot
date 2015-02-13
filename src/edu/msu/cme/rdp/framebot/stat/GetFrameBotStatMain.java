/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.stat;

import edu.msu.cme.rdp.readseq.SequenceFormat;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.*;
import java.util.*;
import org.apache.commons.cli.*;

/**
 *
 * @author wangqion
 */
public class GetFrameBotStatMain {
    
    
    private static Options options = new Options();
    private static final String ID_MAPPING_SHORT_OPT = "i";
    private static final String SAMPLE_MAPPING_SHORT_OPT = "s";
    private static final String STAT_SHORT_OPT = "t";
    private static final String DESCRIPTION_SHORT_OPT = "d";
    private static final String IDENTITY_SHORT_OPT = "e";
    
    private static final String ID_MAPPING_LONG_OPT = "id-mapping";
    private static final String SAMPLE_MAPPING_LONG_OPT = "sample-mapping";
    private static final String STAT_LONG_OPT = "stat-type";
    private static final String DESCRIPTION_LONG_OPT = "seq-desc";
    private static final String IDENTITY_LONG_OPT = "identity";
    
   
    private static final String STAT_DESC =  "stat | hist | summary | matrix | subset "
            + "\n stat ouptuts the # of seqs passed FrameBot, # of frameshifts for each sample"
            + "\n hist outputs a nearest match refseq, description and # of seqs close to the refseq at different identity% ranges"
            + "\n summary outputs a list of subject(refseq), description and  # of seqs close to the subject"
            + "\n matrix outputs the number of sequences to the nearest match. The format is similar to a data matrix used for R analysis"
            + "\n subset outputs the number of sequences to the nearest match for only sequence IDs in sample mapping file ";
            
      
    public static final String DEFAULT_SAMPLE = "default";
    private static HashMap<String, String> descMap = new HashMap<String, String>();
    
    static {
        options.addOption(new Option(ID_MAPPING_SHORT_OPT, ID_MAPPING_LONG_OPT, true, "Id mapping file. Output from Dereplicator (http://fungene.cme.msu.edu/FunGenePipeline/derep/form.spr)."));
        options.addOption(new Option(SAMPLE_MAPPING_SHORT_OPT, SAMPLE_MAPPING_LONG_OPT, true, "Sample mapping file. Output from Dereplicator (http://fungene.cme.msu.edu/FunGenePipeline/derep/form.spr)."));
        options.addOption(new Option(STAT_SHORT_OPT, STAT_LONG_OPT, true, STAT_DESC));
        options.addOption(new Option(DESCRIPTION_SHORT_OPT, DESCRIPTION_LONG_OPT, true, "the description of the reference seq, tab-delimited file or fasta"));
        options.addOption(new Option(IDENTITY_SHORT_OPT, IDENTITY_LONG_OPT, true, "the minimum protein identity, default is 0, range [0-100]"));
    }

    public static class Count{
        String subject;
        int count;
        Count(String s, int c){
            subject = s;
            count = c;
        }
    }
    
    
    public static class ResultComparator implements Comparator {
        public int compare(Object a, Object b){
            if ( ((Count)a).count < ((Count)b).count){
                return 1;
            }else {
                if ( ((Count)a).count > ((Count)b).count){
                    return -1;
                }else {
                    return ((Count)a).subject.compareTo(((Count)b).subject);
                }
            }
        }
    }

    public static HashMap<String, String> readDesc(File descFile) throws IOException{
        HashMap<String, String> descMap = new HashMap<String, String>();
        SequenceFormat format = SeqUtils.guessFileFormat(descFile);
        if ( format == SequenceFormat.UNKNOWN){
            BufferedReader reader = new BufferedReader(new FileReader(descFile));
            String line = null;
            while ((line=reader.readLine()) != null){
                String[] values = line.split("\\t");
                descMap.put(values[0], values[1]);
            }
            reader.close();
        }else {
            SequenceReader reader = new SequenceReader(descFile);       
            Sequence seq ;
            while ( (seq = reader.readNextSequence()) !=  null) {
                descMap.put( seq.getSeqName(), seq.getDesc());
            }
            reader.close();
        }
        return descMap;
    }

    /**
     * take a directory of framebot output, an optional idMapping and sampleMapping files
     * @param framebotResult file or directory
     * @param idMapping
     * @param sampleMapping
     * @throws IOException
     */
    public static void process(File framebotResult, String idMapping, String sampleMapping, String outFile, String type, double identity ) throws IOException{
        HashMap<String, HashSet<String>> sampleMap = SampleMappingReader.getSampleMapping(sampleMapping);
        HashMap<String, HashMap> countMap;
        countMap = IdMappingReader.getIDCount(idMapping, sampleMap);
        
        BufferedWriter outWriter = new BufferedWriter(new FileWriter(new File(outFile)));
        if( type.equalsIgnoreCase("stat")){
            outWriter.write("FileName" + "\t" + "SampleName" + "\t" + "# Seqs passed FrameBot" +"\t" + "# Seqs containing frameshift" +"\t" + "Total frameshifts" +"\n"  );
        }
        if ( framebotResult.isDirectory()){
            for ( File child : framebotResult.listFiles()){
                if ( sampleMap == null){
                    outWriter.write(processOneFramebot(child, countMap.get(DEFAULT_SAMPLE), DEFAULT_SAMPLE, type, identity)+"\n");
                }else {
                    for ( String sample: sampleMap.keySet()){
                        outWriter.write(processOneFramebot(child, countMap.get(sample), sample, type, identity) +"\n");
                    }
                }
            }
        }else {
            if ( sampleMap == null){
                outWriter.write(processOneFramebot(framebotResult, countMap.get(DEFAULT_SAMPLE), DEFAULT_SAMPLE, type, identity)+"\n");
            }else {
                for ( String sample: sampleMap.keySet()){
                    outWriter.write(processOneFramebot(framebotResult, countMap.get(sample), sample, type, identity) +"\n");
                }
            }
            
        }
        outWriter.close();
    }
  
    private static String processOneFramebot(File framebotResult, HashMap countMap, String sampleName, String type, double identity) throws IOException{
        if( type.equalsIgnoreCase("stat")){
            return getOneFramebotStat(framebotResult, countMap, sampleName, identity);
        }else if( type.equalsIgnoreCase("hist")){
            return getOneFramebotHist(framebotResult, countMap, sampleName, identity);
        }else if( type.equalsIgnoreCase("summary")){
            return getOneFramebotSummary(framebotResult, countMap, sampleName, identity);
        }
        return null;
    }
    
    /**
     * This reads in one frambot output and an optional abundance mapping, output the total passed, total frameshifts 
     * # of seqs, Frameshift
     * @param frambotResult file
     * @param countMap
     */
    private static String getOneFramebotStat(File framebotResult, HashMap<String, Integer> countMap, String sampleName, double identity) throws IOException{
        HashMap<Integer, Integer> frameshiftMap = new HashMap<Integer, Integer>();
        FrameBotStatIterator iterator = new FrameBotStatIterator(framebotResult, true);
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            if ( stat.getIdentity() < identity) continue;
            int mappingCount = 1;
            
            if ( countMap != null) {
                if ( countMap.get(stat.queryID)!= null ){
                    mappingCount = countMap.get(stat.queryID);
                }else {  // not in the mapping, that means we should discard this match
                    continue;
                }
            }     
            Integer existingFrameshift = frameshiftMap.get(stat.frameshifts);
            if ( existingFrameshift == null){
                frameshiftMap.put(stat.frameshifts, mappingCount );
            }else {
                frameshiftMap.put(stat.frameshifts, mappingCount + existingFrameshift.intValue() );
            }
        }
        iterator.close();

        int numPassedFramebot = 0;
        int numContainframeshift = 0;
        int totalNumOfFrameshift = 0;
         for ( Integer frameshift: frameshiftMap.keySet()){
            numPassedFramebot += frameshiftMap.get(frameshift);
            if ( frameshift.intValue() > 0) {
                numContainframeshift += frameshiftMap.get(frameshift);
                totalNumOfFrameshift += frameshift.intValue()*frameshiftMap.get(frameshift);
            }
        }
        return (framebotResult + "\t" + sampleName + "\t" + numPassedFramebot + "\t" + numContainframeshift +"\t" + totalNumOfFrameshift);

    }
    
    
    /**
     * This read in one frambot output calculate the seqcount for each identity and alignment length, good for graph
     * bin the mapping count for each reference at the %identity
     *  0-9, 11-19, ... 90-99, 100
     * @param frambotresult
     * @param idMapping
     */
    private static String getOneFramebotHist(File framebotResult, HashMap<String, Integer> countMap, String sampleName, double identity) throws IOException{
        HashMap<String, int[]> identityMap = new HashMap<String, int[]>(); // refseq, int[] of bin range

        TreeSet<Count> orderedSubjectSet = new TreeSet<Count>(new ResultComparator());
        StringBuilder retval = new StringBuilder();
        FrameBotStatIterator iterator = new FrameBotStatIterator(framebotResult, true);
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            if ( stat.getIdentity() < identity) {
                continue;
            }
            int mappingCount = 1;
            if ( countMap != null) {
                if ( countMap.get(stat.queryID)!= null ){
                    mappingCount = countMap.get(stat.queryID);
                }else {  // not in the mapping, that means we should discard this match
                    continue;
                }
            }            
            int[] identityBin = identityMap.get(stat.subjectID);
            if ( identityBin == null){  // every 5%
                identityBin = new int[21];
                identityMap.put(stat.subjectID, identityBin);
            }
            int binIndex = (int)Math.floor(stat.identity)/5;
            
            identityBin[binIndex] += mappingCount;
           
        }
        iterator.close();
       
        retval.append("sampleName\t" + sampleName +"\n");
        for ( String subject: identityMap.keySet()){
            int[] identityBin = identityMap.get(subject);
            int total = 0;
            for ( int i = 0; i < identityBin.length; i++){
                total += identityBin[i];
            }
            orderedSubjectSet.add(new Count(subject,total) );
        }

        retval.append("Percent Identity\tDescription");
        retval.append("\t" + 100 + "~" );
        for ( int i = 19; i >= 0; i--){
            retval.append("\t" + (i*5 + 4) + "~"+ i*5 );
        }
        retval.append("\t" + "Total\n");


        for ( Count c: orderedSubjectSet){
            if (c.count ==0) continue;
            
            retval.append(c.subject +"\t" + descMap.get(c.subject));
            int[] identityBin = identityMap.get(c.subject);
            for ( int i = identityBin.length -1; i >=0; i--){
                retval.append("\t" + identityBin[i]);
            }
            retval.append("\t" + c.count +"\n");
        }
        return retval.toString();
    }

       
    /**
     * This reads in one FrameBot output, calculates the count for each identity and alignment length, good for graph
     * @param framebotResult
     * @param idMapping
     */
    private static String getOneFramebotSummary(File framebotResult, HashMap<String, Integer> countMap, String sampleName, double identity) throws IOException{
        HashMap<Integer, Integer> alignLenthMap = new HashMap<Integer, Integer>();
        HashMap<Integer, Integer> identityMap = new HashMap<Integer, Integer>();
        HashMap<String, Integer> subjectMap = new HashMap<String, Integer>();

        FrameBotStatIterator iterator = new FrameBotStatIterator(framebotResult, true);
        while (iterator.hasNext()){
            FrameBotStat stat = iterator.next();
            if ( stat.getIdentity() < identity) {
                continue;
            }
            int mappingCount = 1;
            if ( countMap != null && countMap.get(stat.queryID)!= null){
                mappingCount = countMap.get(stat.queryID);
            }
            if ( countMap != null) {
                if ( countMap.get(stat.queryID)!= null ){
                    mappingCount = countMap.get(stat.queryID);
                }else {  // not in the mapping, that means we should discard this match
                    continue;
                }
            }     
            Integer existingSubject = subjectMap.get(stat.subjectID);
            if ( existingSubject == null){
                subjectMap.put(stat.subjectID, mappingCount );
            }else {
                subjectMap.put(stat.subjectID, mappingCount + existingSubject.intValue() );
            }
            Integer existingFrameshift = alignLenthMap.get(stat.alignLen);
            if ( existingFrameshift == null){
                alignLenthMap.put(stat.alignLen, mappingCount );
            }else {
                alignLenthMap.put(stat.alignLen, mappingCount + existingFrameshift.intValue() );
            }
            int rounded_identity = (int)Math.floor(stat.identity);
            Integer existingIdentity = identityMap.get(rounded_identity);
            if ( existingIdentity == null){
                identityMap.put( rounded_identity, mappingCount );
            }else {
                identityMap.put(rounded_identity, mappingCount + existingIdentity.intValue() );
            }
        }
        iterator.close();
        StringBuilder retval = new StringBuilder();
        retval.append("sampleName\t" + sampleName +"\n");
        retval.append("subject\tdescription\tcount\n");
        for ( String n: subjectMap.keySet()){
            retval.append(n +"\t" + descMap.get(n) + "\t" + subjectMap.get(n) +"\n");
        }
        retval.append("");
        retval.append("alignment length\tcount\n");
        TreeSet<Integer> orderSet = new TreeSet(alignLenthMap.keySet());
         for ( Integer n: orderSet){
             retval.append(n.intValue() + "\t" + alignLenthMap.get(n) +"\n");
        }
        retval.append("\n");
        retval.append("percent identity\tcount\n");
        orderSet = new TreeSet(identityMap.keySet());
        for ( Integer n:  orderSet){
            retval.append(n.intValue() + "\t" + identityMap.get(n) +"\n");
        }

        return retval.toString();
    }

    
    public static void main(String[] args) throws IOException {
        String framebotResult = null;
        String idMapping = null;
        String sampleMapping = null;
        String stat = null;
        String outFile = null;
        String descFile = null;
        double identity = 0.0;
        
        try {
            CommandLine line = new PosixParser().parse(options, args);
            if (line.hasOption(ID_MAPPING_SHORT_OPT) ) {
                idMapping = line.getOptionValue(ID_MAPPING_SHORT_OPT);                
            }
                       
            if (line.hasOption(SAMPLE_MAPPING_SHORT_OPT) ) {
                sampleMapping = line.getOptionValue(SAMPLE_MAPPING_SHORT_OPT);                
            } 
            if (line.hasOption(STAT_SHORT_OPT) ) {
                stat = line.getOptionValue(STAT_SHORT_OPT); 
                if ( !(stat.equalsIgnoreCase("stat") || stat.equalsIgnoreCase("hist")
                        || stat.equalsIgnoreCase("summary") || stat.equalsIgnoreCase("matrix")
                        || stat.equalsIgnoreCase("subset") )) {
                    throw new IllegalArgumentException("only stat, hist, summary, matrix and subset are allowed for stat-type");
                }
            } 
            if (line.hasOption(DESCRIPTION_SHORT_OPT) ) {
                descFile = line.getOptionValue(DESCRIPTION_SHORT_OPT);
                descMap = readDesc(new File(descFile));
            } 
            if (line.hasOption(IDENTITY_SHORT_OPT) ) {
                identity = Double.parseDouble(line.getOptionValue(IDENTITY_SHORT_OPT));
                if ( identity < 0 || identity > 100) {
                     throw new IllegalArgumentException("identity cutoff should be in the range of 0 and 100");
                }
            } 
                                
            args = line.getArgs();
            if ( args.length != 2){
                throw new Exception("Incorrect number of command line arguments");
            }
            framebotResult = args[0];
            outFile = args[1];
        } catch (Exception e) {
            System.out.println("Command Error: " + e.getMessage());
            new HelpFormatter().printHelp(120, "GetFrameBotStatMain [options] <FrameBot Alignment file or Dir> <out file>", "", options, "");
            return;
        }
        
         
        if ( stat.equals("matrix")){
            GetFrameBotMatchMatrix.getFramebotMatrix(new File(framebotResult), idMapping, sampleMapping, outFile, identity);
        }else if ( stat.equals("subset")){
            GetFrameBotMatchMatrix.getSubsetMatrix(new File(framebotResult), idMapping, sampleMapping, outFile, identity);
        }else {
            process(new File(framebotResult), idMapping, sampleMapping, outFile, stat, identity);
        }

    }
}
