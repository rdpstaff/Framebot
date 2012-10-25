package edu.msu.cme.rdp.framebot.stat;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author wangqion
 */
public class RdmSelectSampleMapping {
    private static Random random = new Random();
    private static final Options options = new Options();
    
    static {
        options.addOption("n", "num-selection", true, "number of sequence IDs for each sample. Default is the smallest sample size");
        options.addOption("x", "exclude-samples", true, "list of sample names to be excluded from selection");
    }
    
    
    public static void select(String sampleMapping, String outFile, int numofSelection , String excludeSampleFile) throws IOException {
        HashMap<String, HashSet<String>> sampleMap = SampleMappingReader.getSampleMapping(sampleMapping);

        HashSet<String> excludeSampleSet = new HashSet<String> ();
        if ( excludeSampleFile != null){
            excludeSampleSet = getExcludeSampleSet(excludeSampleFile);
        }
        
        int minSeqs = Integer.MAX_VALUE;   // minimum number of seqs per sample
        for ( String sample: sampleMap.keySet()){
            if ( excludeSampleSet.contains(sample.toUpperCase())) continue;
            int size = sampleMap.get(sample).size();
            if (  size < minSeqs){
                minSeqs = size;
            }
            if ( size <= numofSelection){
                throw new IllegalArgumentException("sample " + sample + " has only " + size + " sequences, requested for " + numofSelection);
            }
        }
        
        if (numofSelection == 0){  // if not specify numofSelection, use the minSeqs
            numofSelection = minSeqs;
        }
        PrintStream writer = new PrintStream(outFile);
        for ( String sample: sampleMap.keySet()){
            if ( excludeSampleSet.contains(sample.toUpperCase())) continue;
            ArrayList<String> selectedSeqIDs = selectOneSample(sampleMap.get(sample), numofSelection);
            for ( String id: selectedSeqIDs){
                writer.println(id +"\t" + sample);
            }
        }
        writer.close();
    }
    
    public static ArrayList<String> selectOneSample( HashSet<String> seqIDSet, int numofSelection){
        ArrayList<String> remainingSeqIDs = new ArrayList<String>(seqIDSet);
        ArrayList<String> selectedSeqIDs =  new ArrayList<String>();		

        while ( selectedSeqIDs.size() < numofSelection){
            int rdmIndex = random.nextInt(remainingSeqIDs.size());	
            selectedSeqIDs.add(remainingSeqIDs.get(rdmIndex));
            remainingSeqIDs.remove(rdmIndex);			
        }
        
        return selectedSeqIDs;
    }
    
    
            
    public static HashSet<String> getExcludeSampleSet(String file) throws IOException{
        BufferedReader reader = new BufferedReader(new FileReader(new File(file)));
        String line = null;
        HashSet<String> excludeSampleSet = new HashSet<String>();
        while( (line=reader.readLine()) != null){
            excludeSampleSet.add(line.trim().toUpperCase());
        }
        reader.close();
        return excludeSampleSet;
    }

  

    public static void main(String[] args) throws IOException{
        int numOfSelection = 0;
        String excludeSampleFile = null;
        
        try{
            CommandLine line = new PosixParser().parse(options, args);
            if ( line.hasOption("num-selection")){
                numOfSelection = Integer.parseInt(line.getOptionValue("num-selection"));
                if ( numOfSelection <= 0){
                    throw new Exception("num-selection should be greater than 0");
                }
            }
            if ( line.hasOption("exclude-samples")){
                excludeSampleFile = line.getOptionValue("exclude-samples");
            }
            
            args = line.getArgs();
            if ( args.length != 2){
                throw new Exception("Incorrect number of command line arguments");
            }
            
            RdmSelectSampleMapping.select(args[0], args[1], numOfSelection, excludeSampleFile);
            
        }catch(Exception e){
            new HelpFormatter().printHelp(120, "RdmSelectSampleMapping [options] <sampleMapping> <outfile>", "", options, "");
            System.out.println("ERROR: " + e.getMessage());
            return;
        }
        
              
    }
    
}
