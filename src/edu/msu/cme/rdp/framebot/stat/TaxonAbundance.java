/*
 * Copyright (C) 2015 Michigan State University <rdpstaff at msu.edu>
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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 *
 * @author wangqion @msu.edu
 */
public class TaxonAbundance {
    private static final String dformat = "%1$.4f";
    private static Options options = new Options();
    static {
        options.addOption(new Option("c", "seqCoverage", true, "contains the ID and coverage separated by space or tab." +
               " Used to adjust the sequence abundance"));
        options.addOption(new Option("e", "identity", true, "the minimum protein identity, default is 0, range [0-100]"));
    }
    
    public static HashMap<String, Double> parseKmerCoverage(String infile) throws IOException{
        HashMap<String, Double> countMap = new HashMap<String, Double>();
        BufferedReader reader = new BufferedReader(new FileReader(new File(infile)));
        String line ;
        while ( (line = reader.readLine() )!= null) {
            if ( line.startsWith("#")) continue;
            String [] val = line.split("\\s+");
            countMap.put(val[0], Double.parseDouble(val[1]));
        }
        reader.close();
        return countMap;
    }
    
    
    public static void mapAbundance( File framebotResult, File lineagefile, String outfile, HashMap<String, Double> coveragetMap, double identity) throws IOException{
        HashMap<String,String> lineageMap = GetFrameBotStatMain.readDesc(lineagefile);
        HashMap<String, Double> rankMatchMap = new HashMap<String, Double>(); // taxon rank, 
        HashMap<String, Double> matchMap = new HashMap<String, Double>();  // match name
        
        double totalCount = 0.0;
        File[] files ;
        if ( framebotResult.isDirectory()){
            files = framebotResult.listFiles();
        }else {
            files = new File[1];
            files[0] = framebotResult;
        }
        String line ;
        BufferedReader reader ;
        for ( File f: files){
            reader = new BufferedReader(new FileReader(f));
            while ( (line = reader.readLine() )!= null) {
                if ( !line.startsWith("STAT")) continue;
                String[] tokens = line.trim().split("\\t");
                double thisIdentity = Double.parseDouble(tokens[5]);
                if ( thisIdentity < identity) {
                    continue;
                }
                
                String match = tokens[1];

                String[] lineage = lineageMap.get(match).split(";");
                String taxonname = lineage[0].trim();
                if ( lineage.length > 1){
                    taxonname = lineage[1];
                    if ( lineage[1].equalsIgnoreCase("Proteobacteria")){// || lineage[1].equalsIgnoreCase("Bacteroidetes")){
                        taxonname = lineage[2].trim();
                    }
                }
                if ( taxonname.equals("")){
                    taxonname = "NA";
                }

                double count = 1.0;
                if ( coveragetMap != null){
                    count = coveragetMap.get(tokens[2]);
                }
                Double prevCount = rankMatchMap.get(taxonname);
                if ( prevCount != null){
                    rankMatchMap.put(taxonname, count + prevCount );
                }else {
                    rankMatchMap.put(taxonname, count );
                }
                Double prevMatchCount = matchMap.get(match);
                if ( prevMatchCount != null){
                    matchMap.put(match, count + prevMatchCount );
                }else {
                    matchMap.put(match, count );
                }
                totalCount += count;
            }
            reader.close();
        }
        
        PrintStream out = new PrintStream(new File(outfile));
        // print out abundance group by phylum or class (within Proteobacteria)
        out.println("Taxon\tAbundance\tFraction Abundance");
        for ( String taxonname: rankMatchMap.keySet()){
            out.println( taxonname + "\t" + rankMatchMap.get(taxonname) + "\t" + String.format(dformat, rankMatchMap.get(taxonname)/totalCount));
        }
        
        // print out abundance by closest match   
        out.println("\n\nLineage\tMatchName\tAbundance\tFraction Abundance");
        for ( String match: matchMap.keySet()){
            String lineage = lineageMap.get(match);   
            if ( lineage.equals("")){
                lineage = "NA";
            }
            out.println( lineage + "\t" + match + "\t" + matchMap.get(match) + "\t" + String.format(dformat, matchMap.get(match)/totalCount));
        }
        
        out.close();
    }
    
    /**
     * this class group the nearest matches by phylum/class, or by match 
     * @param args
     * @throws Exception 
     */
    public static void main (String[] args) throws Exception {
        HashMap<String, Double> coveragetMap = null;
        double identity = 0.0;
        try {
            CommandLine line = new PosixParser().parse(options, args);
            if (line.hasOption("seqCoverage") ) {
                String coveragefile = line.getOptionValue("seqCoverage");   
                coveragetMap = parseKmerCoverage(coveragefile);
            }
            if (line.hasOption("identity") ) {
                identity = Double.parseDouble(line.getOptionValue("identity"));
                if ( identity < 0 || identity > 100) {
                     throw new IllegalArgumentException("identity cutoff should be in the range of 0 and 100");
                }
            } 
                 
            args = line.getArgs();
            if (args.length != 3) {
                throw new Exception("");
            }
            
        }catch (Exception e) {
            System.out.println("Command Error: " + e.getMessage());
            new HelpFormatter().printHelp(80, "[options] <FrameBot Alignment file or Dir> <seqLineage> <out file> ", "", options, 
                    "seqLineage: a tab-delimited file with ref seqID and lineage, or fasta of ref seq with lineage as the descrption"
                    + "\nframeBot alignment file or Dir: frameBot alignment files "
                    + "\noutfile: output with the nearest match count group by phylum/class; and by match name"
            );
        }
        TaxonAbundance.mapAbundance(new File(args[0]), new File(args[1]), args[2], coveragetMap, identity);
    }
     
}
