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

package edu.msu.cme.rdp.framebot.core;

import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.utils.orientation.ProteinWordGenerator;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;

/**
 *
 * @author wangqion
 */
public class ProteinSeqMatch {
    
    private HashMap<String, HashSet<String>> refWordMap = new HashMap<String, HashSet<String>>();
    private HashMap<String, Sequence> refSeqMap = new HashMap<String, Sequence>();
    private ProteinWordGenerator proteinWordGenerator = null;
    private static final float SabThreshold = 0.6f; //
    
    public static class ResultComparator implements Comparator {
        public int compare(Object a, Object b){
            if ( ((BestMatch)a).sab < ((BestMatch)b).sab){
                return 1;
            }else {
                if ( ((BestMatch)a).sab > ((BestMatch)b).sab){
                    return -1;
                }else {
                    return ((BestMatch)a).match.getSeqName().compareTo(((BestMatch)b).match.getSeqName());
                }
            }
        }
    }
    
    public static class BestMatch{
        Sequence match;
        float sab;
        boolean revComp ;
        
        BestMatch(Sequence s, float score, boolean r){
            match = s;
            sab = score;
            revComp = r;
        }
        
        public Sequence getBestMatch(){
            return match;
        }
        
        public boolean isRevComp(){
            return revComp;
        }
        
        public float getSab(){
            return sab;
        }
        
    }
    
    public ProteinSeqMatch(List<Sequence> refSeqs, int wordSize){
        proteinWordGenerator = new ProteinWordGenerator(wordSize);
        // initialize the protein words 
        for (Sequence seq: refSeqs){
            refSeqMap.put(seq.getSeqName(), seq);
            refWordMap.put(seq.getSeqName(), proteinWordGenerator.parseProtein(seq.getSeqString()) );
        }
    }
    
    
    /*
    * This program take a nucleotide sequence, translates to three protein sequence strings, one for each frame.
    * It then calclates the sab score (shared protein kmers) between query and each reference seq.
    * The above process is repeated for the reverse orientation. 
    * This returns the top k best matching reference sequences 
    */
    public ArrayList<BestMatch> findTopKMatch(Sequence nuclSeq, int k){
        TreeSet<BestMatch> orderedResultSet = new TreeSet<BestMatch>( new ResultComparator());
        HashSet<String> queryWordSet = proteinWordGenerator.parseNuclAllFrames(nuclSeq.getSeqString());
        // we need to divide the query word size by three since we have words from three frames
        float queryWordSize = queryWordSet.size()/3f;
        float tempBestSab = 0;
        HashSet<String> targetWordSet ;
        for ( Sequence target: refSeqMap.values()){
            targetWordSet = refWordMap.get(target.getSeqName());
            HashSet<String> tempSet = new HashSet();
            tempSet.addAll(queryWordSet);
            float minWordCount =  queryWordSize<= targetWordSet.size()? queryWordSize: targetWordSet.size();
            tempSet.retainAll(targetWordSet);
            float sab = tempSet.size()/minWordCount;
            if ( sab >= tempBestSab){
                tempBestSab = sab;
            }
            orderedResultSet.add( new BestMatch(target, sab, false));
        }
        
        if(tempBestSab < this.SabThreshold) {  // check reverse
            queryWordSet = proteinWordGenerator.parseNuclAllFrames( IUBUtilities.reverseComplement(nuclSeq.getSeqString()) );
            for ( Sequence target: refSeqMap.values()){
                targetWordSet = refWordMap.get(target.getSeqName());
                HashSet<String> tempSet = new HashSet();
                tempSet.addAll(queryWordSet);
                float minWordCount =  queryWordSize<= targetWordSet.size()? queryWordSize: targetWordSet.size();
                tempSet.retainAll(targetWordSet);
                float sab = tempSet.size()/minWordCount;                
                orderedResultSet.add( new BestMatch(target, sab, true));
            }
        }
        
        ArrayList<BestMatch> topkMatchList = new ArrayList<BestMatch>();
        for ( BestMatch m: orderedResultSet){
            if ( topkMatchList.size() < k)
                topkMatchList.add(m);
        }
        
        return topkMatchList;        
    }
   
    public static void main (String[] args ) throws IOException{
        String usage = "protein_ref.fa query_nucl.fa word_size k\n word_size 4 is recommended for the best performance. Range from 3 to 6.";
        
        ArrayList<Sequence> targetSeqs = new ArrayList<Sequence>();
        for (Sequence seq : SequenceReader.readFully(new File(args[0]))) {           
            targetSeqs.add(SeqUtils.getUnalignedSeq(seq));
        }
        
        int wordSize = Integer.parseInt(args[2]);
        int k = Integer.parseInt(args[3]);
        ProteinSeqMatch theObj = new ProteinSeqMatch(targetSeqs, wordSize);
        
        SequenceReader queryReader = new SequenceReader(new File(args[1]) );
        Sequence nucl= null;
        while ( ( nucl = queryReader.readNextSequence()) != null){
            ArrayList<BestMatch> results = theObj.findTopKMatch(nucl, k);
            for ( BestMatch m: results){
                System.out.println(nucl.getSeqName() + "\t" + m.getBestMatch().getSeqName() + "\t" + m.getSab());
            }
        }
        queryReader.close();
    }
}
