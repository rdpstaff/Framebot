package edu.msu.cme.rdp.framebot.cli;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAligner;
import edu.msu.cme.rdp.alignment.pairwise.PairwiseAlignment;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.alignment.pairwise.rna.DistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.IdentityDistanceModel;
import edu.msu.cme.rdp.alignment.pairwise.rna.OverlapCheckFailedException;
import edu.msu.cme.rdp.framebot.core.FramebotCore;
import edu.msu.cme.rdp.framebot.core.FramebotResult;
import edu.msu.cme.rdp.framebot.core.FramebotScoringMatrix;
import edu.msu.cme.rdp.framebot.index.FramebotIndex;
import edu.msu.cme.rdp.framebot.index.FramebotIndex.FramebotSearchResult;
import edu.msu.cme.rdp.framebot.output.OutputCoordinator;
import edu.msu.cme.rdp.readseq.readers.QSeqReader;
import java.io.File;
import java.io.IOException;

import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.readseq.utils.kmermatch.ProteinSeqMatch;
import edu.msu.cme.rdp.readseq.utils.orientation.ProteinWordGenerator;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

public class FramebotMain {

    private static final Options options = new Options();
    private static final int DENOVO_ABUND_LIMIT = 10; // minimum size of unique seqs to be considered as abundant
    private static final double DENOVO_ID_CUTOFF = 0.7;
    private static int denovoAbundanceLimit = DENOVO_ABUND_LIMIT;  // minimum size of unique seqs to be considered as abundant
    private static double denovoIdentityCutoff = DENOVO_ID_CUTOFF;
    private static final ProteinUtils proteinUtils = ProteinUtils.getInstance();
    private static final String delims = "[=;]";
    private static DistanceModel dist = new IdentityDistanceModel();

    static {
        options.addOption("o", "result-stem", true, "Result file name stem (default=stem of query nucl file)");
        options.addOption("N", "no-metric-search", false, "Disable metric search (provide fasta file of seeds instead of index file)");
        options.addOption("i", "identity-cutoff", true, "Percent identity cutoff [0-1] (default = .4)");
        options.addOption("l", "length-cutoff", true, "Length cutoff in number of amino acids (default = 80)");
        options.addOption("a", "alignment-mode", true, "Alignment mode: glocal, local or global (default = glocal)");
        options.addOption("q", "quality-file", true, "Sequence quality data");
        options.addOption("m", "max-radius", true, "maximum radius for metric-search ONLY, range [1-" + Integer.MAX_VALUE + "], default uses the maxRadius specified in the index file");
        options.addOption("x", "scoring-matrix", true, "the protein scoring matrix for no-metric-search ONLY. Default is Blosum62");
        options.addOption("g", "gap-open-penalty", true, "gap opening penalty for no-metric-search ONLY. Default is " + ScoringMatrix.DEFAULT_GAP_OPEN_PEALTY);
        options.addOption("e", "gap-ext-penalty", true, "gap extension penalty for no-metric-search ONLY. Default is " + ScoringMatrix.DEFAULT_GAP_EXT_PENALTY);
        options.addOption("f", "frameshift-penalty", true, "frameshift penalty for no-metric-search ONLY. Default is " + ScoringMatrix.DEFAULT_FRAME_SHIFT_PENALTY);
        
        options.addOption("t", "transl-table", true, "Protein translation table to use (integer based on ncbi's translation tables, default=11 bacteria/archaea)");
        options.addOption("w", "word-size", true, "The word size used to find closest protein targets. (default = " + ProteinWordGenerator.WORDSIZE + ", recommended range [3 - 6])");
        options.addOption("k", "knn", true, "The top k closest targets from kmer prefilter step. Set k=0 to disable this step. (default = 10)");
        
        options.addOption("z", "de-novo", false, "Enable de novo mode to add abundant query seqs to refset. Only works for no-metric-search");
        options.addOption("b", "denovo-abund-cutoff", true, "minimum abundance for de-novo mode. default = " + DENOVO_ABUND_LIMIT);
        options.addOption("d", "denovo-id-cutoff", true, "maxmimum aa identity cutoff for query to be added to refset for de-novo mode. default = " + DENOVO_ID_CUTOFF);
    }

    private static void framebotItUp(List<Sequence> targetSeqs, File queryFile, File qualFile, OutputCoordinator outputCoordinator, 
            AlignmentMode mode, ScoringMatrix simMatrix, int translTable, boolean denovo) throws IOException, OverlapCheckFailedException {
        SeqReader queryReader;

        if(qualFile != null) {
            queryReader = new QSeqReader(queryFile, qualFile);
        } else {
            queryReader = new SequenceReader(queryFile);
        }

        Sequence seq;
        int seqCount = 0;
        long totalTime = System.currentTimeMillis();

        while ((seq = queryReader.readNextSequence()) != null) {
            Sequence reverseComplement = new Sequence(seq.getSeqName(), "", IUBUtilities.reverseComplement(seq.getSeqString()));

            FramebotResult bestResult = null;
            boolean bestIsReversed = true;

            long startTime = System.currentTimeMillis();
            for (Sequence targetSeq : targetSeqs) {
                FramebotResult resultForward = FramebotCore.processSequence(seq, targetSeq, true, translTable, mode, simMatrix);
                FramebotResult resultReverse = FramebotCore.processSequence(reverseComplement, targetSeq, true, translTable, mode, simMatrix);

                FramebotResult result;
                if (resultForward.getFrameScore().getMaxScore() < resultReverse.getFrameScore().getMaxScore()) {
                    result = resultReverse;
                } else {
                    result = resultForward;
                }

                if (bestResult == null || bestResult.getFrameScore().getMaxScore() < result.getFrameScore().getMaxScore()) {
                    bestResult = result;
                    bestIsReversed = (result == resultReverse);
                }
            }

            if ( !denovo) {
                outputCoordinator.printResult(bestResult, seq, bestIsReversed);
            } else {
                Sequence targetProSeq = new Sequence( bestResult.getAlignedTarget().getSeqName(), "", SeqUtils.getUnalignedSeqString(bestResult.getAlignedTarget().getSeqString()));
                String protseqToadd = checkDenovo(seq, targetProSeq, bestResult, outputCoordinator, translTable);
                if ( protseqToadd != null){ // we need to redo the frameshift calculation
                    targetSeqs.add(new Sequence(seq.getSeqName(), "", protseqToadd));
                    Sequence tempSeq = new Sequence(seq.getSeqName(), "", protseqToadd);
                    if ( !bestIsReversed){
                        bestResult = FramebotCore.processSequence(seq, tempSeq, true, translTable, mode, simMatrix);
                    }else {
                        bestResult = FramebotCore.processSequence(reverseComplement, tempSeq, true, translTable, mode, simMatrix);
                    }
                }
                outputCoordinator.printResult(bestResult, seq, bestIsReversed);
            }
            
            seqCount++;
        }

        System.out.println("Processed " + seqCount + " sequences in " + (System.currentTimeMillis() - totalTime) + " ms");
    }
    
    /**
     * This method searches the closest target sequences first using the Seqmatch algorithm with protein sequence, then do frameshift correction against the top k target sequences
     * @param targetSeqs
     * @param queryFile
     * @param qualFile
     * @param outputCoordinator
     * @param mode
     * @param simMatrix
     * @param translTable
     * @throws IOException 
     */
    private static void framebotItUp_prefilter(List<Sequence> targetSeqs, File queryFile, File qualFile, 
            OutputCoordinator outputCoordinator, AlignmentMode mode, ScoringMatrix simMatrix, int translTable, int wordSize, int k, boolean denovo) throws IOException, OverlapCheckFailedException {
        SeqReader queryReader = new SequenceReader(queryFile);

        if(qualFile != null) {
            queryReader = new QSeqReader(queryFile, qualFile);
        } else {
            queryReader = new SequenceReader(queryFile);
        }

        Sequence seq;
        ProteinSeqMatch theObj = new ProteinSeqMatch(targetSeqs, wordSize);            

        while ((seq = queryReader.readNextSequence()) != null) {
            Sequence reverseComplement = new Sequence(seq.getSeqName(), "", IUBUtilities.reverseComplement(seq.getSeqString()));

            FramebotResult bestResult = null;
            boolean bestIsReversed = true;               
            ArrayList<ProteinSeqMatch.BestMatch> topKMatches= theObj.findTopKMatch(seq, k);
             
            for (ProteinSeqMatch.BestMatch bestTarget : topKMatches) {
                
                FramebotResult result;
                if ( !bestTarget.isRevComp()){
                    result = FramebotCore.processSequence(seq, bestTarget.getBestMatch(), true, translTable, mode, simMatrix);
                } else {
                    result = FramebotCore.processSequence(reverseComplement, bestTarget.getBestMatch(), true, translTable, mode, simMatrix);
                }

                if (bestResult == null || bestResult.getFrameScore().getMaxScore() < result.getFrameScore().getMaxScore()) {
                    bestResult = result;
                    bestIsReversed = bestTarget.isRevComp();
                }
            }

            if ( !denovo) {
                outputCoordinator.printResult(bestResult, seq, bestIsReversed);
            } else {
                Sequence targetProSeq = theObj.getRefSeq(bestResult.getAlignedTarget().getSeqName());
                String protseqToadd = checkDenovo(seq, targetProSeq, bestResult, outputCoordinator, translTable);
                if ( protseqToadd != null){ // we need to redo the frameshift calculation
                    theObj.addRefSeq(new Sequence(seq.getSeqName(), "", protseqToadd));
                    Sequence tempSeq = new Sequence(seq.getSeqName(), "", protseqToadd);
                    if ( !bestIsReversed){
                        bestResult = FramebotCore.processSequence(seq, tempSeq, true, translTable, mode, simMatrix);
                    }else {
                        bestResult = FramebotCore.processSequence(reverseComplement, tempSeq, true, translTable, mode, simMatrix);
                    }
                }
                outputCoordinator.printResult(bestResult, seq, bestIsReversed);
            }                
        }
    }
    

    private static String checkDenovo( Sequence seq, Sequence targetProSeq, FramebotResult bestResult, 
            OutputCoordinator outputCoordinator, int translTable )throws IOException, OverlapCheckFailedException { 
        // we don't need to get more de novo refs if above certain pct Identity, or below the minmum cutoff
        if ( bestResult.getPercentIdent() < outputCoordinator.getIdentityCutoff() || bestResult.getPercentIdent() >= denovoIdentityCutoff) { 
            return null;
        }  

        //parsing description to get size,              
        String[] elements = seq.getDesc().split(delims);
        int size = Integer.parseInt(elements[elements.length - 1]);
        //abundance check an dpctID check
        if(size < denovoAbundanceLimit ){
            return null;
        }

        Boolean addToRefSet = false;
        String corr_prot_seqstring =  SeqUtils.getUnalignedSeqString(bestResult.getAlignedQuery().getSeqString());
        String bestProSeq = corr_prot_seqstring;
        if(bestResult.getFrameshifts() == 0){                    
            //if do not contain stop codons
            if(corr_prot_seqstring.indexOf('*') == -1){                
                addToRefSet = true;
            }                    
        } else {
            PairwiseAlignment upw = null;                    
            double bestUpwIdent = 0.0;            
            for (int i = 0; i < FramebotScoringMatrix.NUM_FRAMES; i++) {
                String proSeqString = proteinUtils.translateToProtein(seq.getSeqString().substring(i), true, translTable);
                if(proSeqString.indexOf('*') != -1){
                    continue;
                }

                PairwiseAlignment temp_upw = PairwiseAligner.align(targetProSeq.getSeqString(), proSeqString, ScoringMatrix.getDefaultProteinMatrix(), AlignmentMode.glocal);
                double ident = 1 - dist.getDistance(temp_upw.getAlignedSeqi().getBytes(), temp_upw.getAlignedSeqj().getBytes(), 0);
                if(upw == null || ident > bestUpwIdent){
                    upw = temp_upw;
                    bestProSeq = proSeqString;
                    bestUpwIdent = ident;
                }
            }
            if(upw != null){   
                addToRefSet = true;                                   
            }
        }    

        if ( addToRefSet){
            return bestProSeq;
        }else {
            return null;
        }
        
    }

    /**
     * This method should be obsolete now, replace by the pre-filter approach
     * @param index
     * @param maxRadius
     * @param queryFile
     * @param qualFile
     * @param outputCoordinator
     * @throws IOException 
     */
    private static void framebotItUp(FramebotIndex index, int maxRadius, File queryFile, File qualFile, OutputCoordinator outputCoordinator) throws IOException {SeqReader queryReader;

        if(qualFile != null) {
            queryReader = new QSeqReader(queryFile, qualFile);
        } else {
            queryReader = new SequenceReader(queryFile);
        }

        Sequence seq;
        int seqCount = 0;
        long totalTime = System.currentTimeMillis();

        System.out.println("query_id\tseed_id\t%identity\tdistance\treversed?\tseeds_examined\tratio_seeds_examined\ttime(ms)");
        while ((seq = queryReader.readNextSequence()) != null) {
            long startTime = System.currentTimeMillis();

            FramebotSearchResult resultForward = index.getNearestNeighbor(seq, maxRadius);           
            FramebotSearchResult result = resultForward;
                           
            if ( result != null){
                outputCoordinator.printResult(result.getResult(), seq, result.isReversed());

                System.out.println(seq.getSeqName() + "\t"
                    + result.getSeedId() + "\t"
                    + result.getResult().getPercentIdent() + "\t"
                    + result.getDist() + "\t"
                    + result.isReversed() + "\t"
                    + result.getNumComparisons() + "\t"
                    + ((float) result.getNumComparisons() / (index.getSeedCount() * 2)) + "\t"
                    + (System.currentTimeMillis() - startTime));
            }else {
                System.out.println(seq.getSeqName() + "\t" + "Failed to find a nearest neighbor within the radius:" + maxRadius + ". You may need to increase maximum radius.");
                
            }
            seqCount++;
        }

        System.out.println("Processed " + seqCount + " sequences in " + (System.currentTimeMillis() - totalTime) + " ms");
    }

    public static void main(String[] args) throws IOException, OverlapCheckFailedException {
        File indexOrSeedFile = null;
        File queryFile = null;
        File qualFile = null;
        int lengthCutoff = 80;
        double identityCutoff = .4;
        String outputStem = null;
        boolean metricSearch = true;
        int maxRadius = 0;
        int translTable = 11;
        int wordSize = ProteinWordGenerator.WORDSIZE;
        int k = 10;
        AlignmentMode mode = AlignmentMode.glocal;
        File scoringMatrixFile = null;
        int gapPenalty = ScoringMatrix.DEFAULT_GAP_OPEN_PEALTY;
        int gapExtend = ScoringMatrix.DEFAULT_GAP_EXT_PENALTY;
        int frameshiftPenalty = ScoringMatrix.DEFAULT_FRAME_SHIFT_PENALTY;
        boolean useDefaultMatrixParameter = true; 
        boolean denovoMode = false;
        

        File framebotOutput;
        File failedFramebotOutput;
        File correctedProtFile;
        File correctedNuclFile;
        File correctedNuclQualFile = null;
        File failedNuclFile;

        try {
            CommandLine line = new PosixParser().parse(options, args);

            if (line.hasOption("result-stem")) {
                outputStem = line.getOptionValue("result-stem");
            }
            if (line.hasOption("identity-cutoff")) {
                identityCutoff = Double.parseDouble(line.getOptionValue("identity-cutoff"));

                if (identityCutoff < 0 || identityCutoff > 1) {
                    throw new Exception("Identity cutoff must be in the range [0-1]");
                }
            }
            if (line.hasOption("length-cutoff")) {
                lengthCutoff = Integer.parseInt(line.getOptionValue("length-cutoff"));
            }            
            if (line.hasOption("no-metric-search")) {
                metricSearch = false;
            }
            if (line.hasOption("denovo-abund-cutoff")) {
                 denovoAbundanceLimit = Integer.parseInt(line.getOptionValue("denovo-abund-cutoff"));
                 if ( denovoAbundanceLimit < 1){
                     throw new Exception("abund-cutoff must be at least 1");
                 }
            }
            if (line.hasOption("denovo-id-cutoff")) {
                 denovoIdentityCutoff = Double.parseDouble(line.getOptionValue("denovo-id-cutoff"));
                 if ( denovoIdentityCutoff < identityCutoff || denovoIdentityCutoff > 1){
                     throw new Exception("denovo-id-cutoff should be between " + identityCutoff + " and 1.");
                 }
            }
            if (line.hasOption("scoring-matrix")) {
                scoringMatrixFile = new File(line.getOptionValue("scoring-matrix"));
            }
            if (line.hasOption("gap-open-penalty")) {
                gapPenalty = Integer.parseInt(line.getOptionValue("gap-open-penalty"));
                useDefaultMatrixParameter = false;
            }
            if (line.hasOption("gap-ext-penalty")) {
                gapExtend = Integer.parseInt(line.getOptionValue("gap-ext-penalty"));
                useDefaultMatrixParameter = false;
            }
            if (line.hasOption("frameshift-penalty")) {
                frameshiftPenalty = Integer.parseInt(line.getOptionValue("frameshift-penalty"));
                useDefaultMatrixParameter = false;
            } 
            
            if (line.hasOption("max-radius")) {
                maxRadius = Integer.parseInt(line.getOptionValue("max-radius"));
                if(maxRadius < 1) {
                    throw new Exception("Max radius must be at least 1");
                }
            }

            if (line.hasOption("quality-file")) {
                qualFile = new File(line.getOptionValue("quality-file"));
            }

            if (line.hasOption("alignment-mode")) {
                String m = line.getOptionValue("alignment-mode");
                if (m.equalsIgnoreCase(AlignmentMode.glocal.name())) {
                    mode = AlignmentMode.glocal;
                } else if (m.equalsIgnoreCase(AlignmentMode.local.name())) {
                    mode = AlignmentMode.local;
                } else if (m.equalsIgnoreCase(AlignmentMode.global.name())) {
                    mode = AlignmentMode.global;
                } else {
                    throw new IllegalArgumentException("Only " + AlignmentMode.glocal.name() + " or " + AlignmentMode.local.name() + " or " + AlignmentMode.global.name() + " are allowed for the alignment type");
                }
            }
            if (line.hasOption("word-size")) {
                wordSize = Integer.parseInt(line.getOptionValue("word-size"));
                if ( wordSize < 3 ){
                    throw new Exception("Word size must be at least 3");
                }
            }
            if (line.hasOption("knn")) {
                k = Integer.parseInt(line.getOptionValue("knn"));
                // if k==0 means no prfilter
                if ( k < 0 ){
                    throw new Exception("knn must be at least 0. Prefilter step find the top k closest targets based on kmer matching. Set k=0 to disable the kmer prefilter step");
                }
            }
            
            if (line.hasOption("de-novo")) {
                denovoMode = true;
            }

            args = line.getArgs();

            if (outputStem == null && (args.length != 7 || args.length != 8)) {
                throw new Exception("No result-stem specified, must provide result-stem or <framebot_out.txt> <failed_framebot_out.txt> <corr_prot_out.fasta> <corr_nucl_out.fasta> <failed_nucl_out.fasta> [qual_out_file]");
            } else if(outputStem != null && args.length != 2) {
                throw new Exception("Incorrect number of command line arguments");
            }

            indexOrSeedFile = new File(args[0]);
            queryFile = new File(args[1]);

            if (!indexOrSeedFile.exists()) {
                throw new Exception("Index or seed file does not exist");
            }

            if (!queryFile.exists()) {
                throw new Exception("Query file does not exist");
            }
            if ( metricSearch ){
                if ( !useDefaultMatrixParameter){
                    throw new Exception("gap-open-penalty, gap-ext-penalty or frameshift-penalty are not allowed to change for metric search");
                }
                if( scoringMatrixFile != null){
                    throw new Exception( "Scoring matrix is not allowed to change for metric search. Rebuild the index file with the scoring matrix");
                }
            }
            if (line.hasOption("transl-table")) {
                translTable = Integer.parseInt(line.getOptionValue("transl-table"));
            }

            if (outputStem == null) {
                framebotOutput = new File(args[2]);
                failedFramebotOutput = new File(args[3]);
                correctedProtFile = new File(args[4]);
                correctedNuclFile = new File(args[5]);
                failedNuclFile = new File(args[6]);
                if(args.length == 8) {
                    correctedNuclQualFile = new File(args[7]);
                }
            } else {
                framebotOutput = new File(outputStem + "_framebot.txt");
                failedFramebotOutput = new File(outputStem + "_failed_framebot.txt");
                correctedProtFile = new File(outputStem + "_corr_prot.fasta");
                correctedNuclFile = new File(outputStem + "_corr_nucl.fasta");
                failedNuclFile = new File(outputStem + "_failed_nucl.fasta");
                correctedNuclQualFile = new File(outputStem + "_corr_nucl.qual");
            }

        } catch (Exception e) {
            new HelpFormatter().printHelp(120, "FramebotMain [options] <seed or index file> <query file>", "", options, "");
            System.out.println("ERROR: " + e.getMessage());
            return;
        }

        
        OutputCoordinator outputCoordinator = new OutputCoordinator(lengthCutoff, identityCutoff, framebotOutput, failedFramebotOutput, correctedProtFile, correctedNuclFile, correctedNuclQualFile, failedNuclFile);

        try {
            
            if (metricSearch) {
                FramebotIndex index = FramebotIndex.readExternalIndex(indexOrSeedFile, mode);
                if ( maxRadius == 0){
                    maxRadius = index.getMaxRadius();
                }
                FramebotMain.framebotItUp(index, maxRadius, queryFile, qualFile, outputCoordinator);
            } else {
                List<Sequence> targetSeqs = new ArrayList();
                for (Sequence seq : SequenceReader.readFully(indexOrSeedFile)) {
                    if (seq.getSeqName().startsWith("#")) { //Metasequence
                        continue;
                    }
                    targetSeqs.add(SeqUtils.getUnalignedSeq(seq));
                }
                ScoringMatrix scoringMatrix = null;
                if ( scoringMatrixFile!= null){
                    scoringMatrix = new ScoringMatrix(new FileInputStream(scoringMatrixFile), gapPenalty, gapExtend, frameshiftPenalty);
                }else if (!useDefaultMatrixParameter){
                    scoringMatrix = new ScoringMatrix(ScoringMatrix.getDefaultProteinMatrixStream(), gapPenalty, gapExtend, frameshiftPenalty);
                }else {
                    scoringMatrix = ScoringMatrix.getDefaultProteinMatrix();
                }
                if ( k > 0){
                     FramebotMain.framebotItUp_prefilter(targetSeqs, queryFile, qualFile, outputCoordinator, mode, scoringMatrix, translTable, wordSize, k, denovoMode);
                } else {
                    FramebotMain.framebotItUp(targetSeqs, queryFile, qualFile, outputCoordinator, mode, scoringMatrix, translTable, denovoMode);
                }
                
            }
        } finally {
            outputCoordinator.close();
        }
    }
}
