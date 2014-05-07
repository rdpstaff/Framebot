package edu.msu.cme.rdp.framebot.cli;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.framebot.core.ProteinSeqMatch;
import edu.msu.cme.rdp.framebot.core.FramebotCore;
import edu.msu.cme.rdp.framebot.core.FramebotResult;
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
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
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
        options.addOption("k", "knn", true, "The top k closest protein targets. (default = 10)");
        options.addOption("P", "no-prefilter", false, "Disable the pre-filtering step for non-metric search.");
    }

    private static void framebotItUp(List<Sequence> targetSeqs, File queryFile, File qualFile, OutputCoordinator outputCoordinator, AlignmentMode mode, ScoringMatrix simMatrix, int translTable) throws IOException {
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

            outputCoordinator.printResult(bestResult, seq, bestIsReversed);

            System.out.println("Processed " + seq.getSeqName() + " in " + (System.currentTimeMillis() - startTime) + " ms");
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
            OutputCoordinator outputCoordinator, AlignmentMode mode, ScoringMatrix simMatrix, int translTable, int wordSize, int k) throws IOException {
        SeqReader queryReader;

        if(qualFile != null) {
            queryReader = new QSeqReader(queryFile, qualFile);
        } else {
            queryReader = new SequenceReader(queryFile);
        }

        Sequence seq;
        int seqCount = 0;
        long totalTime = System.currentTimeMillis();                
        ProteinSeqMatch theObj = new ProteinSeqMatch(targetSeqs, wordSize);            

        while ((seq = queryReader.readNextSequence()) != null) {
            Sequence reverseComplement = new Sequence(seq.getSeqName(), "", IUBUtilities.reverseComplement(seq.getSeqString()));

            FramebotResult bestResult = null;
            boolean bestIsReversed = true;
            
            long startTime = System.currentTimeMillis();
            ArrayList<ProteinSeqMatch.BestMatch> topKMatches= theObj.findTopKMatch(seq, k);
             
            for (ProteinSeqMatch.BestMatch bestTarget : topKMatches) {
                FramebotResult result;
                if ( !bestTarget.isRevComp()){
                    result = FramebotCore.processSequence(seq, bestTarget.getBestMatch(), true, translTable, mode, simMatrix);
                } else {
                    result = FramebotCore.processSequence(reverseComplement, bestTarget.getBestMatch(), true, translTable, mode, simMatrix);
                }

                //System.err.println(seq.getSeqName() + " target=" + bestTarget.getBestMatch().getSeqName() + " score=" + result.getFrameScore().getMaxScore());
                if (bestResult == null || bestResult.getFrameScore().getMaxScore() < result.getFrameScore().getMaxScore()) {
                    bestResult = result;
                    bestIsReversed = bestTarget.isRevComp();
                }
            }

            outputCoordinator.printResult(bestResult, seq, bestIsReversed);

            System.out.println("Processed " + seq.getSeqName() + " in " + (System.currentTimeMillis() - startTime) + " ms");
            seqCount++;
        }

        System.out.println("Processed " + seqCount + " sequences in " + (System.currentTimeMillis() - totalTime) + " ms");
    }

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

    public static void main(String[] args) throws IOException {
        File indexOrSeedFile = null;
        File queryFile = null;
        File qualFile = null;
        int lengthCutoff = 80;
        double identityCutoff = .4;
        String outputStem = null;
        boolean metricSearch = true;
        boolean usePrefilter = true;
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
            if (line.hasOption("no-prefilter")) {
                usePrefilter = false;
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
                if ( k < 1 ){
                    throw new Exception("knn must be at least 1");
                }
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
                if ( usePrefilter){
                     FramebotMain.framebotItUp_prefilter(targetSeqs, queryFile, qualFile, outputCoordinator, mode, scoringMatrix, translTable, wordSize, k);
                } else {
                    FramebotMain.framebotItUp(targetSeqs, queryFile, qualFile, outputCoordinator, mode, scoringMatrix, translTable);
                }
                
            }
        } finally {
            outputCoordinator.close();
        }
    }
}
