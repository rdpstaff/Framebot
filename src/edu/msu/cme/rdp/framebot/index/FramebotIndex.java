/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.framebot.index;

import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.framebot.core.FramebotCore;
import edu.msu.cme.rdp.framebot.core.FramebotResult;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

/**
 * So this class is mildly complicated because of the parallel arrays and some assumptions made
 * These assumptions if broken will HORRIBLY HORRIBLY MANGLE ANY RESULTS
 * So don't edit the xml index file by hand.  
 * 
 * So the index stores all sorts of information about how/when it was compiled, the alignment mode and the like
 * 
 * In order to improve efficency the index contains the provided seeds and their reverse complements
 * this allows us to search for both forward and reversed at the same time resulting in a modest speed up
 * 
 * It should use the global mode and the same scoring matrix for building index and actually search. 
 * This does add complexity though, so every index in the distance matrix corresponds to the sequence in the same location in the seedSeqs array
 * Also, seedCount is the raw number of seeds, but the total number of seed sequences is seedCount * 2 (the original seeds plus their reverse complement)
 * the array is maintained in such a way that every seed in the array after seedCount is the reverse complement of seedSeq[index - seedCount]
 * 
 * Finally, the metric search is an AESA algorithm and initial tests show that on average only 25% of seed sequences have to be consulted to find 
 * the nearest neighbor
 *
 * @author fishjord
 */
@XmlRootElement
public class FramebotIndex {

    @XmlAttribute
    private int translTable;
    @XmlAttribute
    private boolean dontAllowInitiators;
    @XmlAttribute
    private AlignmentMode mode;
    @XmlAttribute
    private int seedCount;
    @XmlAttribute
    private int maxRadius;   // set the default maximum radius
   
    @XmlElement
    public Sequence[] seedSeqs;
    @XmlElement
    public double[][] distances;
    @XmlElement
    public ScoringMatrix scoringMatrix;

    private Random rand = new Random(42);   //Biologists like replicability...apparently
    public static final float factor = 0.2f;  // allows adjustment of the delta

    public static class FramebotSearchResult {

        private String seedId;
        private int seedNum;
        private FramebotResult result;
        private int dist;
        private boolean reversed;
        private int numComparisons;

        public FramebotSearchResult(String seedId, int seedNum, FramebotResult result, int dist, boolean reversed, int numComparisons) {
            this.seedId = seedId;
            this.result = result;
            this.seedNum = seedNum;
            this.dist = dist;
            this.reversed = reversed;
            this.numComparisons = numComparisons;
        }

        public int getDist() {
            return dist;
        }

        public int getNumComparisons() {
            return numComparisons;
        }

        public FramebotResult getResult() {
            return result;
        }

        public boolean isReversed() {
            return reversed;
        }

        public String getSeedId() {
            return seedId;
        }
    }

    private FramebotIndex() {
    }
    
    public int getMaxRadius(){
        return maxRadius;
    }

    
    private static final Options options = new Options();

    static {
        options.addOption("m", "max-radius", true, "maximum radius for metric-search ONLY, range [1-" + Integer.MAX_VALUE + "]>, default is Integer.MAX_VALUE: " + Integer.MAX_VALUE);
        options.addOption("x", "scoring-matrix", true, "the metric protein scoring matrix. Default is blosum62_metric.txt from Weijia Xu's thesis: On Integrating Biological Sequence Analysis with Metric Distance");
        options.addOption("g", "gap-open-penalty", true, "gap opening penalty. Default is " + ScoringMatrix.DEFAULT_METRIC_GAP_OPEN_PEALTY);
        options.addOption("e", "gap-ext-penalty", true, "gap extension penalty. Default is " + ScoringMatrix.DEFAULT_METRIC_GAP_EXT_PENALTY);
        options.addOption("f", "frameshift-penalty", true, "frameshift penalty. Default is " + ScoringMatrix.DEFAULT_FRAME_SHIFT_PENALTY);
        options.addOption("t", "transl-table", true, "Protein translation table to use (integer based on ncbi's translation tables, default=11 bacteria/archaea)");
    }
    
    
    public static FramebotIndex index(File nuclSeeds, int maxRadius, AlignmentMode mode, boolean dontAllowInitiators, int translTable, ScoringMatrix simMatrix) throws IOException {
        List<Sequence> protSeedSeqs = new ArrayList();
        List<Sequence> nuclSeedSeqs = SequenceReader.readFully(nuclSeeds);
        FramebotIndex ret = new FramebotIndex();
        ProteinUtils protUtils = ProteinUtils.getInstance();
        
        //Seed count is the number of seed sequences, the real number however is seed count * 2 since we reverse complement all the seeds
        //to make life easier (cuts down the number of comparisons and don't have to test both the forward and reverse sequences, can do it
        //at the same time
        ret.seedCount = nuclSeedSeqs.size();
        ret.mode = mode;
        ret.dontAllowInitiators = dontAllowInitiators;
        ret.translTable = translTable;
        ret.maxRadius = maxRadius;
        ret.scoringMatrix = simMatrix;

        //First we reverse complement all the nucl seeds
        //Take note, the forward sense sequence is index - seedCount
        //We'll use this later
        for(int seed = 0;seed < ret.seedCount;seed++) {
            Sequence forward = nuclSeedSeqs.get(seed);
            nuclSeedSeqs.add(new Sequence(forward.getSeqName(), "", IUBUtilities.reverseComplement(forward.getSeqString())));
        }

        for(Sequence nuclSeed : nuclSeedSeqs) {
            protSeedSeqs.add(new Sequence(nuclSeed.getSeqName(), "", protUtils.translateToProtein(nuclSeed.getSeqString(), dontAllowInitiators, translTable)));
        }

        ret.seedSeqs = protSeedSeqs.toArray(new Sequence[protSeedSeqs.size()]);
        ret.distances = new double[ret.seedSeqs.length][ret.seedSeqs.length];
        
        for (int row = 0; row < ret.seedSeqs.length; row++) {
            Sequence nuclSeqi = nuclSeedSeqs.get(row);
            for (int col = 0; col < ret.seedSeqs.length; col++) {
                FramebotSearchResult resultij = ret.getDistance(nuclSeqi, col, 0, simMatrix);
                ret.distances[row][col] = resultij.getDist();
            }
        }

        return ret;
    }

    public static FramebotIndex readExternalIndex(File externalIndexFile, AlignmentMode mode) throws IOException {
        try {
            InputStream is = new GZIPInputStream(new BufferedInputStream(new FileInputStream(externalIndexFile)));
            FramebotIndex ret = (FramebotIndex) JAXBContext.newInstance(FramebotIndex.class).createUnmarshaller().unmarshal(is);
                 
            // Alignment mode to do FrameBot
            ret.mode = mode;
            is.close();
            
            return ret;
        } catch (JAXBException e) {
            throw new RuntimeException(e);
        }
    }

    public void writeExternalIndex(File externalIndexFile) throws IOException {
        try {
            Marshaller m = JAXBContext.newInstance(FramebotIndex.class).createMarshaller();
            m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

            OutputStream out = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(externalIndexFile)));
            m.marshal(this, out);
            out.close();
        } catch (JAXBException e) {
            throw new RuntimeException(e);
        }
    }

    private FramebotSearchResult getDistance(Sequence nuclQuery, int seedId, int numComparisons, ScoringMatrix simMatrix) {
        Sequence protSeed = seedSeqs[seedId];
        FramebotResult result = FramebotCore.processSequence(nuclQuery, protSeed, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);

        int dist = Math.abs(result.getFrameScore().getMaxScore());
        return new FramebotSearchResult(protSeed.getSeqName(), seedId, result, dist, seedId >= seedCount, numComparisons);
    }
    
    public FramebotSearchResult getNearestNeighbor(Sequence nuclQuery) {
        return getNearestNeighbor(nuclQuery,  maxRadius);
    }

    /**
     * This function will return the framebot result for the nearest seed to the query sequence
     * If the query sequence is 3' to 5' the result will be reverse complemented before returning
     *
     * @param nuclQuery
     * @return
     */
    public FramebotSearchResult getNearestNeighbor(Sequence nuclQuery, int radius) {
        List<Integer> candidates = new ArrayList(seedSeqs.length);

        for (int i = 0; i < seedSeqs.length; i++) {
            candidates.add(i);
        }
        FramebotSearchResult result = aesaNNSearchInternal(nuclQuery, radius, candidates, 0);

        //If the result is reversed...well that means we have to reverse complement it
        //This is just to furfill the general contract that framebot should return the sequence in 5' to 3' orientation
        if(result != null && result.isReversed()) {
            int forwardSeed = result.seedNum - seedCount;
            if(forwardSeed < 0 || !seedSeqs[forwardSeed].getSeqName().equals(result.getSeedId())) {
                throw new IllegalStateException("Some how the index got corrupted...gomenne");
            }

            result = getDistance(new Sequence(nuclQuery.getSeqName(), "", IUBUtilities.reverseComplement(nuclQuery.getSeqString())), forwardSeed, result.getNumComparisons(), scoringMatrix);
            result.reversed = true;
        }

        return result;
    }

    private FramebotSearchResult aesaNNSearchInternal(Sequence nuclQuery, int radius, List<Integer> candidates, int numComparisons) {
        if (candidates.isEmpty()) {
            return null;
        }

        numComparisons++;
        int seed = candidates.remove(rand.nextInt(candidates.size()));

        FramebotSearchResult result = getDistance(nuclQuery, seed, numComparisons, scoringMatrix);
        FramebotSearchResult thisHit = null;
        if (result.getDist() <= radius) {
            thisHit = result;
            radius = result.getDist();
        }

        // we optimized the result on psi but it might be possible to optimize with result.dist
        int psi =  result.dist + radius ;
        int lambda = result.dist - (int)(radius + factor*psi);
        
        List<Integer> newCandidates = new ArrayList();
        for (Integer candidate : candidates) {
            if (distances[seed][candidate] <= psi && distances[seed][candidate] >= lambda) {
                newCandidates.add(candidate);
            }
        }
       
        FramebotSearchResult possibleRet = aesaNNSearchInternal(nuclQuery, radius, newCandidates, numComparisons);
                       
        FramebotSearchResult retVal = null;
        if (possibleRet == null) {
            retVal = thisHit ;
        } else if (thisHit == null) {
            retVal = possibleRet ;
        } else if (thisHit.dist < possibleRet.dist) {
            retVal =  thisHit;
        } else {
            retVal= possibleRet;
        }
        
        return retVal;
    }

    public int getSeedCount() {
        return seedCount;
    }

    public static void main(String [] args) throws IOException {
       
        int maxRadius = Integer.MAX_VALUE;
        AlignmentMode mode = AlignmentMode.global;
        File scoringMatrixFile = null;
        int gapPenalty = ScoringMatrix.DEFAULT_METRIC_GAP_OPEN_PEALTY;
        int gapExtend = ScoringMatrix.DEFAULT_METRIC_GAP_EXT_PENALTY;
        int frameshiftPenalty = ScoringMatrix.DEFAULT_FRAME_SHIFT_PENALTY;

        File nuclSeedFile = null;
        File outFile = null;
        boolean dontAllowInitiators = true;
        int translTable = 11;
       
        
        try {
            CommandLine line = new PosixParser().parse(options, args);
                    
            if (line.hasOption("scoring-matrix")) {
                scoringMatrixFile = new File(line.getOptionValue("scoring-matrix"));
            }
            if (line.hasOption("gap-open-penalty")) {
                gapPenalty = Integer.parseInt(line.getOptionValue("gap-open-penalty"));
            }
            if (line.hasOption("gap-ext-penalty")) {
                gapExtend = Integer.parseInt(line.getOptionValue("gap-ext-penalty"));
            }
            if (line.hasOption("frameshift-penalty")) {
                frameshiftPenalty = Integer.parseInt(line.getOptionValue("frameshift-penalty"));
            }
            if (line.hasOption("transl-table")) {
                translTable = Integer.parseInt(line.getOptionValue("transl-table"));
            }
            
            if (line.hasOption("max-radius")) {
                maxRadius = Integer.parseInt(line.getOptionValue("max-radius"));
                if(maxRadius < 1) {
                    throw new Exception("Max radius must be at least 1");
                }
            }
            args = line.getArgs();
            if ( args.length != 2){
                throw new Exception("Incorrect number of command line arguments");
            }
            nuclSeedFile = new File(args[0]); 
            outFile = new File(args[1]);

        } catch (Exception e) {
            new HelpFormatter().printHelp(120, "FramebotIndex [options] <nucl seed file> <out index file>", "", options, "");
            System.out.println("ERROR: " + e.getMessage());
            return;
        }
                     
        ScoringMatrix simMatrix = null;
        if ( scoringMatrixFile!= null){
            simMatrix = new ScoringMatrix(new FileInputStream(scoringMatrixFile), gapPenalty, gapExtend, frameshiftPenalty);
        }else {
            simMatrix = new ScoringMatrix(ScoringMatrix.getDefaultProteinMatrixMetricStream(), gapPenalty, gapExtend, frameshiftPenalty);
        }
        
        FramebotIndex.index(nuclSeedFile, maxRadius, mode, dontAllowInitiators, translTable, simMatrix).writeExternalIndex(outFile);
    }
}

