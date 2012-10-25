/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.framebot.dp;

import edu.msu.cme.rdp.readseq.QSequence;
import edu.msu.cme.rdp.alignment.AlignmentMode;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.framebot.core.FramebotCore;
import edu.msu.cme.rdp.framebot.core.FramebotResult;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wangqion
 */
public class FramebotCoreTest {

    public FramebotCoreTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    @Test
    public void testQual() {
        QSequence expected;
        QSequence query;
        Sequence seed;

        ScoringMatrix simMatrix = ScoringMatrix.getDefaultProteinMatrix();
        query = new QSequence(
                "GOOJNAV07IVXJZ",
                "",
                "gacgcgtctcatcctcagcgccaaggcgcaggacaccgtgctgcaccttgccgcggcgcaaggctcggtcgaggacctcgagctcgaagacgtcctcaagatcggctacaggggcatcaagtgcgtcgagtccggcggacccgagccgggcgtcggctgcgccggccgtggcgtcatcacctcaatcaacttcctcgaagagaacggcgcctacgacgacgtcgactacgtctcctacgacgtgctcggcgacgtggtgtgcggcggcttcgccatgccgatccgcgagaacaaggcgcaggaaatctacatcgtcatg",
                new byte[] { 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 34, 34, 34, 34, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 34, 31, 31, 36, 32, 32, 19, 19, 16, 16, 16, 25, 21, 30, 36, 36, 36, 32, 34, 35, 40, 40, 40, 40, 40, 40, 33, 33, 35, 40, 40, 34, 35, 35, 35, 39, 36, 36, 40, 40, 40, 40, 36, 36, 36, 36, 38, 38, 36, 33, 32, 19, 19, 18, 20, 27, 33, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 36, 36, 36, 33, 24, 26, 26, 30, 24, 26, 26, 34, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 30, 30, 40, 40, 34, 34, 33, 33, 38, 38, 36, 40, 30, 30, 26, 26, 26, 34, 40, 40, 40, 40, 40, 40, 34, 34, 30, 30, 30, 30, 40, 40, 36, 36, 20, 17, 20, 30, 30, 36, 33, 32, 32, 24, 24, 24, 30, 34, 35, 35, 40 }
                );

        seed = new Sequence("NP_435695",
                "",
                "maalrqiafygkggigksttsqntlaalvdlgqkilivgcdpkadstrlilnakaqdtvlhlaategsvedleledvlkvgyrgikcvesggpepgvgcagrgvitsinfleengayndvdyvsydvlgdvvcggfampirenkaqeiyivmsgemmalyaanniakgilkyahaggvrlgglicnerqtdreldlaealaarlnsklihfvprdnivqhaelrkmtviqyapnskqageyralaekihansgrgtvptpitmeeledmlldfgimksdeqmlaelhakeakviaph");

        expected = new QSequence(
                "GOOJNAV07IVXJZ",
                "",
                "acgcgtctcatcctcagcgccaaggcgcaggacaccgtgctgcaccttgccgcggcgcaaggctcggtcgaggacctcgagctcgaagacgtcctcaagatcggctacaggggcatcaagtgcgtcgagtccggcggacccgagccgggcgtcggctgcgccggccgtggcgtcatcacctcaatcaacttcctcgaagagaacggcgcctacgacgacgtcgactacgtctcctacgacgtgctcggcgacgtggtgtgcggcggcttcgccatgccgatccgcgagaacaaggcgcaggaaatctacatcgtcatg",
                new byte[] { 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 34, 34, 34, 34, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 34, 31, 31, 36, 32, 32, 19, 19, 16, 16, 16, 25, 21, 30, 36, 36, 36, 32, 34, 35, 40, 40, 40, 40, 40, 40, 33, 33, 35, 40, 40, 34, 35, 35, 35, 39, 36, 36, 40, 40, 40, 40, 36, 36, 36, 36, 38, 38, 36, 33, 32, 19, 19, 18, 20, 27, 33, 36, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 36, 36, 36, 33, 24, 26, 26, 30, 24, 26, 26, 34, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 30, 30, 40, 40, 34, 34, 33, 33, 38, 38, 36, 40, 30, 30, 26, 26, 26, 34, 40, 40, 40, 40, 40, 40, 34, 34, 30, 30, 30, 30, 40, 40, 36, 36, 20, 17, 20, 30, 30, 36, 33, 32, 32, 24, 24, 24, 30, 34, 35, 35, 40 });

        FramebotResult result = FramebotCore.processSequence(expected, seed, true, 11, AlignmentMode.glocal, simMatrix);

        QSequence resultSeq = (QSequence)result.getCorrectedNucl();
        assertEquals(resultSeq.getSeqName(), expected.getSeqName());
        assertEquals(resultSeq.getSeqString(), expected.getSeqString());
        assertArrayEquals(resultSeq.getQuality(), expected.getQuality());

        query = new QSequence(
                "GOOJNAV07ISLQD",
                "",
                "caccaggctcctccttggcgggctgtcccaaagaccgtgctcgataccctccgtgccgaaggggaagacctcgacctcgacgacgtgatgaagatcggttttcagggcacccgctgcgtggagtcgggagggccggagccgggtgtcggctgcgccggccggggcatcatcacttccatcaacctcctggagcagcttggcgcctactcggaaagtatcggcctcgactatgccttctatgacgtcctcggagacgtggtctgtggcggctttgccatgcccatccgtgacggaaaggccaaggagatctacatcgtggtc",
                new byte[] { 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 33, 33, 33, 34, 34, 34, 28, 19, 19, 19, 24, 22, 30, 36, 22, 27, 27, 26, 31, 32, 40, 40, 40, 40, 40, 36, 35, 30, 30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 29, 29, 29, 29, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 30, 26, 26, 26, 26, 40, 40, 39, 39, 39, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 34, 34, 34, 40, 32, 32, 32, 40, 40, 40, 40, 40, 40, 40, 40, 34, 34, 34, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 30, 29, 29, 31, 31, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 35, 34, 34, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40 }
                );

        seed = new Sequence("YP_001995946",
                "",
                "mrkvaiygkggigkstttqntvaglaemgkkvmvvgcdpkadstrlllggliqktvldtlreegedvelddiikegysatrcvesggpepgvgcagrgiitsvnlleqlgayddewnldyvfydvlgdvvcggfampirdgkaeeiyivvsgemmamyaannickgilkyadaggvrlgglicnsrkvdneqemiqelarqlgtqmihfvprdnmvqraeinrktvidydpthsqadeyrtlakkidenemfvipkpleidalekllvdfgian");

        expected = new QSequence(
                "GOOJNAV07ISLQD",
                "",
                "accaggctcctccttggcgggctgtccaagaccgtgctcgataccctccgtgccgaaggggaagacctcgacctcgacgacgtgatgaagatcggttttcagggcacccgctgcgtggagtcgggagggccggagccgggtgtcggctgcgccggccggggcatcatcacttccatcaacctcctggagcagcttggcgcctactcggaaagtatcggcctcgactatgccttctatgacgtcctcggagacgtggtctgtggcggctttgccatgcccatccgtgacggaaaggccaaggagatctacatcgtggtc",
                new byte[] { 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 33, 33, 33, 34, 34, 34, 28, 19, 19, 22, 30, 36, 22, 27, 27, 26, 31, 32, 40, 40, 40, 40, 40, 36, 35, 30, 30, 30, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 29, 29, 29, 29, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 30, 26, 26, 26, 26, 40, 40, 39, 39, 39, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 34, 34, 34, 40, 32, 32, 32, 40, 40, 40, 40, 40, 40, 40, 40, 34, 34, 34, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 30, 30, 29, 29, 31, 31, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 35, 34, 34, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40 } );

        result = FramebotCore.processSequence(expected, seed, true, 11, AlignmentMode.glocal, simMatrix);

        resultSeq = (QSequence)result.getCorrectedNucl();
        assertEquals(resultSeq.getSeqName(), expected.getSeqName());
        assertEquals(resultSeq.getSeqString(), expected.getSeqString());
        assertArrayEquals(resultSeq.getQuality(), expected.getQuality());

        query = new QSequence(
                "GOOJNAV07IKO60",
                "",
                "tacgcgactccttctcggcggtttggcgcagcgcagcgtcctggacacactccgcgaagagggcgaggacgttgaactggcggacatccgcagtggcggcttctgcaacagcctctgcgtcgagtccggcggccccgaaccggcgtaggctgcgccgggcggggcgatcatcacctccatcaatatgctcgaacatctgggcgcctacgacgagagcgaatgcctggattacgtgttttacgtacgtactgggcgacgtggtttgcggcgggttcgcgatgcctatccgcgatggcaaggtcgaggaaatctatatcgtctgc",
                new byte[] { 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 26, 26, 26, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 35, 34, 33, 36, 33, 33, 28, 28, 15, 15, 15, 15, 14, 20, 20, 12, 12, 15, 15, 21, 21, 31, 31, 32, 32, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 33, 33, 33, 33, 40, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 39, 39, 39, 39, 40, 40, 30, 30, 30, 35, 35, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 21, 21, 21, 30, 24, 34, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 36, 36, 36, 36, 36, 24, 24, 24, 24, 36, 24, 24, 14, 28, 28, 32, 32, 32, 32, 36, 33, 33, 33, 33, 21, 30, 28, 32, 33, 36, 36, 34, 34, 34, 27, 26, 26, 26, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 34, 34, 34, 32, 24, 24, 24, 25, 26, 28, 32, 35, 32, 32, 27, 27, 27, 29, 30, 36, 36, 36, 36, 36, 24, 24, 24, 36, 36, 40, 30, 30, 30, 34, 34, 40, 40, 40, 36, 36, 32 }
                );

        seed = new Sequence("ZP_03723745",
                "",
                "mrkvaiygkggigkstttqntvaglvemgkkvmvvgcdpkadstrlllgglaqrsvldtlreegedvelsdirspgfcnslcvesggpepgvgcagrgiitsinmleqlgaydeseqldyvfydvlgdvvcggfampiregkaeeiyivcsgemmamyaanniskgilkfaktgtvrlgglicnsrkvdneremiekfaeklgtkmihfvprhndvqraeinrktviewnkdceqateyrtlanaianntnfvvpnpltiqeledllmqyglln");

        expected = new QSequence(
                "GOOJNAV07IKO60",
                "",
                "acgcgactccttctcggcggtttggcgcagcgcagcgtcctggacacactccgcgaagagggcgaggacgttgaactggcggacatccgcagtggcggcttctgcaacagcctctgcgtcgagtccggcggccccgaaccgnnngtaggctgcgccgggcggggcatcatcacctccatcaatatgctcgaacatctgggcgcctacgacgagagcgaatgcctggattacgtgttttactacgtactgggcgacgtggtttgcggcgggttcgcgatgcctatccgcgatggcaaggtcgaggaaatctatatcgtctgc",
                new byte[] { 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 26, 26, 26, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 35, 34, 33, 36, 33, 33, 28, 28, 15, 15, 15, 15, 14, 20, 20, 12, 12, 15, 40, 40, 40, 21, 31, 31, 32, 32, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 40, 33, 33, 33, 33, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 39, 39, 39, 39, 40, 40, 30, 30, 30, 35, 35, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 21, 21, 21, 30, 24, 34, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 36, 36, 36, 36, 36, 24, 24, 24, 24, 36, 24, 14, 28, 28, 32, 32, 32, 32, 36, 33, 33, 33, 33, 21, 30, 28, 32, 33, 36, 36, 34, 34, 34, 27, 26, 26, 26, 40, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 36, 34, 34, 34, 32, 24, 24, 24, 25, 26, 28, 32, 35, 32, 32, 27, 27, 27, 29, 30, 36, 36, 36, 36, 36, 24, 24, 24, 36, 36, 40, 30, 30, 30, 34, 34, 40, 40, 40, 36, 36, 32 } );

        result = FramebotCore.processSequence(expected, seed, true, 11, AlignmentMode.glocal, simMatrix);

        resultSeq = (QSequence)result.getCorrectedNucl();
        assertEquals(resultSeq.getSeqName(), expected.getSeqName());
        assertEquals(resultSeq.getSeqString(), expected.getSeqString());
        assertArrayEquals(resultSeq.getQuality(), expected.getQuality());
    }


    /**
     * Test of traceback method, of class FramebotCore.
     */
    @Test
    public void testTraceback() {
        System.out.println("traceback");
        int translTable = 11;
        boolean dontAllowInitiators = true;

        Sequence s1 = new Sequence("s1", "", "GCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAG");
        Sequence s1_ins1t_d_1c = new Sequence("s1_ins1t_d_1c", "", "GCAGTTCTGTAGCGACATGTACCACGCCGGCAGATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAAATTGGATTTGTCGAATGCTCAAG");
        Sequence s1_del_t = new Sequence("s1_del_t", "", "GCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAG");

        Sequence s2 = new Sequence("s2", "", "GCAGTTCTGCTCGGACATGTACCACGCGCCGTTCAGCCATTCCTCGCCAGTGCTGGCCAGCCTGCCGCCGGATATCGACCCCAGCCAGGCCGGCTG");

        Sequence s2_2ins = new Sequence("s2_2ins", "", "GCAGTTCTGCTCGGACATGTACCACGCGCCGTTCAGCCATTCCTCGaCCAGTGCTGGCCAGCCTGCCGCCGGAgTATCGACCCCAGCCAGGCCGGCTG");
        Sequence s2_ins_beg_a = new Sequence("s2_ins_beg_a", "", "GCAGTTCTGCTCGGACATGTAaCCACGCGCCGTTCAGCCATTCCTCGCCAGTGCTGGCCAGCCTGCCGCCGGATATCGACCCCAGCCAGGCCGGCTG");
        Sequence s2_2del_mid = new Sequence("s2_2del_mid", "", "GCAGTTCTGCTCGGACATGTACCACGCGCCGTTCAGCCATTCCTCGCCAGTGCTGGCCACTGCCGCCGGATATCGACCCCAGCCAGGCCGGCTG");

        Sequence smix = new Sequence("mix", "", "AGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCAGCCTGCCGCCgGGATATCGACCCCAGCCAGGCCGGCATGGGA");
        Sequence longseq = new Sequence("longseq", "", "GCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAGTCCCTACCAAAGGGAATCAGTTCCGGGCCAAaTTGGGGCGGGCACGGAACGGGATGGTTTGTTGACGAGCCGGGCATG");
        Sequence extrafront_del_t = new Sequence("extrafront_del_t", "", "CGGCACGATGTCGCATTCCGTGCAACTGGAAGTTTGCCGCCGAGCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAG");


        Sequence subject = new Sequence("subject", "", "PCNWKFAAEQFCSDMYHAGTMSHLSGILAGMPPEMDLSNAQVPTKGNQFRANWGGHGTGW");

        ScoringMatrix simMatrix = ScoringMatrix.getDefaultProteinMatrix();
        // test global local
        FramebotResult result = FramebotCore.processSequence(s1, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);

        assertEquals("QFCSDMYHAGTMSHLSGILAGMPPEMDLSNAQ", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(0, result.getFrameshifts());
        assertEquals(162, result.getFrameScore().getMaxScore() );
        assertEquals(100, result.getPercentIdent() * 100 , 0.5);

     
        result = FramebotCore.processSequence(s1_ins1t_d_1c, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMYHAGXMSHLSGILAGMPPELDLSNAQ", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(2, result.getFrameshifts());
        assertEquals(123, result.getFrameScore().getMaxScore());

        assertEquals(94, result.getPercentIdent() * 100, 0.5);
        assertEquals("CAGTTCTGTAGCGACATGTACCACGCCGGCnnnATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAATTGGATTTGTCGAATGCTCAA", result.getCorrectedNucl().getSeqString());

      
        result = FramebotCore.processSequence(s1_del_t, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMYHAGTMSHLSGXLAGMPPEMDLSNAQ", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(97, result.getPercentIdent() * 100, 0.5);

      
        result = FramebotCore.processSequence(s2, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMYHAPFSHSSPVLASLPPDIDPSQAG", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(0, result.getFrameshifts());
        assertEquals(83, result.getFrameScore().getMaxScore());
        assertEquals(59, result.getPercentIdent() * 100, 0.5);

        result = FramebotCore.processSequence(s2_2ins, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMYHAPFSHSSPVLASLPPEIDPSQAG", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(2, result.getFrameshifts());
        assertEquals("CAGTTCTGCTCGGACATGTACCACGCGCCGTTCAGCCATTCCTCGCCAGTGCTGGCCAGCCTGCCGCCGGAgATCGACCCCAGCCAGGCCGGC", result.getCorrectedNucl().getSeqString());

        result = FramebotCore.processSequence(s2_ins_beg_a, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMNHAPFSHSSPVLASLPPDIDPSQAG", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(56, result.getPercentIdent() * 100, 0.5);
        assertEquals("CAGTTCTGCTCGGACATGAaCCACGCGCCGTTCAGCCATTCCTCGCCAGTGCTGGCCAGCCTGCCGCCGGATATCGACCCCAGCCAGGCCGGC", result.getCorrectedNucl().getSeqString());

        result = FramebotCore.processSequence(s2_2del_mid, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMYHAPFSHSSPVLALPPDIDPSQAG", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(62, result.getFrameScore().getMaxScore());
        assertEquals(59, result.getPercentIdent() * 100, 0.5);
        assertEquals("CAGTTCTGCTCGGACATGTACCACGCGCCGTTCAGCCATTCCTCGCCAGTGCTGGCCCTGCCGCCGGATATCGACCCCAGCCAGGCCGGC", result.getCorrectedNucl().getSeqString());
        assertEquals(1, result.getNumberOfOrigBasesList().get(20).intValue());
     
        result = FramebotCore.processSequence(smix, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("SDMYHAGTMSHLSGSLPPDIDPSQAGMG", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(70, result.getFrameScore().getMaxScore());
        assertEquals(61, result.getPercentIdent() * 100, 0.5);
        assertEquals("AGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCAGCCTGCCGCCgGATATCGACCCCAGCCAGGCCGGCATGGGA", result.getCorrectedNucl().getSeqString());
        assertEquals(0, result.getNumberOfOrigBasesList().get(15).intValue());
        assertEquals(4, result.getNumberOfOrigBasesList().get(21).intValue());
        //now test local
     
        result = FramebotCore.processSequence(smix, subject, dontAllowInitiators, translTable, AlignmentMode.local, simMatrix);
        assertEquals("SDMYHAGTMSHLSGSLPPDIDPSQA", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(83, result.getFrameScore().getMaxScore());
        assertEquals(68, result.getPercentIdent() * 100, 0.5);
        assertEquals("AGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCAGCCTGCCGCCgGATATCGACCCCAGCCAGGCC", result.getCorrectedNucl().getSeqString());

        // now test sequence with extra at the end
        result = FramebotCore.processSequence(longseq, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("QFCSDMYHAGTMSHLSGILAGMPPEMDLSNAQVPTKGNQFRANWGGHGTGWFVDEPGM", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(253, result.getFrameScore().getMaxScore());
        assertEquals(88, result.getPercentIdent() * 100, 0.5);
        assertEquals(58, result.getAlignedQuery().getSeqString().toUpperCase().length());
        assertEquals("CAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAGTCCCTACCAAAGGGAATCAGTTCCGGGCCAaTTGGGGCGGGCACGGAACGGGATGGTTTGTTGACGAGCCGGGCATG", result.getCorrectedNucl().getSeqString());

        result = FramebotCore.processSequence(longseq, subject, dontAllowInitiators, translTable, AlignmentMode.global, simMatrix);
        assertEquals("QFCSDMYHAGTMSHLSGILAGMPPEMDLSNAQVPTKGNQFRANWGGHGTGWFVDEPGM", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(244, result.getFrameScore().getMaxScore());
        assertEquals(76, result.getPercentIdent() * 100, 0.5);
        assertEquals(67, result.getAlignedQuery().getSeqString().toUpperCase().length());
        assertEquals("CAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAGTCCCTACCAAAGGGAATCAGTTCCGGGCCAaTTGGGGCGGGCACGGAACGGGATGGTTTGTTGACGAGCCGGGCATG", result.getCorrectedNucl().getSeqString());

        result = FramebotCore.processSequence(longseq, subject, dontAllowInitiators, translTable, AlignmentMode.local, simMatrix);
        assertEquals("QFCSDMYHAGTMSHLSGILAGMPPEMDLSNAQVPTKGNQFRANWGGHGTGW", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(100, result.getPercentIdent() * 100, 0.5);
        assertEquals(51, result.getAlignedQuery().getSeqString().toUpperCase().length());
        assertEquals("CAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCATTTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAAGTCCCTACCAAAGGGAATCAGTTCCGGGCCAaTTGGGGCGGGCACGGAACGGGATGG", result.getCorrectedNucl().getSeqString());


        // test sequences with extra in the front

        result = FramebotCore.processSequence(extrafront_del_t, subject, dontAllowInitiators, translTable, AlignmentMode.global, simMatrix);
        assertEquals("ARCRIPCNWKFAAEQFCSDMYHAGTMSHLSGXLAGMPPEMDLSNAQ", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(160,result.getFrameScore().getMaxScore());
        assertEquals(62, result.getPercentIdent() * 100, 0.5);
        assertEquals(65, result.getAlignedQuery().getSeqString().toUpperCase().length());
        assertEquals("GCACGATGTCGCATTCCGTGCAACTGGAAGTTTGCCGCCGAGCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCnnnTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAA", result.getCorrectedNucl().getSeqString());


        result = FramebotCore.processSequence(extrafront_del_t, subject, dontAllowInitiators, translTable, AlignmentMode.local, simMatrix);
        assertEquals("PCNWKFAAEQFCSDMYHAGTMSHLSGXLAGMPPEMDLSNAQ", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(204,result.getFrameScore().getMaxScore());
        assertEquals(98, result.getPercentIdent() * 100, 0.5);
        assertEquals(41, result.getAlignedQuery().getSeqString().toUpperCase().length());
        assertEquals("CCGTGCAACTGGAAGTTTGCCGCCGAGCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCnnnTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAA", result.getCorrectedNucl().getSeqString());

        result = FramebotCore.processSequence(extrafront_del_t, subject, dontAllowInitiators, translTable, AlignmentMode.glocal, simMatrix);
        assertEquals("PCNWKFAAEQFCSDMYHAGTMSHLSGXLAGMPPEMDLSNAQ", SeqUtils.getUnalignedSeqString(result.getAlignedQuery().getSeqString()).toUpperCase());
        assertEquals(1, result.getFrameshifts());
        assertEquals(189,result.getFrameScore().getMaxScore());
        assertEquals(98, result.getPercentIdent() * 100, 0.5);
        assertEquals(41, result.getAlignedQuery().getSeqString().toUpperCase().length());
        assertEquals("CCGTGCAACTGGAAGTTTGCCGCCGAGCAGTTCTGTAGCGACATGTACCACGCCGGCACGATGTCGCACCTGTCCGGCnnnTTGGCGGGCATGCCGCCAGAAATGGATTTGTCGAATGCTCAA", result.getCorrectedNucl().getSeqString());

        
    }

}