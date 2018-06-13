// Modified Tassel4 pipeline
// Jesse Poland, Kansas State University
// Narinder Singh, Kansas State University

package pipeline;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import net.maizegenetics.gbs.pipeline.FastqToTBTPlugin;
import net.maizegenetics.gbs.pipeline.FastqToTagCountPlugin;
import net.maizegenetics.gbs.pipeline.MergeMultipleTagCountPlugin;
import net.maizegenetics.gbs.pipeline.MergeTagsByTaxaFilesPlugin;
import net.maizegenetics.gbs.pipeline.QseqToTBTPlugin;
import net.maizegenetics.gbs.pipeline.QseqToTagCountPlugin;
import net.maizegenetics.gbs.pipeline.TagsToSNPsNoAnchor;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBit;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;

public class Tassel4Pipeline {

	private static String projectName = "_170419_tauschii_wgrc"; //project name
	private static String userName = "user"; //user name

	private static String keyFile = "/bulk/"+userName+"/gbs/jobs/" + projectName + ".txt";
	private static String seqDir = "/bulk/user/sequence";
	
	// tagcount defines tags to keep. keyfile tells us which tags to keep
	private static String tagcountsDir = "/bulk/"+userName+"/gbs/projects/" + projectName + "/tagcounts"; //counts of how much each 64bp unique sequence shows up
	private static String masterTagsFile = "/bulk/"+userName+"/gbs/projects/" + projectName + "/MasterTags_" + projectName + ".cnt"; //filtering based on limitations
	private static String tbtDir = "/bulk/"+userName+"/gbs/projects/" + projectName + "/tbt";
	private static String tbtFile = "/bulk/"+userName+"/gbs/projects/" + projectName + "/tbt_" + projectName + "_" + ".bin"; //combine library data files
	private static String tbtFileMerge = "/bulk/"+userName+"/gbs/projects/" + projectName + "/tbtMerge" + "_" + projectName + "_" + ".bin"; // tags by taxa
	private static String hapDir = "/bulk/"+userName+"/gbs/projects/" + projectName + "/" + projectName + "/"; // Output for hapmaps of SNP calls

	private static String enzyme = "PstI-MspI"; // looks for barcode and restriction site

        public static void main(String[] args) throws IOException {

                long start = System.currentTimeMillis();

                new File("/bulk/"+userName+"/gbs/projects/" + projectName).mkdirs();
                new File("/bulk/"+userName+"/gbs/projects/" + projectName + "/tbt").mkdirs();
                new File("/bulk/"+userName+"/gbs/projects/" + projectName + "/tagcounts").mkdirs();
                new File("/bulk/"+userName+"/gbs/projects/" + projectName + "/" + projectName).mkdirs();
                new File("/bulk/"+userName+"/gbs/projects/" + projectName + "/keyPipeAndMisc").mkdirs();
            
                runFastqToTagCountPlugin(); // count tags
                runMergeMultipleTagCountPlugin(); // merge tags

                runFastqToTBTPlugin(); // tags to tbt
                runMergeTagsByTaxaFilesPlugin(); // merge tbt files
                mergeTaxaInTBT(); // combine duplicated sample data

                runTagsToSNPsNoAnchor();

                long time = (System.currentTimeMillis() - start) / 1000 / 60;
                System.out.println();
                System.out.println("Processed in : " + time + " min");
                System.out.println("Done!");
        }

        private static String getDate() {
                DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
                Date date = new Date();
                return dateFormat.format(date);
        }

        public static void runQseqToTagCountPlugin() {

                String[] args = new String[] { 
                                "-i", seqDir, 
                                "-k", keyFile, 
                                "-e", enzyme, // Enzyme used to create the GBS library
                                "-s", "250000000", // Max good reads per lane. (Optional.
                                                                        // Default is 200,000,000)
                                "-c", "1", // Minimum tag count (default is 1)
                                "-o", tagcountsDir };

                QseqToTagCountPlugin plugin = new QseqToTagCountPlugin();
                plugin.setParameters(args);
                plugin.performFunction(null);
        }

        public static void runFastqToTagCountPlugin() {
                String[] testArgs = new String[] { 
                                "-i", seqDir, 
                                "-k", keyFile, 
                                "-e", enzyme, // Enzyme used to create the GBS library
                                "-s", "250000000", // Max good reads per lane. (Optional.
                                                                        // Default is 200,000,000)
                                "-c", "1", // Minimum tag count (default is 1)
                                "-o", tagcountsDir, };

                String[] args = testArgs;
                FastqToTagCountPlugin plugin = new FastqToTagCountPlugin();
                plugin.setParameters(args);
                plugin.performFunction(null);
        }

        public static void runMergeMultipleTagCountPlugin() {

                String[] args = new String[] { 
                                "-i", tagcountsDir, // Input directory
                                                                        // containing .cnt
                                                                        // files
                                "-o", masterTagsFile, // Output file name
                                "-c", "5", // Minimum count of reads to be output (default 1)
                };

                MergeMultipleTagCountPlugin plugin = new MergeMultipleTagCountPlugin();
                plugin.setParameters(args);
                plugin.performFunction(null);
        }

        public static void runQseqToTBTPlugin() {
                String[] testArgs = new String[] { 
                                "-i", seqDir, 
                                "-k", keyFile, 
                                "-e", enzyme, // Enzyme used to create the GBS library
                                "-o", tbtDir, 
                                "-c", "1", // Minimum tag count (default is 1).
                                "-t", masterTagsFile, // master Tags file

                };

                String[] args = testArgs;
                QseqToTBTPlugin plugin = new QseqToTBTPlugin();
                plugin.setParameters(args);
                plugin.performFunction(null);
        }

        public static void runFastqToTBTPlugin() {

                String[] testFastqArgs = new String[] { 
                                "-i", seqDir, 
                                "-k", keyFile,
                                "-e", enzyme, // Enzyme used to create the GBS library
                                "-o", tbtDir, 
                                "-c", "1", // Minimum tag count (default is 1).
                                "-t", masterTagsFile, // master Tags file
                                "-y", };

                String[] args = testFastqArgs;
                FastqToTBTPlugin plugin = new FastqToTBTPlugin();
                plugin.setParameters(args);
                plugin.performFunction(null);
        }

        public static void runMergeTagsByTaxaFilesPlugin() {
                String[] testArgs = new String[] { 
                                "-i", tbtDir, 
                                "-o", tbtFile, };

                String[] args = testArgs;
                MergeTagsByTaxaFilesPlugin plugin = new MergeTagsByTaxaFilesPlugin();
                plugin.setParameters(args);
                plugin.performFunction(null);
        }

        public static void mergeTaxaInTBT() {
                String inputTBTFileS = tbtFile;
                String outputMergedTaxaTBTFileS = tbtFileMerge;
                TagsByTaxaUtils.mergeTaxaByName(inputTBTFileS,
                                outputMergedTaxaTBTFileS, FilePacking.Bit, true);
                TagsByTaxaUtils.streamBinaryToText(outputMergedTaxaTBTFileS, 10000);
        }

        public static void printSumCountsInTBTByTaxa() {
                String inputTBTFileS = tbtFileMerge;
                TagsByTaxaUtils.printSumCounts(inputTBTFileS, FilePacking.Bit, true);
        }

        public static void runTagsToSNPsNoAnchor() {

                System.out.println("Starting TagsToSNPsNoAnchor");

                TagsByTaxa theTBT = new TagsByTaxaBitFileMap(tbtFileMerge);

                String outHapMap = hapDir + projectName + "_1" + ".hap";

                int nDiff = 1;
                double maxMissingData = 0.8;
                double minorAlleleFreq = 0.01;
                double maxHet = 0.1;
                boolean isDHpopulation = false;
                boolean isBiparental = false;
                boolean callHets = true;
                double pVal = 0.001;

                new TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff, maxMissingData,
                                minorAlleleFreq, maxHet, isDHpopulation, isBiparental,
                                callHets, pVal);
                System.gc();

                System.out.println("\n" + "DONE!");

                System.out.println("Starting TagsToSNPsNoAnchor");

                outHapMap = hapDir + projectName + "_2" + ".hap";
                
                nDiff = 2;
                new TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff, maxMissingData,
                                minorAlleleFreq, maxHet, isDHpopulation, isBiparental,
                                callHets, pVal);
                System.gc();

                System.out.println("\n" + "DONE!");

                System.out.println("Starting TagsToSNPsNoAnchor");

                outHapMap = hapDir + projectName + "_3" + ".hap";
                
                nDiff = 3;
                new TagsToSNPsNoAnchor(theTBT, outHapMap, nDiff, maxMissingData,
                                minorAlleleFreq, maxHet, isDHpopulation, isBiparental,
                                callHets, pVal);
                System.gc();

                System.out.println("\n" + "DONE!");

        }

        public static void runMapTagsWithGenetic() throws IOException {

                System.out.println("Starting MapTagsWithGenetic");
                String binMapUpdate = "";
                TagsByTaxa theTBT = new TagsByTaxaBitFileMap(tbtFileMerge);
                double sigThreshold = 0.001;
                int minCountForTesting = 10;
                int initSkimRate = 1;
                String outFile = "";

                System.out.println("\n" + "DONE!");

        }

}
