package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;

import java.io.File;
import java.util.Arrays;

import org.nxtgenutils.impl.ReadCounter;
import org.nxtgenutils.NxtGenUtils;

/**
 * This file is part of NxtGenUtils.
 *
 * NxtGenUtils is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * NxtGenUtils is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with NxtGenUtils.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 14-Oct-2011
 * Time: 15:29:26
 */

/**
 * @author Michael Mueller
 */
public class CountReads extends AbstractCommand {

    private static Logger logger = Logger.getLogger(CountReads.class);

    public CountReads() {
        usage = "Arguments for CountReads\n" +
                "\n" +
                "    -i --mappingInput <path_to_bam_file>\n" +
                "    -r --regionInput <path_to_region_bed>\n" +
                "    -o --coverageOutput <path_to_output_file>\n" +
                "   [-q --minMappingQuality <integer> (default = 0)]\n" +
                "   [-m --mergeReadPairCounts <true|false> (default = false)]\n";
    }

    public void run(String[] args) {

        File mappingInput = null;
        File mappingInputIndex = null;
        File regionInput = null;
        File coverageOutput = null;
        int minMappingQuality = 0;
        boolean mergeOverlappingReadPairs = false;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].equals("--mappingInput") || args[i].equals("-i")) {
                    String filename = args[i + 1];
                    mappingInput = new File(filename);
                }

                if (args[i].startsWith("--regionInput") || args[i].equals("-r")) {
                    String filename = args[i + 1];
                    regionInput = new File(filename);
                }

                if (args[i].startsWith("--coverageOutput") || args[i].equals("-o")) {
                    String filename = args[i + 1];
                    coverageOutput = new File(filename);
                }

                if (args[i].startsWith("--minMappingQuality") || args[i].equals("-q")) {
                    String quality = args[i + 1];
                    minMappingQuality = Integer.parseInt(quality);
                }

                if (args[i].startsWith("--mergeReadPairCounts") || args[i].equals("-m")) {
                    String merge = args[i + 1];
                    mergeOverlappingReadPairs = Boolean.parseBoolean(merge);
                }

            }


        } catch (Exception e) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            System.out.println(Arrays.toString(args));
            logger.error("Exception while reading command line arguments: " + e.getMessage());
            System.exit(1);
        }

        if (mappingInput == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --mappingInput missing.");
            System.exit(1);
        }

        String bamIndexFilePath = mappingInput.getAbsolutePath() + ".bai";
        if (!new File(bamIndexFilePath).exists()) {
            logger.error("Unable to locate index file: " + mappingInput.getAbsolutePath());
            logger.error("CountReads works only on BAM files indexed with Picard Tools." + mappingInput.getAbsolutePath());
            logger.error("Index BAM file first with Picard Tools BuildBamIndex." + mappingInput.getAbsolutePath());
            System.exit(1);
        }

        if (!mappingInput.exists()) {
            logger.error("No such file: " + mappingInput.getAbsolutePath());
            System.exit(1);
        }

        if (regionInput == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --regionInput missing.");
            System.exit(1);
        }

        if (coverageOutput == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --coverageOutput missing.");
            System.exit(1);
        }

        logger.info("Counting reads by region");
        logger.info("Coverage input file                              : " + mappingInput.getAbsolutePath());
        logger.info("Region input file                                : " + regionInput.getAbsolutePath());
        logger.info("Coverage output file                             : " + coverageOutput.getAbsolutePath());
        logger.info("Min. mapping quality                             : " + minMappingQuality);
        logger.info("Merging counts of read pairs mapping to same exon: " + mergeOverlappingReadPairs);

        new ReadCounter(mappingInput, regionInput, coverageOutput, minMappingQuality, mergeOverlappingReadPairs);

    }

    public static void main(String args[]) {
        new CountReads().run(args);
    }

}
