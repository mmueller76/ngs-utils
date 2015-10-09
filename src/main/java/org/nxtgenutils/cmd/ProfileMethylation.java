package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.nxtgenutils.NxtGenUtils;
import org.nxtgenutils.bsseq.impl.MethylationProfiler;

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
 * Date: 26-Oct-2011
 * Time: 13:08:31
 */

/**
 * @author Michael Mueller
 */
public class ProfileMethylation extends AbstractCommand {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(ProfileMethylation.class);

    public ProfileMethylation() {
        usage = "Arguments for ProfileMethylation\n" +
                "\n" +
                "    -i  --pileupInput <path_to_mpileup>\n" +
                "    -o  --profileOutput <path_to_output>\n" +
                "   [-r  --estimateBisulfiteConversionRate [<control_sample_sequence_name>]]\n" +
                "   [-c  --controlSamplePileup <path_to_controle_sample_mpileup>]\n" +
                "   [-s  --sampleNames <sample_name_1,sample_name_2,...,sample_name_n>]";
    }

    public void run(String[] args) {

        File pileupInput = null;
        File profileOutput = null;
        File controleSamplePileupInput = null;
        String estimateBisulfiteConversionRate = null;
        List<String> sampleNames = null;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].equals("--pileupInput") || args[i].equals("-i")) {
                    String filename = args[i + 1];
                    pileupInput = new File(filename);
                }

                if (args[i].equals("--profileOutput") || args[i].equals("-o")) {
                    String filename = args[i + 1];
                    profileOutput = new File(filename);
                }

                if (args[i].equals("--controleSamplePileup") || args[i].equals("-c")) {
                    String filename = args[i + 1];
                    controleSamplePileupInput = new File(filename);
                }

                if (args[i].equals("--estimateBisulfiteConversionRate") || args[i].equals("-r")) {

                    try {
                        if (i + 1 < args.length) {
                            estimateBisulfiteConversionRate = args[i + 1];
                        } else {
                            estimateBisulfiteConversionRate = "";
                        }
                    } catch (ArrayIndexOutOfBoundsException e) {
                        logger.error("Exception while parsing --estimateBisulfiteConversionRate command line argument.", e);
                    }
                    if (estimateBisulfiteConversionRate.startsWith("-")) {
                        estimateBisulfiteConversionRate = "";
                    }
                }

                if (args[i].equals("--sampleNames") || args[i].equals("-s")) {
                    String[] names = args[i + 1].split(",");
                    sampleNames = Arrays.asList(names);
                }

            }

        } catch (Exception e) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            System.out.println(Arrays.toString(args));
            logger.error("Exception while reading command line arguments: " + e.getMessage());
            System.exit(1);
        }

        if (pileupInput == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --pileupInput missing.");
            System.exit(1);
        }

        if (profileOutput == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --profileOutput missing.");
            System.exit(1);
        }

        if (!pileupInput.exists()) {
            logger.error("No such file: " + pileupInput.getAbsolutePath());
            System.exit(1);
        }

        logger.info("Profiling methylation from BS-seq pileup");
        logger.info("Sample pileup input file : " + pileupInput.getAbsolutePath());
        if (controleSamplePileupInput != null) {
            logger.info("Controle pileup input file : " + controleSamplePileupInput.getAbsolutePath());
        }
        logger.info("Writing profile to: " + profileOutput.getAbsolutePath() + ".gz");

        if (estimateBisulfiteConversionRate != null) {
            if (estimateBisulfiteConversionRate.equals("")) {
                logger.info("Estimating bisulfite conversion rate from non-CpG conversion");
            } else {
                logger.info("Estimating bisulfite conversion rate from unmethylated control sample: " + estimateBisulfiteConversionRate);
            }
        } else {
            logger.info("Bisulfite conversion rate will not be estimated.");
        }

        if (sampleNames != null) {
            logger.info("Sample names: " + sampleNames.toString().replace("[", "").replace("]", ""));
        } else {
            logger.info("Sample names not provided. Will be auto-generated");
        }

        logger.info("--------------------------------------------------------");

        new MethylationProfiler(pileupInput, controleSamplePileupInput, profileOutput, estimateBisulfiteConversionRate, sampleNames);

    }

    public static void main(String[] args) {
        new ProfileMethylation().run(args);
    }


}
