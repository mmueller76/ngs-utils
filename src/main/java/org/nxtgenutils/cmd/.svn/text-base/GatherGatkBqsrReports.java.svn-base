package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Gatherer;
import org.broadinstitute.sting.gatk.walkers.bqsr.BQSRGatherer;
import org.nxtgenutils.NxtGenUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


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
 * Created with IntelliJ IDEA.
 * User: mmuelle1
 * Date: 3/2/13
 * Time: 3:47 PM
 */


/**
 * @author Michael Mueller
 */
public class GatherGatkBqsrReports extends AbstractCommand {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(GatherGatkBqsrReports.class);

    public GatherGatkBqsrReports() {
        usage = "Arguments for GatherGatkBqsrReports\n" +
                "\n" +
                "    -i  --reportInput <path_to_basecall_recalibration_report>\n" +
                "          (argument can be specified multiple times for multiple input files)\n" +
                "    -o  --mergedReportOutput <path_to_merged_report_output>\n";
    }


    @Override
    public void run(String[] args) {

        List<File> inputReports = new ArrayList<File>();
        File outputReport = null;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].equals("--reportInput") || args[i].equals("-i")) {
                    String filename = args[i + 1];
                    inputReports.add(new File(filename));
                }

                if (args[i].equals("--mergedReportOutput") || args[i].equals("-o")) {
                    String filename = args[i + 1];
                    outputReport = new File(filename);
                }

            }

        } catch (Exception e) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            System.out.println(Arrays.toString(args));
            logger.error("Exception while reading command line arguments: " + e.getMessage());
            System.exit(1);
        }

        if (inputReports.size() == 0) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --reportInput missing.");
            System.exit(1);
        }

        if (outputReport == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --mergedReportOutput missing.");
            System.exit(1);
        }

        logger.info("Merging GATK reclibration reports...");
        logger.info("Report input file(s) : " + inputReports.toString());
        logger.info("Report output file: " + outputReport);

        new BQSRGatherer().gather(inputReports, outputReport);

        logger.info("done");

    }

    public static void main(String[] args) {
        new GatherGatkBqsrReports().run(args);
    }

}
