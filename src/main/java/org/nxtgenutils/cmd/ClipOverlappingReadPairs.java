package org.nxtgenutils.cmd;

import java.io.File;
import java.util.Arrays;

import org.apache.log4j.Logger;
import org.nxtgenutils.impl.ReadPairClipper;
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
 * Date: 04-Oct-2011
 * Time: 10:28:19
 */


/**
 * @author Michael Mueller
 */
public class ClipOverlappingReadPairs extends AbstractCommand {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(ClipOverlappingReadPairs.class);

    public ClipOverlappingReadPairs() {
        usage = "Arguments for ClipOverlappingReadPairs\n" +
                "\n" +
                "    -i  --mappingInput <path_to_input_bam>\n" +
                "    -o  --clippedMappingOutput <path_to_output_bam>\n" +
                "    -s  --statsOutput <path_to_stats_output>\n";
    }

    public void run(String[] args) {

        File mappingInput = null;
        File mappingOutput = null;
        File statsOutput = null;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].equals("--mappingInput") || args[i].equals("-i")) {
                    String filename = args[i + 1];
                    mappingInput = new File(filename);
                }

                if (args[i].equals("--clippedMappingOutput") || args[i].equals("-o")) {
                    String filename = args[i + 1];
                    mappingOutput = new File(filename);
                }

                if (args[i].equals("--statsOutput") || args[i].equals("-s")) {
                    String filename = args[i + 1];
                    statsOutput = new File(filename);
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

        if (!mappingInput.exists()) {
            logger.error("No such file: " + mappingInput.getAbsolutePath());
            System.exit(1);
        }

        //create stats mappingOutput file
        if (mappingOutput == null) {
            mappingOutput = new File(mappingInput.getAbsolutePath().replace(".bam", ".overlapclipped.bam"));
        }

        if (statsOutput == null) {
            statsOutput = new File(mappingInput.getAbsolutePath().replace(".bam", ".overlapclipped.bam.stats"));
        }

        logger.info("Clipping overlapping read pairs, retaining higher quality end.");
        logger.info("Mapping input file : " + mappingInput.getAbsolutePath());
        logger.info("Mapping output file: " + mappingOutput);
        logger.info("Clipping stats file: " + statsOutput);

        new ReadPairClipper(mappingInput, mappingOutput, statsOutput);

    }

    public static void main(String[] args) {
        new ClipOverlappingReadPairs().run(args);
    }

}
