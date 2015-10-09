package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;
import org.nxtgenutils.NxtGenUtils;
import org.nxtgenutils.impl.PrimerClipper;
import org.nxtgenutils.impl.ReadPairClipper;

import java.io.File;
import java.util.Arrays;

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
 * Date: 09-Aug-2013
 * Time: 10:28:19
 */

/**
 * @author Michael Mueller
 */
public class ClipPrimerSequences extends AbstractCommand {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(ClipPrimerSequences.class);

    public ClipPrimerSequences() {
        usage = "Arguments for ClipPrimerSequences\n" +
                "\n" +
                "    -i  --mappingInput <path_to_input_bam>\n" +
                "   [-o  --clippedMappingOutput <path_to_output_bam>]\n" +
                "    -p  --primerCoordinates <path_to_primer_genomic_coordinates_bed>\n" +
                "   [-s  --primerCoordinatesOffSet <offset_bases>; default = 10\n";
    }

    public void run(String[] args) {

        File mappingInput = null;
        File mappingOutput = null;
        File primerCoordinates = null;
        int primerCoordinatesOffSet = 10;

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

                if (args[i].equals("--primerCoordinates") || args[i].equals("-p")) {
                    String filename = args[i + 1];
                    primerCoordinates = new File(filename);
                }

                if (args[i].equals("--primerCoordinatesOffSet") || args[i].equals("-s")) {
                    primerCoordinatesOffSet = Integer.parseInt(args[i + 1]);
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
            mappingOutput = new File(mappingInput.getAbsolutePath().replace(".bam", ".clipped.bam"));
        }

        if (primerCoordinates == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --primerCoordinates missing.");
            System.exit(1);
        }

        logger.info("Clipping primer sequences from reads.");
        logger.info("Mapping input file : " + mappingInput.getAbsolutePath());
        logger.info("Mapping output file: " + mappingOutput);
        logger.info("Primer genomic coordinates: " + primerCoordinates);
        logger.info("Primer coordinates offset: +/-" + primerCoordinatesOffSet + " bases");

        new PrimerClipper(mappingInput, mappingOutput, primerCoordinates, primerCoordinatesOffSet);

    }

    public static void main(String[] args) {
        new ClipOverlappingReadPairs().run(args);
    }

}
