package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;

import java.io.File;
import java.util.Arrays;

import org.nxtgenutils.bsseq.impl.PairedEndBSMappingProcessor;
import org.nxtgenutils.bsseq.impl.SingleEndBSMappingProcessor;
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
 * Date: 18-Feb-2011
 * Time: 09:47:52
 */

/**
 * @author Michael Mueller
 */
public class PostProcessBSMapping extends AbstractCommand {


    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(PostProcessBSMapping.class);

    public PostProcessBSMapping() {
        usage = "Arguments for PostProcessBSMapping\n" +
                "\n" +
                "    -f --forwardMapping <path_to_forward_strand_mapping_bam>\n" +
                "    -r --reverseMapping <path_to_reverse_strand_mapping_bam>\n" +
                "    -s --referenceSequence <path_to_reference_sequence_fasta>\n" +
                "   [-o --outputDirectory <path_to_output_directory>]" +
                "   [-p --pairedEnd <TRUE|FALSE> (default TRUE)]";
    }

    public void run(String[] args) {

        File mappingForwardStrand = null;
        File mappingReverseStrand = null;
        File referenceSequence = null;
        File outputDirectory = null;
        boolean isPairedEndMapping = true;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].startsWith("-f") || args[i].startsWith("--forwardMapping")) {
                    String filename = args[i + 1];
                    mappingForwardStrand = new File(filename);
                }

                if (args[i].startsWith("-r") || args[i].startsWith("--reverseMapping")) {
                    String filename = args[i + 1];
                    mappingReverseStrand = new File(filename);
                }

                if (args[i].startsWith("-s") || args[i].startsWith("--referenceSequence")) {
                    String filename = args[i + 1];
                    referenceSequence = new File(filename);
                }

                if (args[i].startsWith("-o") || args[i].startsWith("--outputDirectory")) {
                    String filename = args[i + 1];
                    outputDirectory = new File(filename);
                }

                if (args[i].startsWith("-p") || args[i].startsWith("--pairedEnd")) {
                    String value = args[i + 1];
                    isPairedEndMapping = Boolean.getBoolean(value);
                }


            }

        } catch (Exception e) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            System.out.println(Arrays.toString(args));
            logger.error("Exception while reading command line arguments: " + e.getMessage());
            System.exit(1);
        }


        if (mappingForwardStrand == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --forwardMapping missing.");
            System.exit(1);
        }

        if (mappingReverseStrand == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --reverseMapping missing.");
            System.exit(1);
        }

        if (referenceSequence == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --referenceSequence missing.");
            System.exit(1);
        }


        if (isPairedEndMapping) {
            new PairedEndBSMappingProcessor(mappingForwardStrand, mappingReverseStrand, referenceSequence, outputDirectory);
        } else {
            new SingleEndBSMappingProcessor(mappingForwardStrand, mappingReverseStrand, referenceSequence, outputDirectory);
        }

    }

    public static void main(String[] args) {
        new PostProcessBSMapping().run(args);
    }

}
