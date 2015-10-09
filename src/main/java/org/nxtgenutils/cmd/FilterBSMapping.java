package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;

import java.io.File;
import java.util.Arrays;

import org.nxtgenutils.bsseq.BSMappingFiltrator;
import org.nxtgenutils.NxtGenUtils;

/**
 * This file is part of NxtGenUtils.
 * <p/>
 * NxtGenUtils is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p/>
 * NxtGenUtils is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p/>
 * You should have received a copy of the GNU General Public License
 * along with NxtGenUtils.  If not, see <http://www.gnu.org/licenses/>.
 * <p/>
 * BT:A:<value> BS converted template strand the mapping record
 * refers to (character attribute)
 * BT:A:F BS-converted forward strand
 * BT:A:R BS-converted reverse strand
 * <p/>
 * BS:A:<value> strand specificity of mapping (character attribute)
 * BS:A:F read maps to BS-converted forward strand only
 * BS:A:R read maps to BS-converted reverse strand only
 * BS:A:B read maps to both BS-converted strands
 * BS:A:N read does not map to either strand
 * <p/>
 * BP:A:<value> strand with proper read pairing (character attribute)
 * BP:A:F read properly paired on BS-converted forward strand only
 * BP:A:R read properly paired on BS-converted reverse strand only
 * BP:A:B read properly paired on both BS-converted strands
 * BP:A:N read not properly paired on either strand
 * <p/>
 * BW:i:<value>
 * Phred scaled probablitity of read mapped to incorrect strand
 * based on C residue frequency (integer attribute)
 * <p/>
 * BV:i:<frequency>
 * number of valid BS-mismatches to the reference, i.e. T in read and C in reference
 * (forward strand: unconverted reference = C && unconverted read = T,
 * reverse strand: unconverted reference = G && unconverted read = A)
 * <p/>
 * BI:i:<frequency>
 * number of BS-invalid mismatches to the reference, i.e. C in read and T in reference
 * i.e. due to sequencing errors/SNPs
 * <p/>
 * BN:i:<frequency>
 * number of non-BS mismatches to the reference, i.e. any mismatch other than T/C mismatches on forward strand and G/A mismatches on reverse strand
 * because of sequencing errors/SNPs
 * <p/>
 * BC:Z:<frequency>
 * ratio(number of unconverted Cs/total number of Cs) outside CpG context
 * <p/>
 * BM:Z:<frequency>
 * ratio(number of unconverted Cs/total number of Cs) in CpG context
 * BO:Z:<value> BS-mapping orientation
 * BO:Z:V valid:
 * - forward strand template: read1 mapping to (+) strand and read2 mapping to (-) strand
 * - reverse strand template: read1 mapping to (-) strand and read2 mapping to (+) strand
 * BO:Z:I invalid
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 03-Feb-2011
 * Time: 13:05:33
 */

/**
 * @author Michael Mueller
 */
public class FilterBSMapping extends AbstractCommand {


    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(FilterBSMapping.class);

    public FilterBSMapping() {
        usage = "Arguments for FilterBSMapping\n" +
                "\n" +
                "    -i --mappingInput <path_to_input_bam>\n" +
                "   [-o --filteredMappingOutput <path_to_filtered_output_mapping_bam>]\n" +
                "    -f --filterSet <path_to_filter_set_file>";
    }

    public void run(String[] args) {

        File mappingInput = null;
        File mappingOutput = null;
        File filterSet = null;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].startsWith("--mappingInput") || args[i].equals("-i")) {
                    String filename = args[i + 1];
                    mappingInput = new File(filename);
                }

                if (args[i].startsWith("--filteredMappingOutput") || args[i].equals("-o")) {
                    String filename = args[i + 1];
                    mappingOutput = new File(filename);
                }

                if (args[i].startsWith("--filterSet") || args[i].equals("-f")) {
                    String filename = args[i + 1];
                    filterSet = new File(filename);
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

        if (filterSet == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --filterSet missing.");
            System.exit(1);
        }

        new BSMappingFiltrator(mappingInput, mappingOutput, filterSet);

    }

    public static void main(String args[]) {
        new CountReads().run(args);
    }


}
