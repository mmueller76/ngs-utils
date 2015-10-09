package org.nxtgenutils.cmd;

import org.apache.log4j.Logger;

import java.io.File;
import java.util.Arrays;

import org.nxtgenutils.ReadPhaser;
import org.nxtgenutils.impl.SingleEndReadPhaser;
import org.nxtgenutils.impl.PairedEndReadPhaser;
import org.nxtgenutils.impl.BSPairedEndReadPhaser;
import org.nxtgenutils.impl.BSSingleEndReadPhaser;
import org.nxtgenutils.NxtGenUtils;
import org.nxtgenutils.Strand;

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
 * Date: 10-Oct-2011
 * Time: 12:00:34
 */

/**
 * @author Michael Mueller
 */
public class PhaseMapping extends AbstractCommand {

    private static Logger logger = Logger.getLogger(PhaseMapping.class);

    public PhaseMapping() {
        usage = "Arguments for PhaseMapping" +
                "\n" +
                "    -i --mappingInput <path_to_input_bam>\n" +
                "    -g --genotypeVCF <path_to_genotype_vcf>\n" +
                "    -a --sampleA <name_sample_a>\n" +
                "    -b --sampleB <name_sample_b>\n" +
                "[   -p --phaseMate <TRUE|FALSE> (default TRUE)]\n" +
                "[   -m --minorAlleleFrequencyCutoff <rel. frequency>]\n" +
                "[   -c --isBisulphiteConverted <TRUE|FALSE (default FALSE)>]\n" +
                "[   -s --bisulphiteConvertedStrand <FORWARD|REVERSE (default FORWARD)>]";
    }

    public void run(String[] args) {

        File mappingInput = null;
        File genotypeVCF = null;
        String sampleA = null;
        String sampleB = null;
        double minorAlleleFrequencyCutoff = -1;
        boolean phaseMate = true;
        boolean isBisulphiteConverted = false;
        Strand bisulphiteConvertedStrand = Strand.FORWARD;

        try {

            for (int i = 0; i < args.length; i++) {

                if (args[i].startsWith("--mappingInput") || args[i].startsWith("-i")) {
                    String filename = args[i + 1];
                    mappingInput = new File(filename);
                }

                if (args[i].startsWith("--genotypeVCF") || args[i].startsWith("-g")) {
                    String filename = args[i + 1];
                    genotypeVCF = new File(filename);
                }

                if (args[i].startsWith("--sampleA") || args[i].startsWith("-a")) {
                    sampleA = args[i + 1];
                }

                if (args[i].startsWith("--sampleB") || args[i].startsWith("-b")) {
                    sampleB = args[i + 1];
                }

                if (args[i].startsWith("--minorAlleleFrequencyCutoff") || args[i].startsWith("-m")) {
                    minorAlleleFrequencyCutoff = Double.parseDouble(args[i + 1]);
                }

                if (args[i].startsWith("--phaseMate") || args[i].startsWith("-p")) {
                    phaseMate = Boolean.parseBoolean(args[i + 1]);
                }

                if (args[i].startsWith("--isBisulphiteConverted") || args[i].startsWith("-c")) {
                    isBisulphiteConverted = Boolean.parseBoolean(args[i + 1]);
                }

                if (args[i].startsWith("--bisulphiteConvertedStrand") || args[i].startsWith("-s")) {
                    bisulphiteConvertedStrand = Strand.valueOf(args[i + 1]);
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

        if (genotypeVCF == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --genotypeVCF missing.");
            System.exit(1);
        }

        if (sampleA == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --sampleA missing.");
            System.exit(1);
        }

        if (sampleB == null) {
            System.out.println(NxtGenUtils.usage);
            System.out.println(usage);
            logger.error("Required input argument --sampleB missing.");
            System.exit(1);
        }

        logger.info("Phasing read mapping...");
        logger.info("Input BAM file   : " + mappingInput.getAbsolutePath());
        logger.info("Genotype VCF file: " + genotypeVCF.getAbsolutePath());
        logger.info("sample A         : " + sampleA);
        logger.info("sample B         : " + sampleB);
        logger.info("phase mate       : " + phaseMate);

        if (minorAlleleFrequencyCutoff != -1) {
            logger.info("minum minor allele frequency heterozygous SNPs: " + minorAlleleFrequencyCutoff);
        }

        if (isBisulphiteConverted) {
            logger.info("bisulphite converted DNA: " + isBisulphiteConverted);
        }

        if (isBisulphiteConverted) {
            logger.info("genomic template strand: " + bisulphiteConvertedStrand);
        }

        ReadPhaser newPhaser;
        if (isBisulphiteConverted) {
            if (phaseMate) {
                newPhaser = new BSPairedEndReadPhaser(mappingInput, genotypeVCF, bisulphiteConvertedStrand);
            } else {
                newPhaser = new BSSingleEndReadPhaser(mappingInput, genotypeVCF, bisulphiteConvertedStrand);
            }
        } else {
            if (phaseMate) {
                newPhaser = new PairedEndReadPhaser(mappingInput, genotypeVCF);
            } else {
                newPhaser = new SingleEndReadPhaser(mappingInput, genotypeVCF);
            }
        }

        newPhaser.setHeterozygousMinorAlleleFrequencyCutoff(minorAlleleFrequencyCutoff);
        newPhaser.phase(sampleA, sampleB);

    }

    public static void main(String[] args) {
        new PhaseMapping().run(args);
    }

}
