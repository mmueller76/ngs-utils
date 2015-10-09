package org.nxtgenutils.playground;

import org.nxtgenutils.io.VCFParser;
import org.nxtgenutils.io.VCFRecord;
import org.nxtgenutils.io.VCFSampleRecord;
import org.nxtgenutils.io.impl.SimpleVCFParser;
import org.apache.log4j.Logger;
import net.sf.samtools.SAMRecord;

import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;

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
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 25-Jan-2012
 * Time: 15:06:45
 */
public class SnpExtractor {

    private static Logger logger = Logger.getLogger(SnpExtractor.class);


    public static void main(String[] args) {

        String sampleA = "BN.Lx";
        String sampleB = "SHR";
        double minorAlleleFrequencyCutoff = 0.3;
        String genotypeVCF = "/home/mmuelle1/dev/java/nxtgen-utils/testdata/22genome_snp_het06_mmq30_mbq20_soft_filter_t30_ts99_g8.vcf";
        String genotypeBed = "/home/mmuelle1/dev/java/nxtgen-utils/testdata/22genome_snp_het06_mmq30_mbq20_soft_filter_t30_ts99_g8." + sampleA + "." + sampleB + ".bed";

        VCFParser parser = new SimpleVCFParser(new File(genotypeVCF));
        parser.setHeterozygousMinorAlleleFrequencyCutoff(minorAlleleFrequencyCutoff);

        //iterate over SNPs in VCF file
        Iterator<VCFRecord> vcfIterator = parser.iterator();

        int snpCount = 0;
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(genotypeBed);
        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        while (vcfIterator.hasNext()) {

            VCFRecord vcfRecord = vcfIterator.next();

            snpCount++;
            if (snpCount % 100000 == 0) {
                logger.info(snpCount + " SNPs processed...");
            }

            VCFSampleRecord recordSampleA = vcfRecord.getSampleRecord(sampleA);
            VCFSampleRecord recordSampleB = vcfRecord.getSampleRecord(sampleB);

            //if the genotype of sample A and sample B differ
            //AND the genotype has been called for both samples
            //AND both genotypes are homozygous
            if (vcfRecord.getFilterValue().equals("PASS") &&
                    !recordSampleA.getGenotype().equals(recordSampleB.getGenotype())
                    && !recordSampleA.getGenotype().equals("./.") && !recordSampleB.getGenotype().equals("./.")
                    && !recordSampleA.isHeterozygous() && !recordSampleB.isHeterozygous()) {

                String alleleSampleA = vcfRecord.getReferenceAllele();
                String alleleSampleB = vcfRecord.getReferenceAllele();

                if(recordSampleA.isAlternative()){
                    alleleSampleA = vcfRecord.getAlternativeAllele();
                }

                if(recordSampleB.isAlternative()){
                    alleleSampleB = vcfRecord.getAlternativeAllele();
                }

                 pw.println(vcfRecord.getChromsome() + "\t" + ( vcfRecord.getPosition() - 1 ) + "\t" + vcfRecord.getPosition() + "\t" + alleleSampleA + vcfRecord.getPosition() + alleleSampleB);
                 pw.flush();

            }

        }

        pw.close();

    }

}
