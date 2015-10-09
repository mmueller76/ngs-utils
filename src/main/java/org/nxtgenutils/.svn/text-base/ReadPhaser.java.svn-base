package org.nxtgenutils;

import java.io.File;

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
 * Date: 06-Jan-2012
 * Time: 16:54:59
 * To change this template use File | Settings | File Templates.
 */
public interface ReadPhaser {

    /**
     * Phases reads according to genotype information for sample A and sample B in
     * the VCF file. The names of the output BAM files will be generated automatically
     * from the input filename.
     *
     * @param sampleA the name of sample A
     * @param sampleB the name of sample B
     */
    void phase(String sampleA, String sampleB);

    /**
     * Phases reads according to genotype information for sample A and sample B in
     * the VCF file.
     *
     * @param sampleA the name of sample A
     * @param sampleB the name of sample B
     * @param outputSampleA path to BAM file for reads matching genotype of sample A
     * @param outputSampleB path to BAM file for reads matching genotype of sample B
     * @param outputAmbiguous path to BAM file for reads with ambiguous genotype information
     * @param outputPhasingSNPs path to VCF file for SNPs used to phase reads
     * @param outputRead2SNP2Sample path to TSV file for read to SNP to sample mapping
     * @param outputRead2SampleSummary
     */
    void phase(String sampleA, String sampleB, File outputSampleA, File outputSampleB, File outputAmbiguous, File outputPhasingSNPs, File outputRead2SNP2Sample, File outputRead2SampleSummary);
                           
    /**
     * Returns the minimum minor allele frequency at which a SNP called as heterozygous
     * in the VCF file will be returned as heterozyogous.
     *
     * @param cutoff the relative frequency
     */
    void setHeterozygousMinorAlleleFrequencyCutoff(double cutoff);

    /**
     * Sets the minimum minor allele frequency at which a SNP called as heterozygous
     * in the VCF file will be returned as heterozyogous.
     *
     * @return the relative frequency
     */
    double getHeterozygousMinorAlleleFrequencyCutoff();

}
