package org.nxtgenutils.io;

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
 * Date: 09-Nov-2011
 * Time: 13:44:01
 */

/**
 * Interface to access information in a sample record of a VCF record.
 *
 * @author Michael Mueller
 */
public interface VCFSampleRecord {

    /**
     * Returns the sample name.
     *
     * @return the sample name
     */
    String getSampleName();

    /**
     * Sets the sample name.
     *
     * @param sampleName the sample name
     */
    void setSampleName(String sampleName);

    /**
     * Returns the genotype.
     * 0/0 homozygous reference allele
     * 1/1 homozygous alternative allele
     * 0/1 heterozygous allele
     *
     * @return the genotype
     */
    String getGenotype();

    /**
     * Sets the genotype.
     *
     * @param genotype the genotype
     */
    void setGenotype(String genotype);

    /**
     * Returns true of the genotype of the sample at the given position
     * is heterozyogous.
     *
     * @return true of heterozyogous
     */
    boolean isHeterozygous();

    /**
     * Returns true if the sample has the reference genotype at the given position.
     *
     * @return true if sample has reference genotype
     */
    boolean isReference();

    /**
     * Returns true if the sample has the alternative genotype at the given position.
     *
     * @return true of the sample has the alternative genotype
     */
    boolean isAlternative();

}
