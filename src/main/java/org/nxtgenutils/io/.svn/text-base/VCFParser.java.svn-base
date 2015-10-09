package org.nxtgenutils.io;

import java.util.Set;
import java.util.Iterator;

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
 * Time: 16:07:57
 */

/**
 * Interface to access information in a VCF genotype formatted file.
 */
public interface VCFParser {

    /**
     * Returns the header lines of the VCF file.
     *
     * @return a set of header lines
     */
    Set<String> getHeaderLines();

    /**
     * Returns a VCFRecord iterator.
     *
     * @return the BedRecord iterator
     */
    Iterator<VCFRecord> iterator();

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
