package org.nxtgenutils.io;

import java.util.Map;
import java.util.NoSuchElementException;

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
 * Time: 13:44:39
 */

/**
 * Interface to access information in a VCF (variant calling format)
 * formatted record.
 *
 * @author Michael Mueller
 */
public interface VCFRecord {

    /**
     * Returns the chromosome.
     *
     * @return the chromosome name
     */
    String getChromsome();

    /**
     * Sets the chromosome.
     *
     * @param chromsome the chromosome name
     */
    void setChromsome(String chromsome);

    /**
     * Returns the sequence position.
     *
     * @return the sequence position
     */
    int getPosition();

    /**
     * Sets the sequence position.
     *
     * @param position the sequence position
     */
    void setPosition(int position);

    /**
     * Returns the reference allele.
     *
     * @return the reference allele
     */
    String getReferenceAllele();

    /**
     * Sets the reference allele.
     *
     * @param referenceAllele the reference allele
     */
    void setReferenceAllele(String referenceAllele);

    /**
     * Returns the alternative allele.
     *
     * @return the alternative allele
     */
    String getAlternativeAllele();

    /**
     * Sets the alternative allele.
     *
     * @param alternativeAllele the alternative allele
     */
    void setAlternativeAllele(String alternativeAllele);

    /**
     * Returns the sample records.
     *
     * @return a map of sample names and records
     */
    Map<String, VCFSampleRecord> getSampleRecords();

    /**
     * Sets the sample records
     *
     * @param sampleRecords a map of sample names and records
     */
    void setSampleRecords(Map<String, VCFSampleRecord> sampleRecords);

    /**
     * Returns the record for the specified sample.
     *
     * @param sampleName the name of the sample
     * @return the sample record
     * @throws NoSuchElementException if the specified sample name does not exist
     */
    VCFSampleRecord getSampleRecord(String sampleName) throws NoSuchElementException;

    /**
     * Adds a sample record.
     *
     * @param sampleRecord the sample record
     */
    void addSampleRecord(VCFSampleRecord sampleRecord);

    /**
     * Returns the VCF formatted string representation of the record.
     *
     * @return the VCF string
     */
    String getVcfString();

    /**
     * Sets the VCF formatted string representation of the record.
     *
     * @param vcfString the VCF string
     */
    void setVcfString(String vcfString);

    /**
     * Returns the allele of the specified sample
     *
     * @param sampleName the sample name
     * @return the allele
     * @throws NoSuchElementException if the specified sample name does not exist
     */
    String getSampleAllele(String sampleName) throws NoSuchElementException;

    /**
     * Returns the value of the filter column.
     *
     * @return the filter value
     */
    String getFilterValue();


    /**
     * Sets the value of the filter column.
     *
     * @param filterValue the value
     */
    void setFilterValue(String filterValue);
    
}
