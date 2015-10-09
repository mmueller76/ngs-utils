package org.nxtgenutils.io;

import org.nxtgenutils.Strand;


import java.util.Set;

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
 * Date: 05-Jan-2012
 * Time: 17:41:59
 */

/**
 *  Interface to access information of a methylation profile record.
 *
 * @author Michael Mueller
 */
public interface MethylationProfileRecord {

    /**
     * Returns the sequence name of the profiled position.
     *
     * @return the sequence name
     */
    String getSequenceName();

    /**
     * Sets the sequence name of the profiled position.
     *
     * @param sequenceName the sequence name
     */
    void setSequenceName(String sequenceName);

    /**
     * Returns the start coordinate of the profiled position.
     *
     * @return the start coordinate
     */
    int getStartPosition();

    /**
     * Sets the start coordinate of the profiled position.
     *
     * @param startPosition the start coordinate
     */
    void setStartPosition(int startPosition);

    /**
     * Returns the end coordinate of the profiled position.
     *
     * @return the end coordinate
     */
    int getEndPosition();

    /**
     * Sets the end coordinate of the profiled position.
     *
     * @param endPosition the end coordinate
     */
    void setEndPosition(int endPosition);

    /**
     * Returns the strand of the profiled position.
     *
     * @return the strand
     */
    Strand getStrand();

    /**
     * Sets the strand of the profiled position.
     *
     * @param strand the strand
     */
    void setStrand(Strand strand);

    /**
     * Returns the reference nucleotide context of the profiled position.
     *
     * @return the reference nucleotide context
     */
    String getReferenceContext();

    /**
     * Sets the reference nucleotide context of the profile position.
     *
     * @param referenceContext the reference nucleotide context
     */
    void setReferenceContext(String referenceContext);

    /**
     * Returns true if the profiled position is repeat masked.
     *
     * @return true if the profiled position is repeat masked
     */
    boolean isRepeatMasked();

    /**
     * Sets wether or not the profiled position is repeat masked.
     *
     * @param repeatMasked wether or not the position is repeat masked
     */
    void setRepeatMasked(boolean repeatMasked);

    /**
     * Returns the sample records for the profiled position.
     *
     * @return set of profile records
     */
    Set<MethylationProfileSampleRecord> getSampleRecords();

    /**
     * Adds a sample record to the profile record.
     *
     * @param sampleRecord the sample record
     */
    void addSampleRecord(MethylationProfileSampleRecord sampleRecord);

    /**
     * Creates a string representation of the profile record.
     *
     * @return string representation of profile record
     */
    String formatProfileRecord();

    /**
     * Returns the number of samples in the record.
     *
     * @return the sample count
     */
    int getSampleCount();
}
