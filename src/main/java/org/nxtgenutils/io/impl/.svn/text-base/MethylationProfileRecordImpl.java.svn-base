package org.nxtgenutils.io.impl;

import org.nxtgenutils.Strand;
import org.nxtgenutils.io.MethylationProfileRecord;
import org.nxtgenutils.io.MethylationProfileSampleRecord;

import java.util.Set;
import java.util.LinkedHashSet;

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
 * Date: 03-Nov-2011
 * Time: 15:24:50
 */

/**
 * Implementation of the MethylationProfileRecord interface.
 *
 * @author Michael Mueller
 */
public class MethylationProfileRecordImpl implements MethylationProfileRecord {

    private String sequenceName;
    private int startPosition;
    private int endPosition;
    private Strand strand;
    private String referenceContext;
    private boolean isRepeatMasked;
    private static String TAB = "\t";

    static private String forwardStrandString = "+";
    static private String reverseStrandString = "-";

    private Set<MethylationProfileSampleRecord> sampleRecords = new LinkedHashSet<MethylationProfileSampleRecord>();

    /**
     * Constructs a MethylationProfileRecord instance for a profiled sequence position.
     *
     * @param sequenceName the name of the sequence
     * @param startPosition the start coordinate of the profiled position
     * @param endPosition the end coordinate of the profiled position
     * @param strand the strand of the profiled position
     * @param referenceContext the reference sequence nucleotide context of the profiled position
     * @param repeatMasked wether or not the position is repeat masked
     */
    public MethylationProfileRecordImpl(String sequenceName, int startPosition, int endPosition, Strand strand, String referenceContext, boolean repeatMasked) {
        this.sequenceName = sequenceName;
        this.startPosition = startPosition;
        this.endPosition = endPosition;
        this.strand = strand;
        this.referenceContext = referenceContext;
        isRepeatMasked = repeatMasked;
    }

    /**
     * Returns the sequence name of the profiled position.
     *
     * @return the sequence name
     */
    @Override
    public String getSequenceName() {
        return sequenceName;
    }

    /**
     * Sets the sequence name of the profiled position.
     *
     * @param sequenceName the sequence name
     */
    @Override
    public void setSequenceName(String sequenceName) {
        this.sequenceName = sequenceName;
    }

    /**
     * Returns the start coordinate of the profiled position.
     *
     * @return the start coordinate
     */
    @Override
    public int getStartPosition() {
        return startPosition;
    }

    /**
     * Sets the start coordinate of the profiled position.
     *
     * @param startPosition the start coordinate
     */
    @Override
    public void setStartPosition(int startPosition) {
        this.startPosition = startPosition;
    }

    /**
     * Returns the end coordinate of the profiled position.
     *
     * @return the end coordinate
     */
    @Override
    public int getEndPosition() {
        return endPosition;
    }

    /**
     * Sets the end coordinate of the profiled position.
     *
     * @param endPosition the end coordinate
     */
    @Override
    public void setEndPosition(int endPosition) {
        this.endPosition = endPosition;
    }

    /**
     * Returns the strand of the profiled position.
     *
     * @return the strand
     */
    @Override
    public Strand getStrand() {
        return strand;
    }

    /**
     * Sets the strand of the profiled position.
     *
     * @param strand the strand
     */
    @Override
    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    /**
     * Returns the reference nucleotide context of the profiled position.
     *
     * @return the reference nucleotide context
     */
    @Override
    public String getReferenceContext() {
        return referenceContext;
    }

    /**
     * Sets the reference nucleotide context of the profile position.
     *
     * @param referenceContext the reference nucleotide context
     */
    @Override
    public void setReferenceContext(String referenceContext) {
        this.referenceContext = referenceContext;
    }

    /**
     * Returns true if the profiled position is repeat masked.
     *
     * @return true if the profiled position is repeat masked
     */
    @Override
    public boolean isRepeatMasked() {
        return isRepeatMasked;
    }

    /**
     * Sets wether or not the profiled position is repeat masked.
     *
     * @param repeatMasked wether or not the position is repeat masked
     */
    @Override
    public void setRepeatMasked(boolean repeatMasked) {
        isRepeatMasked = repeatMasked;
    }

    /**
     * Returns the sample records for the profiled position.
     *
     * @return set of profile records
     */
    @Override
    public Set<MethylationProfileSampleRecord> getSampleRecords() {
        return sampleRecords;
    }

    /**
     * Adds a sample record to the profile record.
     *
     * @param sampleRecord the sample record
     */
    @Override
    public void addSampleRecord(MethylationProfileSampleRecord sampleRecord) {
        this.sampleRecords.add(sampleRecord);
    }

    /**
     * Creates a string representation of the profile record.
     *
     * @return string representation of profile record
     */
    @Override
    public String formatProfileRecord() {

        StringBuffer retVal = new StringBuffer();

        String strandString = forwardStrandString;
        if (strand == Strand.REVERSE) {
            strandString = reverseStrandString;
        }

        int isRepeatMaskedInt = 0;
        if(isRepeatMasked){
            isRepeatMaskedInt = 1;
        }

        retVal.append(sequenceName).append(TAB)
                .append(startPosition).append(TAB)
                .append(endPosition).append(TAB)
                .append(strandString).append(TAB)
                .append(referenceContext).append(TAB)
                .append(isRepeatMaskedInt);

        for (MethylationProfileSampleRecord record : this.getSampleRecords()) {

            int hasIndels = 0;
            if(record.hasIndels()){
                hasIndels = 1;
            }

            retVal.append(TAB).append(record.getSampleContext()).append(TAB)
                    .append(record.getCCount()).append(TAB)
                    .append(record.getCtCoverage()).append(TAB)
                    .append(record.getNonCtCoverage()).append(TAB)
                    .append(record.getMethylationLevel()).append(TAB)
                    .append(record.getScore()).append(TAB)
                    .append(record.getSnpType()).append(TAB)
                    .append(hasIndels);

        }

        return retVal.toString();

    }

    /**
     * Returns the number of samples in the record.
     *
     * @return the sample count
     */
    @Override
    public int getSampleCount(){
        return sampleRecords.size();
    }

}
