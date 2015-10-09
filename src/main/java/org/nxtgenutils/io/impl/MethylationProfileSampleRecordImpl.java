package org.nxtgenutils.io.impl;

import org.nxtgenutils.io.MethylationProfileRecord;
import org.nxtgenutils.io.MethylationProfileSampleRecord;

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
 * Time: 15:25:10
 */

/**
 * Implementation of MethylationProfileSampleRecord interface.
 *
 * @author Michael Mueller
 */
public class MethylationProfileSampleRecordImpl implements MethylationProfileSampleRecord {

    private String sampleContext;
    private int cCount;
    private int ctCoverage;
    private int nonCtCoverage;
    private int snpType = -1;
    private boolean hasIndels = false;
    private MethylationProfileRecord methylationProfileRecord;
    private String sampleName = "";
    private int score = -1;

    /**
     * Constructs a MethylationProfileSampleRecord instance.
     *
     * @param parentRecord the parent MethylationProfileRecord instance
     * @param sampleContext the sample nucleotide context of the profiled position
     * @param cCount the number of cytosine base calls at the profiled position
     * @param ctCoverage the number of cytosine and thymidine base calls at the profiled position
     * @param nonCtCoverage the number of non-cytosine/thymidine base calls at the profiled position
     */
    public MethylationProfileSampleRecordImpl(MethylationProfileRecord parentRecord, String sampleContext, int cCount, int ctCoverage, int nonCtCoverage) {
        this.methylationProfileRecord = parentRecord;
        this.sampleContext = sampleContext;
        this.cCount = cCount;
        this.ctCoverage = ctCoverage;
        this.nonCtCoverage = nonCtCoverage;
    }

    /**
     * Constructs a MethylationProfileSampleRecord instance.
     *
     * @param parentRecord the parent MethylationProfileRecord instance
     * @param sampleContext the sample nucleotide context of the profiled position
     * @param cCount the number of cytosine base calls at the profiled position
     * @param ctCoverage the number of cytosine and thymidin base calls at the profiled position
     * @param nonCtCoverage the number of non-cytosine/thymidine base calls at the profiled position
     * @param sampleName the sample name
     */
    public MethylationProfileSampleRecordImpl(MethylationProfileRecord parentRecord, String sampleContext, int cCount, int ctCoverage, int nonCtCoverage, String sampleName) {
        this.sampleName = sampleName;
        this.methylationProfileRecord = parentRecord;
        this.sampleContext = sampleContext;
        this.cCount = cCount;
        this.ctCoverage = ctCoverage;
        this.nonCtCoverage = nonCtCoverage;
    }

    /**
     * Constructs a MethylationProfileSampleRecord instance.
     *
     * @param parentRecord the parent MethylationProfileRecord instance
     * @param sampleContext the sample nucleotide context of the profiled position
     * @param cCount the number of cytosine base calls at the profiled position
     * @param ctCoverage the number of cytosine and thymidin base calls at the profiled position
     * @param nonCtCoverage the number of non-cytosine/thymidine base calls at the profiled position
     * @param sampleName the sample name
     * @param score the score assigned to the methylation call
     * @param snpType the SNP type code
     * @param hasIndels whether or not there is effidence for an indel at the profiled position
     */
    public MethylationProfileSampleRecordImpl(MethylationProfileRecord parentRecord, String sampleContext, int cCount, int ctCoverage, int nonCtCoverage, String sampleName, int score, int snpType, boolean hasIndels) {
        this.sampleName = sampleName;
        this.methylationProfileRecord = parentRecord;
        this.sampleContext = sampleContext;
        this.cCount = cCount;
        this.ctCoverage = ctCoverage;
        this.nonCtCoverage = nonCtCoverage;
        this.score=score;
        this.snpType=snpType;
        this.hasIndels=hasIndels;
    }

    /**
     * Returns the sample nucleotide context at the profiled position.
     *
     * @return the sample nucleotide context
     */
    @Override
    public String getSampleContext() {
        return sampleContext;
    }

    /**
     * Sets the sample nucleotide context at the profiled position.
     *
     * @param sampleContext the sample nucleotide context
     */
    @Override
    public void setSampleContext(String sampleContext) {
        this.sampleContext = sampleContext;
    }

    /**
     * Returns the number of cytosine base calls at the the profiled position.
     *
     * @return the nuber of ctosin base calls
     */
    @Override
    public int getCCount() {
        return cCount;
    }

    /**
     * Sets the number of cytosine base calls at the profiled position.
     *
     * @param cCount the nuber of ctosin base calls
     */
    @Override
    public void setCCount(int cCount) {
        this.cCount = cCount;
    }

    /**
     * Returns the number of cytosine and thymidine base calls at the the profiled position.
     *
     * @return the number of cytosine and thymidine base calls
     */
    @Override
    public int getCtCoverage() {
        return ctCoverage;
    }

    /**
     * Sets the number of cytosine and thymidine base calls at the the profiled position.
     *
     * @param ctCoverage the number of cytosine and thymidine base calls
     */
    @Override
    public void setCtCoverage(int ctCoverage) {
        this.ctCoverage = ctCoverage;
    }

    /**
     * Returns the number of base calls other than cytosine and thymidine at the the profiled position.
     *
     * @return the number of base calls
     */
    @Override
    public int getNonCtCoverage() {
        return nonCtCoverage;
    }

    /**
     * Returns the number of base calls other than cytosine and thymidine at the the profiled position.
     *
     * @param nonCtCoverage the number of base calls
     */
    @Override
    public void setNonCtCoverage(int nonCtCoverage) {
        this.nonCtCoverage = nonCtCoverage;
    }

    /**
     * Returns the methylation level in % at the profiled position.
     *
     * @return the methylation level in %
     */
    @Override
    public int getMethylationLevel() {
        return calculateMethylationLevel(cCount, (ctCoverage - cCount));
    }

    /**
     * Returns the SNP type code at the profiled position.
     *
     * @return 0 if there is no evidence for a SNP, -1 if there is insufficient data
     *         and values >0 specifying the SNP type.
     */
    @Override
    public int getSnpType() {
        return snpType;
    }

    /**
     * Sets the SNP type code at the profiled position.
     *
     * @param snpType the SNP type code
     */
    @Override
    public void setSnpType(int snpType) {
        this.snpType = snpType;
    }

    /**
     * Returns whether or not there is evidence for an indel at the profiled position.
     *
     * @return true if there is evidence for an indel
     */
    @Override
    public boolean hasIndels() {
        return hasIndels;
    }

    /**
     * Sets whether or not there is evidence for an indel at the profiled position.
     *
     * @param hasIndels boolean specifying if there is evidence for an indel
     */
    @Override
    public void setHasIndels(boolean hasIndels) {
        this.hasIndels = hasIndels;
    }

    /**
     * Returns the parent MethylationProfileRecord instance
     *
     * @return the parent MethylationProfileRecord instance
     */
    @Override
    public MethylationProfileRecord getMethylationProfileRecord() {
        return methylationProfileRecord;
    }

    /**
     * Returns the sample name.
     *
     * @return the sample name
     */
    @Override
    public String getSampleName() {
        return sampleName;
    }

    /**
     * Sets the sample name.
     *
     * @param sampleName the sample name
     */
    @Override
    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    /**
     * Returns the score assigned to the profiled methylation level.
     *
     * @return the score
     */
    @Override
    public int getScore() {
        return score;
    }

    /**
     * Sets the score assigned to the profiled methylation level.
     *
     * @param score the score
     */
    @Override
    public void setScore(int score) {
        this.score = score;
    }

    /**
     * Calculates the % methylation = C count/(C count + T count) *100.
     *
     * @param cCount the number of Cs
     * @param tCount the number of Ts
     * @return the methylation level in percent
     */
    private int calculateMethylationLevel(int cCount, int tCount) {
        return (int) Math.round((double) cCount / (double) (cCount + tCount) * 100.0);
    }

}
