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
 * Date: 06-Jan-2012
 * Time: 13:46:16
 * To change this template use File | Settings | File Templates.
 */

/**
 * Interface to access information in a methylation profile sample record.
 *
 * @author Michael Mueller
 */
public interface MethylationProfileSampleRecord {

    /**
     * Returns the sample nucleotide context at the profiled position.
     *
     * @return the sample nucleotide context
     */
    String getSampleContext();

    /**
     * Sets the sample nucleotide context at the profiled position.
     *
     * @param sampleContext the sample nucleotide context
     */
    void setSampleContext(String sampleContext);

    /**
     * Returns the number of cytosine base calls at the the profiled position.
     *
     * @return the nuber of ctosin base calls
     */
    int getCCount();

    /**
     * Sets the number of cytosine base calls at the profiled position.
     *
     * @param cCount the nuber of ctosin base calls
     */
    void setCCount(int cCount);

    /**
     * Returns the number of cytosine and thymidine base calls at the the profiled position.
     *
     * @return the number of cytosine and thymidine base calls
     */
    int getCtCoverage();

    /**
     * Sets the number of cytosine and thymidine base calls at the the profiled position.
     *
     * @param ctCoverage the number of cytosine and thymidine base calls
     */
    void setCtCoverage(int ctCoverage);

    /**
     * Returns the number of base calls other than cytosine and thymidine at the the profiled position.
     *
     * @return the number of base calls
     */
    int getNonCtCoverage();

    /**
     * Returns the number of base calls other than cytosine and thymidine at the the profiled position.
     *
     * @param nonCtCoverage the number of base calls
     */
    void setNonCtCoverage(int nonCtCoverage);

    /**
     * Returns the methylation level in % at the profiled position.
     *
     * @return the methylation level in %
     */
    int getMethylationLevel();

    /**
     * Returns the SNP type code at the profiled position.
     *
     * @return 0 if there is no evidence for a SNP, -1 if there is insufficient data
     *         and values >0 specifying the SNP type.
     */
    int getSnpType();

    /**
     * Sets the SNP type code at the profiled position.
     *
     * @param snpType the SNP type code
     */
    void setSnpType(int snpType);

    /**
     * Returns whether or not there is evidence for an indel at the profiled position.
     *
     * @return true if there is evidence for an indel
     */
    boolean hasIndels();

    /**
     * Sets whether or not there is evidence for an indel at the profiled position.
     *
     * @param hasIndels boolean specifying if there is evidence for an indel
     */
    void setHasIndels(boolean hasIndels);

    /**
     * Returns the parent MethylationProfileRecord instance
     *
     * @return the parent MethylationProfileRecord instance
     */
    MethylationProfileRecord getMethylationProfileRecord();

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
     * Returns the score assigned to the profiled methylation level.
     *
     * @return the score
     */
    int getScore();

    /**
     * Sets the score assigned to the profiled methylation level.
     *
     * @param score the score
     */
    void setScore(int score);
}
