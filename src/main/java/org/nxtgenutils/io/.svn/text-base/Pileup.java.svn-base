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
 * Time: 13:39:03
 */

/**
 * Interface to access information in a pileup record of a multi-pileup record in
 * a multi-pileup file generated with the samtools mpileup command.
 *
 * @author Michael Mueller
 */
public interface Pileup {

    /**
     * Returns the depth of coverage at the given sequence position.
     *
     * @return the depth of coverage
     */
    public int getCoverage();

    /**
     * Returns the base called at the given sequence position.
     *
     * @return the base calls
     */
    public String getBases();

    /**
     * Returns the number of base calls that match the reference sequence.
     *
     * @return the number of base calls
     */
    public int getReferenceMatchCount();

    /**
     * Returns the number of times the specified base has been called at the
     * given sequence position.
     *
     * @param base the base
     * @return the number of base calls
     */
    public int getBaseCallCount(String base);

    /**
     * Returns the consensus base call at the given sequence position.
     *
     * @return the consensus base call
     */
    public String getConsensusBaseCall();

    /**
     * Returns the consensus base call at the given position using
     * an allele frequency cutoff to call heterozygous positions.
     *
     * @param alleleFrequencyCutoff minimum frequency at which a base will
     *                              be considered a second allele at the
     *                              given position
     * @return the consensus base call
     */
    public String getConsensusBaseCall(double alleleFrequencyCutoff);

    /**
     * Returns the parent MultiPileup instance.
     *
     * @return the parent MultiPileup
     */
    public MultiPileup getPileupRecord();

    /**
     * Returns the ASCII encoded base call quality scores.
     *
     * @return the ASCII encoded quality scores
     */
    public String getQualities();

    /**
     * Returns true if any of the reads pileup at the given positions contains
     * evidence for an indel with respect to the reference sequence.
     *
     * @return true if there is evidence for an indel
     */
    public boolean hasIndels();

    /**
     * Returns the base called at the given sequence position omitting indel calls.
     *
     * @return base calls without indel calls
     */
    public String getBasesWithoutIndels();

}
