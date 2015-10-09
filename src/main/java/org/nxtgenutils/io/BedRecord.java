package org.nxtgenutils.io;

import org.nxtgenutils.Strand;

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
 * Time: 13:30:47
 */

/**
 * Interface for accessing information in a BED formatted record.
 *
 * @author Michael Mueller
 */
public interface BedRecord {

    /**
     * Returns chromosome.
     *
     * @return the chromosome name
     */
    String getChromosome();

    /**
     * Sets chromosome.
     *
     * @param chromosome the chromosome name
     */
    void setChromosome(String chromosome);

    /**
     * Returns feature start (0 based).
     *
     * @return the feature start
     */
    int getStart();

    /**
     * Sets feature start (0 based).
     *
     * @param start the feature start
     */
    void setStart(int start);

    /**
     * Returns feature end (1 based).
     *
     * @return the feature end
     */
    int getEnd();

    /**
     * Sets feature end (1 based).
     *
     * @param end the feature end
     */
    void setEnd(int end);

    /**
     * Returns feature name.
     *
     * @return the feature name
     */
    String getName();

    /**
     * Sets feature name.
     *
     * @param name the feature name
     */
    void setName(String name);

    /**
     * Returns feature score.
     *
     * @return the score
     */
    int getScores();

    /**
     * Sets feature score.
     *
     * @param scores the score
     */
    void setScores(int scores);

    /**
     * Returns feature strand.
     *
     * @return the strand
     */
    Strand getStrand();

    /**
     * Sets feature strand.
     *
     * @param strand the strand
     */
    void setStrand(Strand strand);
    
}
