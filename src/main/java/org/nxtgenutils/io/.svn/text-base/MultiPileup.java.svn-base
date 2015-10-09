package org.nxtgenutils.io;

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
 * Date: 09-Nov-2011
 * Time: 13:37:02
 */

/**
 * Interface to access information of a multi-pileup
 * record in a multi-pileup file generated with the
 * samtools mpileup command.
 *
 * @author Michael Mueller
 */
public interface MultiPileup {

    /**
     * Adds a pileup.
     *
     * @param pileup the pileup record.
     */
    void addPileup(Pileup pileup);

    /**
     * Returns the number of pileups in the record.
     *
     * @return the pileup count
     */
    int getPileupCount();

    /**
     * Returns the pileups in the record.
     *
     * @return a set of pileup records
     */
    Set<Pileup> getPileups();

    /**
     * Return the name of the reference sequence.
     *
     * @return the sequence name
     */
    String getSequenceName();

    /**
     *  Returns the position in the refernce sequence.
     *
     * @return the sequence position
     */
    int getPosition();

    /**
     * Returns the reference base at the position.
     *
     * @return the reference base
     */
    String getReferenceBase();

    /**
     * Returns the mpileup record as a string.
     *
     * @return the record string
     */
    String getRecordString();

    /**
     * Sets the mpilueup record string.
     *
     * @param recordString the record string
     */
    void setRecordString(String recordString);

    /**
     * Returns true if the position is repeat masked.
     *
     * @return true if the position is repeat masked
     */
    boolean isRepeatMasked();
    
}
