package org.nxtgenutils.io;

import java.util.List;
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
 * Date: 18-Nov-2011
 * Time: 16:21:47
 */

/**
 * An interface to access in formation in a methylation profile.
 *
 * @author Michael Mueller
 */
public interface MethylationProfileParser {

    /**
     * Returns the line containig the column
     * headers.
     *
     * @return the column header line
     */
    String getHeaderLine();

    /**
     * Returns the names of the samples in the order they
     * appear in the profile file.
     *
     * @return a list of sample names
     */
    List<String> getSampleNames();

    /**
     * Sets the names of the samples in the profile file.
     * The order of the names must be the same as the order
     * of the samples in the profile.
     *
     * @param sampleNames a list of sample names
     */
    void setSampleNames(List<String> sampleNames);

    /**
     * Returns the number of samples in the profile.
     *
     * @return the number of samples
     */
    int getSampleCounts();

    /**
     * Returns iterator for methylation profile records.
     *
     * @return the record iterator
     */
    Iterator<MethylationProfileRecord> iterator();

}
