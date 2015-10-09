package org.nxtgenutils.bsseq;

import net.sf.samtools.SAMRecord;

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
 * Date: 24-Feb-2011
 * Time: 13:17:16
 */

/**
 * @author Michael Mueller
 * @param <V>
 */
public interface BSReadFilter <V> {

    public boolean passesFilter(SAMRecord bsRead);

    public V getFilterValue();

    public void setFilterValue(V value);

    public void setFilterValueFromString(String value);

}

