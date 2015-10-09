package org.nxtgenutils.bsseq;

import net.sf.samtools.SAMRecord;

import java.util.Set;
import java.util.Map;

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
 * Date: 03-Feb-2011
 * Time: 13:06:39
 */

/**
 * @author Michael Mueller
 */
public interface BSReadFilterSet extends Set<BSReadFilter> {

    boolean passesFilterSet(SAMRecord record);

    Map<String, Boolean> getFilterResults();

}
