package org.nxtgenutils.bsseq;

import net.sf.samtools.SAMRecord;

import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;

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
 * Time: 17:00:04

/**
 * @author Michael Mueller
 */
public abstract class AbstractBSReadFilterSet extends HashSet<BSReadFilter> implements BSReadFilterSet  {

    private Map<String, Boolean> filterResults = new HashMap<String, Boolean>();

    public AbstractBSReadFilterSet() {
    }

    @Override
    public boolean passesFilterSet(SAMRecord record) {

        boolean retVal = true;

        for(BSReadFilter filter : this){

            boolean passesFilter = filter.passesFilter(record);
            filterResults.put(filter.toString(), passesFilter);
            if(!passesFilter){
                retVal = false;
            }

        }

        return retVal;
    }

    @Override
    public Map<String, Boolean> getFilterResults() {
        return filterResults;
    }
}