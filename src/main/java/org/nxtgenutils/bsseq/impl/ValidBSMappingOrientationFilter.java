package org.nxtgenutils.bsseq.impl;

import org.nxtgenutils.bsseq.AbstractBSReadFilter;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;

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
 * Time: 13:31:45
 */

/**
 * @author Michael Mueller
 */
public class ValidBSMappingOrientationFilter extends AbstractBSReadFilter<Boolean> {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(ValidBSMappingOrientationFilter.class);


    @Override
    public boolean passesFilter(SAMRecord bsRead) {

        //passes filter if read is not mapped
        if(bsRead.getReadUnmappedFlag()) return true;

        //if valid mapping orientation not required
        //read passes filter
        if(this.getFilterValue()  == false) return true;

        //else return true if read has valid BS mapping orientation
        return bsRead.getCharacterAttribute("BO") == 'V';

    }

    @Override
    public void setFilterValueFromString(String value) {
        try {
            this.setFilterValue(Boolean.parseBoolean(value));
        } catch (Exception e) {
            logger.error("Unable to set filter value: " + value, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof ValidBSMappingOrientationFilter)) return false;

        return super.equals(o);
    }

    @Override
    public int hashCode() {
        return super.hashCode();
    }

    @Override
    public String toString(){

           return this.getClass().toString().replace("class org.nxtgenutils.impl.","") + " = " + this.getFilterValue();

    }

}
