package org.nxtgenutils.bsseq.impl;

import net.sf.samtools.SAMRecord;
import org.nxtgenutils.bsseq.BSPairingStrandSpecificity;
import org.nxtgenutils.bsseq.AbstractBSReadFilter;
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
 * Time: 13:20:47
 */

/**
 * @author Michael Mueller
 */
public class BSPairingStrandSpecificityFilter extends AbstractBSReadFilter<BSPairingStrandSpecificity> {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(BSPairingStrandSpecificityFilter.class);


    @Override
    public boolean passesFilter(SAMRecord bsRead) {

        //passes filter if read is not mapped
        if(bsRead.getReadUnmappedFlag()) return true;

        //passes filter if not paired read
        if(!bsRead.getReadPairedFlag()) return true;

        //passes filer if filter value set to BOTH strands
        if(this.getFilterValue() == BSPairingStrandSpecificity.BOTH) return true;

        //if filter value is set to UNIQUE
        //get mapping template strand
        char mappingTemplateStrand = bsRead.getCharacterAttribute("BT");
        
        //...and pairing strand specificity
        char properPairingStrand = bsRead.getCharacterAttribute("BP");

        //...and mapping strand specificity
        char mappingStrandSpecificity = bsRead.getCharacterAttribute("BS");

        //passes filter if read is not properly paired on either strand
        //but maps only to one strand
        if(properPairingStrand == 'N' && mappingStrandSpecificity != 'B') return true;

        //return true if mapping template strand and pairing strand specificity are equal
        return (mappingTemplateStrand == properPairingStrand);

    }

    @Override
    public void setFilterValueFromString(String value) {
        try {
            this.setFilterValue(BSPairingStrandSpecificity.valueOf(value.toUpperCase()));
        } catch (Exception e) {
            logger.error("Unable to set filter value: " + value, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BSPairingStrandSpecificityFilter)) return false;

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

