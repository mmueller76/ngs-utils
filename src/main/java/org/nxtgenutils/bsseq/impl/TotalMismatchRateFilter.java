package org.nxtgenutils.bsseq.impl;

import org.nxtgenutils.bsseq.AbstractBSReadFilter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
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
 * Time: 13:22:10
 * To change this template use File | Settings | File Templates.
 *
 * Filters SAM records by the mismatch rate in the aligned part of
 * read, i.e. excluding clipped parts of the read.
 *
 */

public class TotalMismatchRateFilter extends AbstractBSReadFilter<Double> {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(TotalMismatchRateFilter.class);

    @Override
    public boolean passesFilter(SAMRecord bsRead) {

        //passes filter if read is not mapped
        if(bsRead.getReadUnmappedFlag()) return true;

        //get clipped read bases
        int clippedBases = 0;
        for(CigarElement element : bsRead.getCigar().getCigarElements()){

            if(element.getOperator() == CigarOperator.SOFT_CLIP ||
                element.getOperator() == CigarOperator.HARD_CLIP){
                    clippedBases = clippedBases + element.getLength();
            }

        }

        //get non-BS and BS-invalid mismatches
        int nonBSMismatchCount = bsRead.getIntegerAttribute("BN");
        int bsInvalidMismatchCount = bsRead.getIntegerAttribute("BI");
        double totalMismatchRate =  (double)(nonBSMismatchCount + bsInvalidMismatchCount)/(double)(bsRead.getReadLength() - clippedBases);

        return totalMismatchRate <= this.getFilterValue();

    }

    @Override
    public void setFilterValueFromString(String value) {
        try {
            this.setFilterValue(Double.parseDouble(value));
        } catch (Exception e) {
            logger.error("Unable to set filter value: " + value, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof TotalMismatchRateFilter)) return false;

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
