package org.nxtgenutils.io.impl;

import org.nxtgenutils.io.VCFSampleRecord;

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
 * Date: 06-Oct-2011
 * Time: 16:31:01
 */

/**
 * @author Michael Mueller
 */
public class VCFSampleRecordImpl implements VCFSampleRecord {

    private String sampleName;
    private String genotype;

    public VCFSampleRecordImpl(String sampleName, String genotype) {
        this.sampleName = sampleName;
        this.genotype = genotype;
    }

    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public String getGenotype() {
        return genotype;
    }

    public void setGenotype(String genotype) {
        this.genotype = genotype;
    }

    public boolean isHeterozygous(){
        return genotype.equals("0/1") || genotype.equals("1/0");
    }

    public boolean isReference(){
        return genotype.equals("0/0");
    }

    public boolean isAlternative(){
        return genotype.equals("1/1");
    }

    @Override
    public String toString() {
        return "VCFSampleRecord{" +
                "sampleName='" + sampleName + '\'' +
                ", genotype='" + genotype + '\'' +
                '}';
    }
}
