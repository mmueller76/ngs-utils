package org.nxtgenutils.io.impl;

import org.nxtgenutils.io.VCFSampleRecord;
import org.nxtgenutils.io.VCFRecord;

import java.util.Map;
import java.util.LinkedHashMap;
import java.util.NoSuchElementException;

/**
 * This file is part of NxtGenUtils.
 * <p/>
 * NxtGenUtils is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p/>
 * NxtGenUtils is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p/>
 * You should have received a copy of the GNU General Public License
 * along with NxtGenUtils.  If not, see <http://www.gnu.org/licenses/>.
 * <p/>
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 06-Oct-2011
 * Time: 16:29:45
 */

/**
 * @author Michael Mueller
 */
public class VCFRecordImpl implements VCFRecord {

    private String chromsome;
    private int position;
    private String referenceAllele;
    private String alternativeAllele;
    private Map<String, VCFSampleRecord> sampleRecords = new LinkedHashMap<String, VCFSampleRecord>();
    private String vcfString;
    private String filterValue;

    public String getChromsome() {
        return chromsome;
    }

    public void setChromsome(String chromsome) {
        this.chromsome = chromsome;
    }

    public int getPosition() {
        return position;
    }

    public void setPosition(int position) {
        this.position = position;
    }

    public String getReferenceAllele() {
        return referenceAllele;
    }

    public void setReferenceAllele(String referenceAllele) {
        this.referenceAllele = referenceAllele;
    }

    public String getAlternativeAllele() {
        return alternativeAllele;
    }

    public void setAlternativeAllele(String alternativeAllele) {
        this.alternativeAllele = alternativeAllele;
    }

    public Map<String, VCFSampleRecord> getSampleRecords() {
        return sampleRecords;
    }

    public void setSampleRecords(Map<String, VCFSampleRecord> sampleRecords) {
        this.sampleRecords = sampleRecords;
    }

    public VCFSampleRecord getSampleRecord(String sampleName) throws NoSuchElementException {
        if (this.sampleRecords.get(sampleName) != null) {
            return this.sampleRecords.get(sampleName);
        } else {
            throw new NoSuchElementException("VCF file does not contain sample " + sampleName);
        }
    }

    public void addSampleRecord(VCFSampleRecord sampleRecord) {
        this.sampleRecords.put(sampleRecord.getSampleName(), sampleRecord);
    }

    public String getVcfString() {
        return vcfString;
    }

    public void setVcfString(String vcfString) {
        this.vcfString = vcfString;
    }

    /**
     * Returns the value of the filter column.
     *
     * @return the filter value
     */
    @Override
    public String getFilterValue() {
        return filterValue;
    }

    /**
     * Sets the value of the filter column.
     *
     * @param filterValue the value
     */
    @Override
    public void setFilterValue(String filterValue) {
        this.filterValue=filterValue;
    }

    public String getSampleAllele(String sampleName) throws NoSuchElementException {

        if (this.sampleRecords.get(sampleName) != null) {
            VCFSampleRecord record = this.sampleRecords.get(sampleName);
            if (record.getGenotype().equals("0/0")) {
                return this.getReferenceAllele();
            } else if (record.getGenotype().equals("1/1")) {
                return this.getAlternativeAllele();
            } else if (record.getGenotype().equals("0/1") || record.getGenotype().equals("1/0")) {
                return this.getReferenceAllele() + "/" + this.getAlternativeAllele();
            } else {
                return ".";
            }
        } else {
            throw new NoSuchElementException("VCF file does not contain sample " + sampleName);
        }

    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof VCFRecordImpl)) return false;

        VCFRecordImpl vcfRecord = (VCFRecordImpl) o;

        if (position != vcfRecord.position) return false;
        if (!alternativeAllele.equals(vcfRecord.alternativeAllele)) return false;
        if (!chromsome.equals(vcfRecord.chromsome)) return false;
        if (!referenceAllele.equals(vcfRecord.referenceAllele)) return false;
        if (sampleRecords != null ? !sampleRecords.equals(vcfRecord.sampleRecords) : vcfRecord.sampleRecords != null)
            return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = chromsome.hashCode();
        result = 31 * result + position;
        result = 31 * result + referenceAllele.hashCode();
        result = 31 * result + alternativeAllele.hashCode();
        result = 31 * result + (sampleRecords != null ? sampleRecords.hashCode() : 0);
        return result;
    }


}
