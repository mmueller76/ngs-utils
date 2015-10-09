package org.nxtgenutils.impl;

import org.nxtgenutils.Strand;
import org.nxtgenutils.io.VCFRecord;

import java.io.File;

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
 * Date: 13-Jan-2012
 * Time: 11:07:57
 */

/**
 * @author Michael Mueller
 */
public class BSSingleEndReadPhaser extends SingleEndReadPhaser {

    private Strand genomicTemplateStrand;

    public BSSingleEndReadPhaser(File mappingInput, File genotypeVCF, Strand genomicTemplateStrand) {
        super(mappingInput, genotypeVCF);
        this.genomicTemplateStrand=genomicTemplateStrand;
    }

    @Override
    protected String determineSample(VCFRecord record, String allele, String sampleA, String sampleB) {

        String retVal = null;
        //get sample alleles
        String alleleSampleA;
        String alleleSampleB;
        if(record.getSampleRecord(sampleA).isReference()){
            alleleSampleA = record.getReferenceAllele();
        } else {
            alleleSampleA = record.getAlternativeAllele();
        }

        if(record.getSampleRecord(sampleB).isReference()){
            alleleSampleB = record.getReferenceAllele();
        } else {
            alleleSampleB = record.getAlternativeAllele();
        }

        //in case of a C <-> T SNP on we cannot determine the sample
        //on the bisulphite converted forward strand in case of a
        //G<->A SNP we cannot determine the sample on the bisulphite
        //converted reverse strand
        if((genomicTemplateStrand == Strand.FORWARD &&
                (alleleSampleA.equalsIgnoreCase("C") && alleleSampleB.equalsIgnoreCase("T") ||
                alleleSampleA.equalsIgnoreCase("T") && alleleSampleB.equalsIgnoreCase("C"))) ||
            (genomicTemplateStrand == Strand.REVERSE &&
                (alleleSampleA.equalsIgnoreCase("G") && alleleSampleB.equalsIgnoreCase("A") ||
                alleleSampleA.equalsIgnoreCase("A") && alleleSampleB.equalsIgnoreCase("G")))){

            return retVal;

        }

        //if we are not dealing with a C on the forward strand or G on the reverse strand
        //we don't need to do anything special
        if((genomicTemplateStrand == Strand.FORWARD && !(alleleSampleA.equalsIgnoreCase("C") || alleleSampleB.equalsIgnoreCase("C"))) ||
           (genomicTemplateStrand == Strand.REVERSE && !(alleleSampleA.equalsIgnoreCase("G") || alleleSampleB.equalsIgnoreCase("G")))){

            return super.determineSample(record, allele, sampleA, sampleB);

        }

        if(genomicTemplateStrand == Strand.FORWARD){

            if(alleleSampleA.equalsIgnoreCase("C") && (allele.equalsIgnoreCase("C") || allele.equalsIgnoreCase("T"))){
                retVal = sampleA;
            } else if (alleleSampleB.equalsIgnoreCase("C") && (allele.equalsIgnoreCase("C") || allele.equalsIgnoreCase("T"))){
                retVal = sampleB;
            } else {
                retVal = null;
            }

        } else {

            if(alleleSampleA.equalsIgnoreCase("G") && (allele.equalsIgnoreCase("G") || allele.equalsIgnoreCase("A"))){
                retVal = sampleA;
            } else if (alleleSampleB.equalsIgnoreCase("G") && (allele.equalsIgnoreCase("G") || allele.equalsIgnoreCase("A"))){
                retVal = sampleB;
            } else {
                retVal = null;
            }

        }

        return retVal;

    }

}
