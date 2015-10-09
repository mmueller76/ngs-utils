package org.nxtgenutils.impl;

import org.nxtgenutils.ReadPhaser;
import org.nxtgenutils.io.VCFRecord;
import org.apache.log4j.Logger;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.AlignmentBlock;

import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.io.*;

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
 * Date: 09-Jan-2012
 * Time: 16:45:30
 */

/**
 * @author Michael Mueller
 */
public abstract class AbstractReadPhaser implements ReadPhaser {

    protected File mappingInput;
    protected File mappingIndex;
    protected File genotypeVCF;
    protected double minorAlleleFrequencyCutoff = -1;

    private static Logger logger = Logger.getLogger(AbstractReadPhaser.class);


    public AbstractReadPhaser(File genotypeVCF, File mappingInput) {
        this.genotypeVCF = genotypeVCF;
        this.mappingInput = mappingInput;
    }

    protected Map<SAMRecord, Integer> getOverlappingReads(VCFRecord vcfRecord) {

        Map<SAMRecord, Integer> retVal = new HashMap<SAMRecord, Integer>();

        //get all reads that overlap with the SNP
        SAMFileReader samFileReader = new SAMFileReader(mappingInput, mappingIndex);
        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        String chr = vcfRecord.getChromsome();
        int snpPos = vcfRecord.getPosition();
        Iterator<SAMRecord> it = samFileReader.queryOverlapping(chr, snpPos, snpPos);

        //for each read...
        while (it.hasNext()) {

            SAMRecord samRecord = it.next();

            //...determine the position of the SNP in the read
            int readOffset = -1;

//                    if (samRecord.getCigarString().equals("100M")) {
//
//                        readOffset = ReadUtils.getReadCoordinateForReferenceCoordinate(samRecord, snpPos);

//                        System.out.println(samRecord.getReadName() + " " + snpPos + " " + readOffset + " " + samRecord.getReadString().substring(readOffset, readOffset+1));
//                        System.out.println(samRecord.getSAMString());

//                    } else if (!samRecord.getCigarString().contains("N")) {
//
//                        readOffset = ReadUtils.getReadCoordinateForReferenceCoordinate(samRecord, snpPos);
//
//                        System.out.println(samRecord.getReadName() + " " + snpPos + " " + readOffset + " " + samRecord.getReadString().substring(readOffset, readOffset+1));
//                        System.out.println(samRecord.getSAMString());

            //} else {

            for (AlignmentBlock b : samRecord.getAlignmentBlocks()) {

                int referenceStart = b.getReferenceStart();
                int length = b.getLength();
                int referenceEnd = referenceStart + length - 1;

                if (snpPos >= referenceStart && snpPos <= referenceEnd) {

                    int snpOffset = snpPos - referenceStart;
                    int readStart = b.getReadStart();
                    readOffset = readStart + snpOffset - 1;

                }

            }

            //...if the SNP position has been found (i.e. is not in intronic sequence span by the read)
            if (readOffset != -1) {

                samRecord.setReadPairedFlag(true);
                retVal.put(samRecord, readOffset);

            }

        }

        //close BAM file
        samFileReader.close();


        return retVal;
    }

    /**
     * @param vcfRecord
     * @param overlappingReads
     * @param sampleA
     * @param sampleB
     * @param read2SNP2Sample
     */
    protected void assignReadsToSample(VCFRecord vcfRecord, Map<SAMRecord, Integer> overlappingReads, String sampleA, String sampleB, Map<SAMRecord, Map<VCFRecord, String>> read2SNP2Sample) {

        //...get the SNP allele in the read
        for (SAMRecord samRecord : overlappingReads.keySet()) {

            int readOffset = overlappingReads.get(samRecord);

            String allele = samRecord.getReadString().substring(readOffset, readOffset + 1);

            //if the read has not been
            //seen before
            if (!read2SNP2Sample.containsKey(samRecord)) {
                read2SNP2Sample.put(samRecord, new HashMap<VCFRecord, String>());
            }

            //...determine the sample from the allele
            String sample = determineSample(vcfRecord, allele, sampleA, sampleB);

            //if sample could not be determined
            //because there was a discrepency between
            //SNP information and read information
            //because of e.g. sequencing errors
            //skip read

            if (sample != null) {

                //...and store the read-to-SNP-association
                read2SNP2Sample.get(samRecord).put(vcfRecord, sample);

            }

        }

    }

    /**
     * @param vcfRecord
     * @param overlappingReads
     * @param sampleA
     * @param sampleB
     * @param pwRead2SNP2Sample
     */
    protected void writeReadsToSampleInformation(VCFRecord vcfRecord, Map<SAMRecord, Integer> overlappingReads, String sampleA, String sampleB, PrintWriter pwRead2SNP2Sample) {

        //...get the SNP allele in the read
        for (SAMRecord samRecord : overlappingReads.keySet()) {

            int readOffset = overlappingReads.get(samRecord);

            String allele = samRecord.getReadString().substring(readOffset, readOffset + 1);

            //if the read has not been
            //seen before
//            if (!pwRead2SNP2Sample.containsKey(samRecord)) {
//                pwRead2SNP2Sample.put(samRecord, new HashMap<VCFRecord, String>());
//            }

            //...determine the sample from the allele
            String sample = determineSample(vcfRecord, allele, sampleA, sampleB);


            //if sample could not be determined
            //because there was a discrepency between
            //SNP information and read information
            //because of e.g. sequencing errors
            //skip read
            //if (sample != null) {
            //...and store the read-to-SNP-association
            pwRead2SNP2Sample.println(samRecord.getReadName() + "\t" +
                    vcfRecord.getChromsome() + "\t" +
                    vcfRecord.getPosition() + "\t" +
                    vcfRecord.getSampleAllele(sampleA) + "/" + vcfRecord.getSampleAllele(sampleB) + "\t" +
                    allele + "\t" +
                    sample);

            pwRead2SNP2Sample.flush();
            //}

        }

    }

    protected void writeReadToSampleInformationSummary(Map<String, Integer[]> read2SampleCount, File outputRead2SampleSummary, String sampleA, String sampleB) {

        try {

            PrintWriter pw = new PrintWriter(outputRead2SampleSummary);

            //print header
            pw.println("read_name\t" + sampleA + "_allele\t" + sampleA + "_allele\t" + "sample_assignment");

            for (String readName : read2SampleCount.keySet()) {

                int snpCountSampleA = read2SampleCount.get(readName)[0];
                int snpCountSampleB = read2SampleCount.get(readName)[1];
                String sample = null;

                if (snpCountSampleA > snpCountSampleB) {

                    sample = sampleA;

                } else if (snpCountSampleA < snpCountSampleB) {

                    sample = sampleB;

                } else {

                    sample = "null";

                }

                pw.println(readName + "\t" + snpCountSampleA + "\t" + snpCountSampleB + "\t" + sample);
                pw.flush();

            }

            pw.close();

        } catch (FileNotFoundException e) {
            logger.error("Exception while writing sample assignment summary.", e);  //To change body of catch statement use File | Settings | File Templates.
        }

    }

    protected Map<String, Integer[]> readReadsToSampleInformation(File inputPath, String sampleA, String sampleB) {

        Map<String, Integer[]> retVal = new HashMap<String, Integer[]>();

        try {

            BufferedReader br = new BufferedReader(new FileReader(inputPath));

            String line;
            while ((line = br.readLine()) != null) {

                String[] cols = line.split("\t");
                String readName = cols[0];
                String sample = cols[5];

                if (!retVal.containsKey(readName)) {
                    retVal.put(readName, new Integer[2]);
                    retVal.get(readName)[0] = 0;
                    retVal.get(readName)[1] = 0;
                }

                if (sample.equals(sampleA)) {
                    retVal.get(readName)[0]++;
                }

                if (sample.equals(sampleB)) {
                    retVal.get(readName)[1]++;
                }

            }

            br.close();

        } catch (IOException e) {
            logger.error("Unable to open file containing read to SNP to sample information.", e);
        }

        return retVal;
    }

    protected String determineSample(VCFRecord record, String allele, String sampleA, String sampleB) {

        String retVal = null;

        if (allele.equalsIgnoreCase(record.getReferenceAllele())) {

            if (record.getSampleRecord(sampleA).isReference()) {
                retVal = sampleA;
            } else if (record.getSampleRecord(sampleB).isReference()) {
                retVal = sampleB;
            }

        } else if (allele.equalsIgnoreCase(record.getAlternativeAllele())) {

            if (record.getSampleRecord(sampleA).isAlternative()) {
                retVal = sampleA;
            } else if (record.getSampleRecord(sampleB).isAlternative()) {
                retVal = sampleB;
            }

        }

        return retVal;
    }

    /**
     * Returns the minimum minor allele frequency at which a SNP called as heterozygous
     * in the VCF file will be returned as heterozyogous.
     *
     * @param cutoff the relative frequency
     */
    @Override
    public void setHeterozygousMinorAlleleFrequencyCutoff(double cutoff) {
        this.minorAlleleFrequencyCutoff = cutoff;
    }

    /**
     * Sets the minimum minor allele frequency at which a SNP called as heterozygous
     * in the VCF file will be returned as heterozyogous.
     *
     * @return the relative frequency
     */
    @Override
    public double getHeterozygousMinorAlleleFrequencyCutoff() {
        return minorAlleleFrequencyCutoff;
    }
}
