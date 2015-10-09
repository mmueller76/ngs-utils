package org.nxtgenutils.impl;

import net.sf.picard.util.CigarUtil;
import net.sf.samtools.*;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Map;
import java.util.HashMap;
import java.util.TreeMap;

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
 * Date: 28-Sep-2011
 * Time: 15:33:16
 * To change this template use File | Settings | File Templates.
 * <p/>

 */

/**
 * Clips 3' end of overlapping read pairs. Higher quality end is retained.
 * Output BAM is unsorted.
 *
 * @author Michael Mueller
 */
public class ReadPairClipper {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(ReadPairClipper.class);

    private int inconsistantStrandedness = 0;
    private int totalClippedBases = 0;
    private int totalClippedReadPairs = 0;
    private int totalReads = 0;

    private Map<Integer, Integer> overlapDistribution = new TreeMap<Integer, Integer>();
    private File mappingInput;
    private File mappingOutput;


    public ReadPairClipper(File mapping) throws IllegalStateException {
        this(mapping, null, null);
    }

    public ReadPairClipper(File mappingInput, File mappingOutput, File statsOutput) throws IllegalStateException {

        boolean pickFirst = true;

        this.mappingInput = mappingInput;
        this.mappingOutput = mappingOutput;

        //create SAM reader & set validation stringency
        SAMFileReader samFileReader = new SAMFileReader(mappingInput);
        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMFileHeader header = samFileReader.getFileHeader();
        header.addComment("overlapping read pairs clipped retaining higher-quality end");
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        SAMFileWriter bamFileWriter;

        if (mappingOutput != null) {

            bamFileWriter = new SAMFileWriterFactory().makeBAMWriter(samFileReader.getFileHeader(),
                    true, mappingOutput);

        } else {
            bamFileWriter = new SAMFileWriterFactory().makeSAMWriter(samFileReader.getFileHeader(), true, System.out);
        }


        //iterate over input file and store reads that
        //are mapped as proper pairs to check for
        //overlap
        Map<String, SAMRecord> potentiallyOverlappingMates = new HashMap<String, SAMRecord>();

        for (SAMRecord read : samFileReader) {


            //count reads
            totalReads++;

            if (totalReads % 1000000 == 0) {
                System.out.println(totalReads + " reads processed");
            }

            String readName = read.getReadName();

            //check if mate has been seen already
            if (potentiallyOverlappingMates.containsKey(readName)) {

                //if so, get mate
                SAMRecord mate = potentiallyOverlappingMates.get(readName);

                SAMRecord firstRead = read;
                SAMRecord secondRead = mate;

                if (!read.getReadUnmappedFlag() && read.getAlignmentStart() != 0) {

                    //order reads by mappingInput location
                    if (mate.getAlignmentStart() < read.getAlignmentStart()) {
                        firstRead = mate;
                        secondRead = read;
                    }

                    //if first read in genomic order on negative strand OR
                    //second read on positive strand there is something
                    //wrong -> skip pair
                    if (firstRead.getReadNegativeStrandFlag() || !secondRead.getReadNegativeStrandFlag()) {
                        inconsistantStrandedness++;
                    } else {

                        int overlapBeforeClipping = firstRead.getAlignmentEnd() - secondRead.getAlignmentStart();

                        //if first read in genomic order on negative strand OR
                        //second read on positive strand there is something
                        //wrong -> skip pair
                        if (overlapBeforeClipping > 0) {

                            totalClippedReadPairs++;

                            if (!overlapDistribution.containsKey(overlapBeforeClipping)) {
                                overlapDistribution.put(overlapBeforeClipping, 0);
                            }

                            int count = overlapDistribution.get(overlapBeforeClipping);
                            count++;
                            overlapDistribution.put(overlapBeforeClipping, count);

                            int clippingPositionFirstRead = getReadClippingPosition(firstRead, secondRead, true);
                            int clippingPositionSecondRead = getReadClippingPosition(firstRead, secondRead, false);

                            double averageQualityClippedBasesFirstRead = getAverageBaseQualityClippedBases(firstRead, clippingPositionFirstRead);
                            double averageQualityClippedBasesSecondRead = getAverageBaseQualityClippedBases(secondRead, clippingPositionSecondRead);

                            //choose read to clip based on quality
                            //of clipped bases
                            if (averageQualityClippedBasesFirstRead > averageQualityClippedBasesSecondRead) {
                                CigarUtil.softClip3PrimeEndOfRead(secondRead, clippingPositionSecondRead);
                                totalClippedBases = totalClippedBases + secondRead.getReadLength() - clippingPositionSecondRead;
                            } else if (averageQualityClippedBasesFirstRead < averageQualityClippedBasesSecondRead) {
                                CigarUtil.softClip3PrimeEndOfRead(firstRead, clippingPositionFirstRead);
                                totalClippedBases = totalClippedBases + firstRead.getReadLength() - clippingPositionFirstRead;
                            }

                            //if quality is equal pick "random" read
                            //since read sampling is random alternating
                            //between first and second read will pick a random read
                            //and at the same time ensures reproducability of
                            //the clipping results
                            else {
                                if (pickFirst) {
                                    CigarUtil.softClip3PrimeEndOfRead(firstRead, clippingPositionFirstRead);
                                    totalClippedBases = totalClippedBases + firstRead.getReadLength() - clippingPositionFirstRead;
                                    pickFirst = false;
                                } else {
                                    CigarUtil.softClip3PrimeEndOfRead(secondRead, clippingPositionSecondRead);
                                    totalClippedBases = totalClippedBases + secondRead.getReadLength() - clippingPositionSecondRead;
                                    pickFirst = true;
                                }
                            }

                        }

                    }

                }
                //write clipped pair to mappingOutput
                bamFileWriter.addAlignment(read);
                bamFileWriter.addAlignment(mate);

                //remove read from hash map
                potentiallyOverlappingMates.remove(readName);

            } else if (read.getProperPairFlag() && !read.getReadUnmappedFlag() && read.getAlignmentStart() != 0) {
                potentiallyOverlappingMates.put(readName, read);
            } else {
                bamFileWriter.addAlignment(read);
            }

        }

        //flush potentially overlapping read buffer
        for (String readId : potentiallyOverlappingMates.keySet()) {
            SAMRecord read = potentiallyOverlappingMates.get(readId);
            bamFileWriter.addAlignment(read);
        }

        samFileReader.close();
        bamFileWriter.close();

        writeStats(statsOutput);

    }

    private int get3PrimeSoftClippedBaseCount(SAMRecord samRecord) {

        int retVal = 0;

        //if mapping on (-) strand
        try {

            if (samRecord.getReadNegativeStrandFlag()) {

                //3' clipped bases are encoded at start of Cigar string
                if (samRecord.getCigar().getCigarElement(0).getOperator() == CigarOperator.SOFT_CLIP) {
                    retVal = samRecord.getCigar().getCigarElement(0).getLength();
                }

            } else {

                //3' clipped bases are encoded at the end of Cigar string
                int cigarLength = samRecord.getCigarLength();
                int lastCigarIndex = cigarLength - 1;
                if (samRecord.getCigar().getCigarElement(lastCigarIndex).getOperator() == CigarOperator.SOFT_CLIP) {
                    retVal = samRecord.getCigar().getCigarElement(lastCigarIndex).getLength();
                }

            }
        } catch (Exception e) {
            System.out.println(e.getMessage());
            System.out.println(samRecord.getSAMString());
        }

        return retVal;

    }

    private int getReadClippingPosition(SAMRecord firstRead, SAMRecord secondRead, boolean forFirstRead) {


        //get soft clipped bases...
        int softClippedBasesFirstRead = get3PrimeSoftClippedBaseCount(firstRead);
        int softClippedBasesSecondRead = get3PrimeSoftClippedBaseCount(secondRead);

        int overlap = firstRead.getAlignmentEnd() - secondRead.getAlignmentStart();

        int retVal;

        if (forFirstRead) {

            //calculate first read clipping position
            retVal = firstRead.getReadLength() - softClippedBasesFirstRead - overlap;
            int clippingReferencePosition = firstRead.getReferencePositionAtReadPosition(retVal);

            //in case of clipped end containing insertions or deletions
            //the clipping position might have to be adjusted
            if (clippingReferencePosition > secondRead.getAlignmentStart()) {

                while (clippingReferencePosition > secondRead.getAlignmentStart() && retVal > 1 && retVal < firstRead.getReadLength()) {
                    retVal--;
                    clippingReferencePosition = firstRead.getReferencePositionAtReadPosition(retVal);
                }

            } else if (clippingReferencePosition < secondRead.getAlignmentStart()) {

                while (clippingReferencePosition < secondRead.getAlignmentStart() && retVal > 1 && retVal < firstRead.getReadLength()) {
                    retVal++;
                    clippingReferencePosition = firstRead.getReferencePositionAtReadPosition(retVal);
                }

            }

        } else {

            //calculate second read clipping position
            retVal = softClippedBasesSecondRead + overlap + 1;
            int clippingReference = secondRead.getReferencePositionAtReadPosition(retVal);

            //in case of clipped end containing insertions or deletions
            //the clipping position might have to be adjusted
            if (clippingReference > firstRead.getAlignmentEnd()) {

                while (clippingReference > firstRead.getAlignmentEnd() && retVal > 1 && retVal < secondRead.getReadLength()) {
                    retVal--;
                    clippingReference = secondRead.getReferencePositionAtReadPosition(retVal);
                }

            } else if (clippingReference < firstRead.getAlignmentEnd()) {

                while (clippingReference < firstRead.getAlignmentEnd() && retVal > 1 && retVal < secondRead.getReadLength()) {
                    retVal++;
                    clippingReference = secondRead.getReferencePositionAtReadPosition(retVal);
                }

            }

            //as read position on negative strand is given relative to reference alignment and
            //not to read direction we have to inverse the value
            retVal = secondRead.getReadLength() - retVal + 1;

        }

        return retVal;

    }

    private double getAverageBaseQualityClippedBases(SAMRecord record, int clippingPosition) {

        //get soft clipped bases...
        int softClippedBases = get3PrimeSoftClippedBaseCount(record);

        if(clippingPosition < 1 || clippingPosition > record.getReadLength() - softClippedBases){
            logger.error("Clipping position provided to getAverageBaseQualityClippedBases() not within valid range: " + clippingPosition);
            logger.error(record.getSAMString());
            return -1;
        }

        //get average quality of clipped bases
        byte[] quals = record.getBaseQualities();
        double qualsSum = 0;
        double retVal = 0;
        int bases = 0;

        if (!record.getReadNegativeStrandFlag()) {
            for (int i = clippingPosition - 1; i < record.getReadLength() - softClippedBases; i++) {

                int qual = quals[i];
                qualsSum = qualsSum + qual;

                bases++;
            }

        } else {

            for (int i = softClippedBases; i < record.getReadLength() - clippingPosition; i++) {

                int qual = quals[i];
                qualsSum = qualsSum + qual;

                bases++;
            }

        }
        retVal = qualsSum / (double) bases;
        return retVal;

    }


    private void writeStats(File statsOut) {

        try {

            PrintWriter pw = new PrintWriter(statsOut);

            pw.println("input BAM file: " + mappingInput.getAbsolutePath());
            if (mappingOutput != null) {
                pw.println("output BAM file: " + mappingOutput.getAbsolutePath());
            } else {
                pw.println("output BAM file: STDOUT");
            }

            pw.flush();

            //calculate average overlap length
            int n = 0;
            int sum = 0;

            for (Integer length : overlapDistribution.keySet()) {

                int frequency = overlapDistribution.get(length);
                n = n + frequency;
                sum = sum + length * frequency;

            }

            double averageOverlap = (double) sum / (double) n;


            pw.println(totalReads / 2 + " read pairs processed");
            pw.println(inconsistantStrandedness / 2 + " read pairs skipped because of inconsistant strandedness");
            pw.println(totalClippedReadPairs + " read pairs clipped");
            pw.println(totalClippedBases + " bases clipped");
            pw.println(averageOverlap + " average overlap");
            pw.println();
            pw.flush();

            pw.println("overlap\tcount");
            for (Integer length : overlapDistribution.keySet()) {

                int frequency = overlapDistribution.get(length);
                pw.println(length + "\t" + frequency);
                pw.flush();

            }


            pw.close();

        } catch (FileNotFoundException e) {

            logger.error("Unable to write filter stats to " + statsOut.getAbsolutePath() + ".", e);

        }


    }

    public static void main(String[] args) {

        //new ReadPairClipper(new File("/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.rn4.ct.fa.sorted.nophix.processed.head.bam"));
        //new ReadPairClipper(new File("/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.rn4.ga.fa.sampe.uncon.read.readname.sorted.bam"));
        new ReadPairClipper(new File("/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.rn4.ga.fa.sampe.uncon.read.readname.sorted.bam"));


    }

}
