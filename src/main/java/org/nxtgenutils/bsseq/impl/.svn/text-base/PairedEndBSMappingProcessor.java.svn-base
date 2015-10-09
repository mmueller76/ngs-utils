package org.nxtgenutils.bsseq.impl;

import net.sf.samtools.*;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;


import org.apache.log4j.Logger;
import org.nxtgenutils.bsseq.AbstractBSMappingProcessor;
import org.nxtgenutils.bsseq.BSTemplateStrand;
import org.nxtgenutils.bsseq.BSPairingStrandSpecificity;

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
 * for a pair of BAM files of BS-seq reads aligned to the
 * forward (C->T converted reference) and reverse strand
 * (G->A converted reference) adds BS-seq specific
 * information as attributes in the 'OPT' column:
 * <p/>
 * BA:A:<coordinates> coordinates of amplicon if BS-seq data
 * generated from PCR amplified genomic region
 * <p/>
 * BT:A:<value> BS converted template strand the mapping record
 * refers to (character attribute)
 * BT:A:F BS-converted forward strand
 * BT:A:R BS-converted reverse strand
 * <p/>
 * BS:A:<value> strand specificity of mapping (character attribute)
 * BS:A:F read maps to BS-converted forward strand only
 * BS:A:R read maps to BS-converted reverse strand only
 * BS:A:B read maps to both BS-converted strands
 * BS:A:N read does not map to either strand
 * <p/>
 * BP:A:<value> strand with proper read pairing (character attribute)
 * BP:A:F read properly paired on BS-converted forward strand only
 * BP:A:R read properly paired on BS-converted reverse strand only
 * BP:A:B read properly paired on both BS-converted strands
 * BP:A:N read not properly paired on either strand
 * <p/>
 * BW:i:<value>
 * Phred scaled probablitity of read mapped to incorrect strand
 * based on C residue frequency (integer attribute)
 * <p/>
 * BV:i:<frequency>
 * number of valid BS-mismatches to the reference, i.e. T in read and C in reference
 * (forward strand: unconverted reference = C && unconverted read = T,
 * reverse strand: unconverted reference = G && unconverted read = A)
 * <p/>
 * BI:i:<frequency>
 * number of BS-invalid mismatches to the reference, i.e. C in read and T in reference
 * i.e. due to sequencing errors/SNPs
 * <p/>
 * BN:i:<frequency>
 * number of non-BS mismatches to the reference, i.e. any mismatch other than T/C mismatches on forward strand and G/A mismatches on reverse strand
 * because of sequencing errors/SNPs
 * <p/>
 * BC:Z:<frequency>
 * BS conversion rate outside CpG context = ratio(number of unconverted Cs/total number of Cs)
 * <p/>
 * BM:Z:<frequency>
 * BS conversion rate in CpG context = ratio(number of unconverted Cs/total number of Cs)
 * BO:Z:<value> BS-mapping orientation
 * BO:Z:V valid:
 * - forward strand template: read1 mapping to (+) strand and read2 mapping to (-) strand
 * - reverse strand template: read1 mapping to (-) strand and read2 mapping to (+) strand
 * BO:Z:I invalid
 * <p/>
 * BWA custom tags
 * <p/>
 * NM:i	Edit distance
 * MD:Z	Mismatching positions/bases
 * AS:i	Alignment score
 * X0:i	Number of best hits
 * X1:i	Number of suboptimal hits found by BWA
 * XN:i	Number of ambiguous bases in the referenece
 * XM:i	Number of mismatches in the alignment
 * XO:i	Number of gap opens
 * XG:i	Number of gap extentions
 * XT:A	Type: Unique/Repeat/N/Mate-sw
 * XA:Z	Alternative hits; format: (chr,pos,CIGAR,NM;)*
 * XS:i	Suboptimal alignment score
 * XF	Support from forward/reverse alignment
 * XE	Number of supporting seeds
 * <p/>
 * <p/>
 * A summary statistic of the processing is printed to the
 * STDOUT
 * <p/>
 * <p/>
 * <p/>
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 26-Aug-2010
 * Time: 17:04:35
 * To change this template use File | Settings | File Templates.
 */
public class PairedEndBSMappingProcessor extends AbstractBSMappingProcessor {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(PairedEndBSMappingProcessor.class);

    public PairedEndBSMappingProcessor(File mappingForwardStrand, File mappingReverseStrand, File referenceSequenceFile) throws IllegalStateException {
        this(mappingForwardStrand, mappingReverseStrand, referenceSequenceFile, null);
    }

    /**
     * @param mappingForwardStrand
     * @param mappingReverseStrand
     * @param referenceSequenceFile
     * @param outputDirectory
     * @throws IllegalStateException
     */
    public PairedEndBSMappingProcessor(File mappingForwardStrand, File mappingReverseStrand, File referenceSequenceFile, File outputDirectory) throws IllegalStateException {

        if (outputDirectory == null) {
            outputDirectory = new File(mappingForwardStrand.getAbsolutePath().replace(mappingForwardStrand.getName(), ""));
        }

        File mappingForwardStrandOutput = new File(outputDirectory.getAbsolutePath() + "/" + mappingForwardStrand.getName().replace(".bam", ".processed.bam"));
        File mappingReverseStrandOutput = new File(outputDirectory.getAbsolutePath() + "/" + mappingReverseStrand.getName().replace(".bam", ".processed.bam"));
        File statsOutput = new File(outputDirectory.getAbsolutePath() + "/" + mappingForwardStrand.getName().replace(".bam", ".processed.bam.stats"));

        logger.info("post-processing paired end bisulphite-seq mapping");
        logger.info("forward strand mapping input           :" + mappingForwardStrand.getAbsolutePath());
        logger.info("reverse strand mapping input           :" + mappingReverseStrand.getAbsolutePath());
        logger.info("reference sequence                     :" + referenceSequenceFile.getAbsolutePath());
        logger.info("processed forward strand mapping output:" + mappingForwardStrandOutput.getAbsolutePath());
        logger.info("processed reverse strand mapping output:" + mappingReverseStrandOutput.getAbsolutePath());
        logger.info("processing statistics output:" + mappingReverseStrandOutput.getAbsolutePath());

        if (referenceSequenceFile.getName().endsWith(".gz")) {
            logger.error("Compressed reference sequence not supported. Decompress first.");
            System.exit(1);
        }

        //create reference sequence file
        loadReferenceSequences(ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceSequenceFile));

        //create SAM reader for forward and reverse strand
        SAMFileReader samFileReaderForward = new SAMFileReader(mappingForwardStrand);
        SAMFileReader samFileReaderReverse = new SAMFileReader(mappingReverseStrand);
        
        SAMFileWriter bamFileWriterForward = new SAMFileWriterFactory().makeBAMWriter(samFileReaderForward.getFileHeader(),
                true, mappingForwardStrandOutput);
        SAMFileWriter bamFileWriterReverse = new SAMFileWriterFactory().makeBAMWriter(samFileReaderReverse.getFileHeader(),
                true, mappingReverseStrandOutput);

        //set validation stringency
        samFileReaderForward.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        samFileReaderReverse.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        //create SAM record iterator
        //if SAM files were not sorted with Picard tools but samtools
        //validator will throw exceptions because of wrong sort order
//        SAMRecordIterator iteratorForwardStrand = samFileReaderForward.iterator().assertSorted(SAMFileHeader.SortOrder.queryname);
//        SAMRecordIterator iteratorReverseStrand = samFileReaderReverse.iterator().assertSorted(SAMFileHeader.SortOrder.queryname);
        SAMRecordIterator iteratorForwardStrand = samFileReaderForward.iterator();
        SAMRecordIterator iteratorReverseStrand = samFileReaderReverse.iterator();


        //iterate synchronously over entries in
        //forward and reverse strand mapping
        while (iteratorForwardStrand.hasNext() && iteratorReverseStrand.hasNext()) {

            //read next forward and reverse strand
            // read pair mapping
            SAMRecord read1ForwardStrandMapping = null;
            SAMRecord read1ReverseStrandMapping = null;

            SAMRecord read2ForwardStrandMapping = null;
            SAMRecord read2ReverseStrandMapping = null;

            //count read pair
            readPairCount++;

            //count read
            readCount++;
            try {
                read1ForwardStrandMapping = iteratorForwardStrand.next();
            } catch (Exception e) {

                logger.error("exception while reading forward strand mapping for read1 of read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }


            try {
                read1ReverseStrandMapping = iteratorReverseStrand.next();
            } catch (Exception e) {

                logger.error("exception while reading reverse strand mapping for read1 of read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }

            //count record
            readCount++;
            try {
                read2ForwardStrandMapping = iteratorForwardStrand.next();
            } catch (Exception e) {

                logger.error("exception while reading forward strand mapping for read2 of read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }

            try {
                read2ReverseStrandMapping = iteratorReverseStrand.next();
            } catch (Exception e) {

                logger.error("exception while reading reverse strand mapping for read2 of read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }


            //make sure forward and reverse strand
            //mapping are not null
            if (read1ForwardStrandMapping == null) {
                logger.error("read 1 forward strand mapping is NULL");
                System.exit(1);
            }
            if (read1ReverseStrandMapping == null) {
                logger.error("read1 reverse strand mapping is NULL");
                System.exit(1);
            }
            if (read2ForwardStrandMapping == null) {
                logger.error("read2 forward strand mapping is NULL");
                System.exit(1);
            }
            if (read2ReverseStrandMapping == null) {
                logger.error("read2 reverse strand mapping is NULL");
                System.exit(1);
            }

            //make sure read1 and read2 are in correct order
            if (!read1ForwardStrandMapping.getFirstOfPairFlag()) {
                SAMRecord readTemp = read1ForwardStrandMapping;
                read1ForwardStrandMapping = read2ForwardStrandMapping;
                read2ForwardStrandMapping = readTemp;
            }

            if (!read1ReverseStrandMapping.getFirstOfPairFlag()) {
                SAMRecord readTemp = read1ReverseStrandMapping;
                read1ReverseStrandMapping = read2ReverseStrandMapping;
                read2ReverseStrandMapping = readTemp;
            }

            //make sure forward and reverse strand mapping are for the same read
            if (!readIsEqual(read1ForwardStrandMapping, read1ReverseStrandMapping)) {

                String readForwardString = "";
                String readReverseString = "";

                if (read1ForwardStrandMapping.getReadNegativeStrandFlag()) {
                    readForwardString = getReverseComplementReadString(read1ForwardStrandMapping);
                } else {
                    readForwardString = read1ForwardStrandMapping.getReadString();
                }

                if (read1ReverseStrandMapping.getReadNegativeStrandFlag()) {
                    readReverseString = getReverseComplementReadString(read1ForwardStrandMapping);
                } else {
                    readReverseString = read1ReverseStrandMapping.getReadString();
                }

                logger.error(readCount + " reads NOT equal:");
                logger.error(read1ForwardStrandMapping.getReadName() + " " + readForwardString);
                logger.error(read1ReverseStrandMapping.getReadName() + " " + readReverseString);
                logger.error("Make sure both SAM files are sorted by read name");
                System.exit(1);

            }

            //fix mapping flag/CIGAR string
            //some records (in read pairs where one read maps and the other one does not)
            //will have the read unmapped flag set to true but still have a CIGAR string
            fixRecord(read1ForwardStrandMapping);
            fixRecord(read2ForwardStrandMapping);
            fixRecord(read1ReverseStrandMapping);
            fixRecord(read2ReverseStrandMapping);

            //append BS template strand tag
            appendBSTemplateStrandTag(read1ForwardStrandMapping, BSTemplateStrand.FORWARD);
            appendBSTemplateStrandTag(read2ForwardStrandMapping, BSTemplateStrand.FORWARD);
            appendBSTemplateStrandTag(read1ReverseStrandMapping, BSTemplateStrand.REVERSE);
            appendBSTemplateStrandTag(read2ReverseStrandMapping, BSTemplateStrand.REVERSE);

            //append BS mapping strand specificity tag
            appendMappingStrandSpecificityTag(read1ForwardStrandMapping, read1ReverseStrandMapping);
            appendMappingStrandSpecificityTag(read2ForwardStrandMapping, read2ReverseStrandMapping);

            //append BS proper pairing strand tag
            appendProperPairingStrandTag(read1ForwardStrandMapping,
                    read2ForwardStrandMapping,
                    read1ReverseStrandMapping,
                    read2ReverseStrandMapping);

            //append BS mapping orientation tag
            appendBSMappingOrientationTag(read1ForwardStrandMapping);
            appendBSMappingOrientationTag(read2ForwardStrandMapping);
            appendBSMappingOrientationTag(read1ReverseStrandMapping);
            appendBSMappingOrientationTag(read2ReverseStrandMapping);

            //append mismatch tags (non-BS mismatches, BS-valid/BS-invalid mismatches)
            appendMismatchTags(read1ForwardStrandMapping);
            appendMismatchTags(read2ForwardStrandMapping);
            appendMismatchTags(read1ReverseStrandMapping);
            appendMismatchTags(read2ReverseStrandMapping);

            //append BS conversion tag (inside and outside CpG context)
            appendBSConversionTags(read1ForwardStrandMapping);
            appendBSConversionTags(read2ForwardStrandMapping);
            appendBSConversionTags(read1ReverseStrandMapping);
            appendBSConversionTags(read2ReverseStrandMapping);

            //unless both reads are not mapped to forward strand template...
            if (!(read1ForwardStrandMapping.getReadUnmappedFlag() && read2ForwardStrandMapping.getReadUnmappedFlag())) {
                //..write to output
                bamFileWriterForward.addAlignment(read1ForwardStrandMapping);
                bamFileWriterForward.addAlignment(read2ForwardStrandMapping);
            }

            //unless both reads are not mapped to reverse strand template...
            if (!(read1ReverseStrandMapping.getReadUnmappedFlag() && read2ReverseStrandMapping.getReadUnmappedFlag())) {
                //..write to output
                bamFileWriterReverse.addAlignment(read1ReverseStrandMapping);
                bamFileWriterReverse.addAlignment(read2ReverseStrandMapping);
            }

            if (readCount > 0 && readCount % 1000000 == 0) {
                logger.info(readCount + " records processed");
//                if (readCount == 1000000) {
//                    break;
//                }
            }


        }

        samFileReaderForward.close();
        samFileReaderReverse.close();

        bamFileWriterForward.close();
        bamFileWriterReverse.close();

        printStatistics(statsOutput);

        logger.info("done");

    }


    /**
     * @param outputPath
     */
    protected void printStatistics(File outputPath) {

        try {

            PrintWriter pw = new PrintWriter(outputPath);

            pw.println(readCount + " reads processed");
            pw.println("---------------------------------------");
            pw.println("BS mapping strand");
            pw.println(unmapped + "(" + oneDecimalPlace.format((double) unmapped / (double) readCount * 100) + "%) unmapped reads");
            pw.println(mappingForwardStrandOnly + "(" + oneDecimalPlace.format((double) mappingForwardStrandOnly / (double) readCount * 100) + "%) reads mapping to forward strand only");
            pw.println(mappingReverseStrandOnly + "(" + oneDecimalPlace.format((double) mappingReverseStrandOnly / (double) readCount * 100) + "%) reads mapping to reverse strand only");
            pw.println(mapingBothStrands + "(" + oneDecimalPlace.format((double) mapingBothStrands / (double) readCount * 100) + "%) reads mapping to both strands");
            //pw.println((unmapped + mappingForwardStrandOnly + mappingReverseStrandOnly + mapingBothStrands));
            pw.println(invalidMappingOrientation + "(" + oneDecimalPlace.format((double) invalidMappingOrientation / (double) (readCount - unmapped) * 100) + "%) invalid mapping orientation");
            pw.println("---------------------------------------");
            pw.println("BS proper pairing strand mapped reads");
            pw.println((notProperlyPaired - unmapped) + "(" + oneDecimalPlace.format((double) (notProperlyPaired - unmapped) / (double) (readCount - unmapped) * 100) + "%) improperly paired reads");
            pw.println(properlyPairedForwardStrandOnly + "(" + oneDecimalPlace.format((double) properlyPairedForwardStrandOnly / (double) (readCount - unmapped) * 100) + "%) properly paired forward strand only");
            pw.println(properlyPairedReverseStrandOnly + "(" + oneDecimalPlace.format((double) properlyPairedReverseStrandOnly / (double) (readCount - unmapped) * 100) + "%) properly paired reverse strand only");
            pw.println(properlyPairedBothStrands + "(" + oneDecimalPlace.format((double) properlyPairedBothStrands / (double) (readCount - unmapped) * 100) + "%) properly paired on both strands");
            pw.println(properlyPairedBothStrandsAndRepeat + "(" + oneDecimalPlace.format((double) properlyPairedBothStrandsAndRepeat / (double) (readCount - unmapped) * 100) + "%) properly paired on both strands and repeat read pair");
            pw.println((notProperlyPaired + properlyPairedForwardStrandOnly + properlyPairedReverseStrandOnly + properlyPairedBothStrands));

            pw.println("---------------------------------------");
            pw.println("non-BS mismatches");
            pw.println("count\tfreq. forward\tfreq. reverse");

            int maxNonBSMismatchForwardStrand = Collections.max(nonBSMismatchDistributionForwardStrand.keySet());
            int maxNonBSMismatchReverseStrand = Collections.max(nonBSMismatchDistributionReverseStrand.keySet());
            int maxNonBSMismatch = Math.max(maxNonBSMismatchForwardStrand, maxNonBSMismatchReverseStrand);
            int totalForwardStrand = 0;
            int totalReverseStrand = 0;
            for (int i = 0; i <= maxNonBSMismatch; i++) {

                int freqForwardStrand = 0;
                if (nonBSMismatchDistributionForwardStrand.containsKey(i)) {
                    freqForwardStrand = nonBSMismatchDistributionForwardStrand.get(i);
                }
                int freqReverseStrand = 0;
                if (nonBSMismatchDistributionReverseStrand.containsKey(i)) {
                    freqReverseStrand = nonBSMismatchDistributionReverseStrand.get(i);
                }

                pw.println(i + "\t" +
                        freqForwardStrand + "\t" +
                        "(" + oneDecimalPlace.format((double) freqForwardStrand / (double) (readCount - unmapped) * 100) + "%)" + "\t" +
                        freqReverseStrand + "\t" +
                        "(" + oneDecimalPlace.format((double) freqReverseStrand / (double) (readCount - unmapped) * 100) + "%)");

                totalForwardStrand = totalForwardStrand + freqForwardStrand;
                totalReverseStrand = totalReverseStrand + freqReverseStrand;

            }

            pw.println("---------------------------------------");
            pw.println("BS invalid mismatches");
            pw.println("count\tfreq. forward\tfreq. reverse");
            int maxBsInvalidMismatchForwardStrand = Collections.max(bsInvalidMismatchDistributionForwardStrand.keySet());
            int maxBsInvalidMismatchReverseStrand = Collections.max(bsInvalidMismatchDistributionReverseStrand.keySet());
            int maxBsInvalidMismatch = Math.max(maxBsInvalidMismatchForwardStrand, maxBsInvalidMismatchReverseStrand);
            totalForwardStrand = 0;
            totalReverseStrand = 0;
            for (int i = 0; i <= maxBsInvalidMismatch; i++) {

                int freqForwardStrand = 0;
                if (bsInvalidMismatchDistributionForwardStrand.containsKey(i)) {
                    freqForwardStrand = bsInvalidMismatchDistributionForwardStrand.get(i);
                }
                int freqReverseStrand = 0;
                if (bsInvalidMismatchDistributionReverseStrand.containsKey(i)) {
                    freqReverseStrand = bsInvalidMismatchDistributionReverseStrand.get(i);
                }

                pw.println(i + "\t" +
                        freqForwardStrand + "\t" +
                        "(" + oneDecimalPlace.format((double) freqForwardStrand / (double) (readCount - unmapped) * 100) + "%)" + "\t" +
                        freqReverseStrand + "\t" +
                        "(" + oneDecimalPlace.format((double) freqReverseStrand / (double) (readCount - unmapped) * 100) + "%)");

                totalForwardStrand = totalForwardStrand + freqForwardStrand;
                totalReverseStrand = totalReverseStrand + freqReverseStrand;

            }



            pw.println("---------------------------------------");
            pw.println("BS conversion rate");
            pw.println("    \tforward\tforward\tforward\tforward\treverse\treverse\treverse\treverse");
            pw.println("rate\tread1  \tread2  \tread2  \tread1  \tread1  \tread2  \tread2  \tread1");
            pw.println("    \t(+)    \t(-)    \t(+)    \t(-)    \t(+)    \t(-)    \t(+)    \t(-)");

            Set<String> keys = new LinkedHashSet<String>();
            keys.add(BSTemplateStrand.FORWARD + " " + false + " " + false);
            keys.add(BSTemplateStrand.FORWARD + " " + true + " " + true);
            keys.add(BSTemplateStrand.FORWARD + " " + true + " " + false);
            keys.add(BSTemplateStrand.FORWARD + " " + false + " " + true);
            keys.add(BSTemplateStrand.REVERSE + " " + false + " " + false);
            keys.add(BSTemplateStrand.REVERSE + " " + true + " " + true);
            keys.add(BSTemplateStrand.REVERSE + " " + true + " " + false);
            keys.add(BSTemplateStrand.REVERSE + " " + false + " " + true);

            for (int i = 0; i <= 100; i++) {

                StringBuffer line = new StringBuffer();

                line.append(i);

                for (String key : keys) {
                    int freq = 0;
                    if (bsConversionRate.containsKey(key) && bsConversionRate.get(key).containsKey(i)) {
                        freq = bsConversionRate.get(key).get(i);

                    }
                    line.append("\t");
                    line.append(freq);
                }

                pw.println(line.toString());


            }

            pw.flush();
            pw.close();

        } catch (Exception e) {
            logger.error("Error while writing processing stats.");
            logger.error(e.getMessage());
        }

    }

    protected void appendBSMappingOrientationTag(SAMRecord read) {

        boolean reverseStrandTemplate = false;
        if (read.getAttribute("BT").toString().equals("R")) {
            reverseStrandTemplate = true;
        }

        if (!read.getReadUnmappedFlag()) {
            if (!reverseStrandTemplate) {

                if (read.getFirstOfPairFlag() && !read.getReadNegativeStrandFlag() ||
                        !read.getFirstOfPairFlag() && read.getReadNegativeStrandFlag()) {
                    read.setAttribute("BO", 'V');
                } else {
                    read.setAttribute("BO", 'I');
                    invalidMappingOrientation++;
                }

            } else {

                if (read.getFirstOfPairFlag() && read.getReadNegativeStrandFlag() ||
                        !read.getFirstOfPairFlag() && !read.getReadNegativeStrandFlag()) {
                    read.setAttribute("BO", 'V');
                } else {
                    read.setAttribute("BO", 'I');
                    invalidMappingOrientation++;
                }

            }
        }
    }

    private BSPairingStrandSpecificity appendProperPairingStrandTag(SAMRecord read1ForwardStrandMapping,
                                                               SAMRecord read2ForwardStrandMapping,
                                                               SAMRecord read1ReverseStrandMapping,
                                                               SAMRecord read2ReverseStrandMapping) {

        //check if read maps only to one strand
        BSPairingStrandSpecificity retVal;

        boolean properlyPairedForwardStrand = read1ForwardStrandMapping.getProperPairFlag() && read2ForwardStrandMapping.getProperPairFlag();
        boolean properlyPairedReverseStrand = read1ReverseStrandMapping.getProperPairFlag() && read2ReverseStrandMapping.getProperPairFlag();

        if (properlyPairedForwardStrand && !properlyPairedReverseStrand) {

            //assign return value
            retVal = BSPairingStrandSpecificity.FORWARD;

            //set read attribute
            read1ForwardStrandMapping.setAttribute("BP", 'F');
            read2ForwardStrandMapping.setAttribute("BP", 'F');
            read1ReverseStrandMapping.setAttribute("BP", 'F');
            read2ReverseStrandMapping.setAttribute("BP", 'F');

            //increment count
            properlyPairedForwardStrandOnly = properlyPairedForwardStrandOnly + 2;

        } else if (!properlyPairedForwardStrand && properlyPairedReverseStrand) {

            //assign return value
            retVal = BSPairingStrandSpecificity.REVERSE;

            //set read attribute
            read1ForwardStrandMapping.setAttribute("BP", 'R');
            read2ForwardStrandMapping.setAttribute("BP", 'R');
            read1ReverseStrandMapping.setAttribute("BP", 'R');
            read2ReverseStrandMapping.setAttribute("BP", 'R');

            //increment count
            properlyPairedReverseStrandOnly = properlyPairedReverseStrandOnly + 2;

        } else if (properlyPairedForwardStrand && properlyPairedReverseStrand) {

            //assign return value
            retVal = BSPairingStrandSpecificity.BOTH;

            //set read attribute
            read1ForwardStrandMapping.setAttribute("BP", 'B');
            read2ForwardStrandMapping.setAttribute("BP", 'B');
            read1ReverseStrandMapping.setAttribute("BP", 'B');
            read2ReverseStrandMapping.setAttribute("BP", 'B');

            //increment count
            properlyPairedBothStrands = properlyPairedBothStrands + 2;

        } else {

            //assign return value
            retVal = BSPairingStrandSpecificity.NONE;

            //set read attribute
            read1ForwardStrandMapping.setAttribute("BP", 'N');
            read2ForwardStrandMapping.setAttribute("BP", 'N');
            read1ReverseStrandMapping.setAttribute("BP", 'N');
            read2ReverseStrandMapping.setAttribute("BP", 'N');

            //increment count
            notProperlyPaired = notProperlyPaired + 2;

        }

        //check if read has proper pairing on both strand
        //AND is repeat read
        if (retVal == BSPairingStrandSpecificity.BOTH) {

            String read1TypeForwardStrand = read1ForwardStrandMapping.getCharacterAttribute("XT").toString();
            String read2TypeForwardStrand = read2ForwardStrandMapping.getCharacterAttribute("XT").toString();
            String read1TypeReverseStrand = read1ReverseStrandMapping.getCharacterAttribute("XT").toString();
            String read2TypeReverseStrand = read2ReverseStrandMapping.getCharacterAttribute("XT").toString();

            if ((read1TypeForwardStrand.equals("R") || read2TypeForwardStrand.equals("R")) ||
                    (read1TypeReverseStrand.equals("R") || read2TypeReverseStrand.equals("R"))) {
                properlyPairedBothStrandsAndRepeat = properlyPairedBothStrandsAndRepeat + 2;
            }

        }

        return retVal;

    }

    /**
     *
     */
    protected void appendBSConversionTags(SAMRecord read) {

        if (!read.getReadUnmappedFlag()) {

            String referenceSequenceName = read.getReferenceName().replace(".ct", "").replace(".ga", "");
            ReferenceSequence rs = referenceSequence.get(referenceSequenceName);

            int[] bsConversionInCpGContext = getBSConversion(read, rs.getBases(), 0, true);
            int[] bsConversionOutsideCpGContext = getBSConversion(read, rs.getBases(), 0, false);

            read.setAttribute("BC", bsConversionOutsideCpGContext[0] + "/" + (bsConversionOutsideCpGContext[0] + bsConversionOutsideCpGContext[1]));
            read.setAttribute("BM", bsConversionInCpGContext[0] + "/" + (bsConversionInCpGContext[0] + bsConversionInCpGContext[1]));

            if ((bsConversionOutsideCpGContext[0] + bsConversionOutsideCpGContext[1]) > 0) {

                int unconvertedCs = bsConversionOutsideCpGContext[0];
                int convertedCs = bsConversionOutsideCpGContext[1];

                int conversionRate = (int) Math.round((double) convertedCs / (double) (convertedCs + unconvertedCs) * 100);
//                System.out.println(read.format());
//                System.out.println(bsConversionOutsideCpGContext[0] + "|" + bsConversionOutsideCpGContext[1]);
//                System.out.println(conversionRate);

                BSTemplateStrand bsTemplateStrand = BSTemplateStrand.FORWARD;
                if (read.getAttribute("BT").toString().equals("R")) {
                    bsTemplateStrand = BSTemplateStrand.REVERSE;
                }


                String key = bsTemplateStrand + " " + read.getSecondOfPairFlag() + " " + read.getReadNegativeStrandFlag();
                if (!bsConversionRate.containsKey(key)) {
                    bsConversionRate.put(key, new TreeMap<Integer, Integer>());
                }
                if (!bsConversionRate.get(key).containsKey(conversionRate)) {
                    bsConversionRate.get(key).put(conversionRate, 0);
                }
                int freq = bsConversionRate.get(key).get(conversionRate);
                freq++;
                bsConversionRate.get(key).put(conversionRate, freq);


            }

        }

    }

}
