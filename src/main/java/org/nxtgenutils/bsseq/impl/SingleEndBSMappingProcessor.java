package org.nxtgenutils.bsseq.impl;

import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Set;
import java.util.LinkedHashSet;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.nxtgenutils.bsseq.AbstractBSMappingProcessor;
import org.nxtgenutils.bsseq.BSTemplateStrand;

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
 * Date: 18-Feb-2011
 * Time: 09:58:51
 */

/**
 * @author Michael Mueller
 */
public class SingleEndBSMappingProcessor extends AbstractBSMappingProcessor {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(PairedEndBSMappingProcessor.class);

    
    public SingleEndBSMappingProcessor(File mappingForwardStrand, File mappingReverseStrand, File referenceSequenceFile) throws IllegalStateException {
        this(mappingForwardStrand, mappingReverseStrand, referenceSequenceFile, null);
    }

    /*
     * @param mappingForwardStrand
     * @param mappingReverseStrand
     * @param referenceSequenceFile
     * @param outputDirectory
     * @throws IllegalStateException
     */
    public SingleEndBSMappingProcessor(File mappingForwardStrand, File mappingReverseStrand, File referenceSequenceFile, File outputDirectory) throws IllegalStateException {

        if (outputDirectory == null) {
            outputDirectory = new File(mappingForwardStrand.getAbsolutePath().replace(mappingForwardStrand.getName(), ""));
        }

        File mappingForwardStrandOutput = new File(outputDirectory.getAbsolutePath() + "/" + mappingForwardStrand.getName().replace(".bam", ".processed.bam"));
        File mappingReverseStrandOutput = new File(outputDirectory.getAbsolutePath() + "/" + mappingReverseStrand.getName().replace(".bam", ".processed.bam"));
        File statsOutput = new File(outputDirectory.getAbsolutePath() + "/" + mappingForwardStrand.getName().replace(".bam", ".processed.bam.stats"));

        logger.info("post-processing single end bisulphite-seq mapping");
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

            //if(read1ForwardStrandMapping.getReferenceName().startsWith("chr1.") && read1ReverseStrandMapping.getReferenceName().startsWith("chr1.") ){

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
            fixRecord(read1ReverseStrandMapping);


            //append BS template strand tag
            appendBSTemplateStrandTag(read1ForwardStrandMapping, BSTemplateStrand.FORWARD);
            appendBSTemplateStrandTag(read1ReverseStrandMapping, BSTemplateStrand.REVERSE);

            //append BS mapping strand specificity tag
            appendMappingStrandSpecificityTag(read1ForwardStrandMapping, read1ReverseStrandMapping);

            //append BS mapping orientation tag
            appendBSMappingOrientationTag(read1ForwardStrandMapping);
            appendBSMappingOrientationTag(read1ReverseStrandMapping);

            //append mismatch tags (non-BS mismatches, BS-valid/BS-invalid mismatches)
            appendMismatchTags(read1ForwardStrandMapping);
            appendMismatchTags(read1ReverseStrandMapping);

            //append BS conversion tag (inside and outside CpG context)
            appendBSConversionTags(read1ForwardStrandMapping);
            appendBSConversionTags(read1ReverseStrandMapping);

            //unless read is not mapped to forward strand template...
            if (!read1ForwardStrandMapping.getReadUnmappedFlag()) {

                //..write to output
                bamFileWriterForward.addAlignment(read1ForwardStrandMapping);

            }

            //unless both reads are not mapped to reverse strand template...
            if (!read1ReverseStrandMapping.getReadUnmappedFlag()) {
                //..write to output
                bamFileWriterReverse.addAlignment(read1ReverseStrandMapping);

            }

            if (readCount > 0 && readCount % 1000000 == 0) {
                logger.info(readCount + " records processed");
//                if (readCount == 1000000) {
//                    break;
//                }
            }

            }

        //}

        samFileReaderForward.close();
        samFileReaderReverse.close();

        bamFileWriterForward.close();
        bamFileWriterReverse.close();

        printStatistics(statsOutput);

        logger.info("done");

    }


    @Override
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
            pw.println("    \tforward\tforward\treverse\treverse");
            pw.println("    \t(+)    \t(-)    \t(+)    \t(-)    ");

            Set<String> keys = new LinkedHashSet<String>();
            keys.add(BSTemplateStrand.FORWARD + " " + false);
            keys.add(BSTemplateStrand.FORWARD + " " + true);
            keys.add(BSTemplateStrand.REVERSE + " " + false);
            keys.add(BSTemplateStrand.REVERSE + " " + true);


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

    @Override
    protected void appendBSMappingOrientationTag(SAMRecord read) {

        boolean reverseStrandTemplate = false;
        if (read.getAttribute("BT").toString().equals("R")) {
            reverseStrandTemplate = true;
        }

        if (!read.getReadUnmappedFlag()) {

            if (!reverseStrandTemplate) {

                if (!read.getReadNegativeStrandFlag()) {
                    read.setAttribute("BO", 'V');
                } else {
                    read.setAttribute("BO", 'I');
                    invalidMappingOrientation++;
                }

            } else {

                if (read.getReadNegativeStrandFlag()) {
                    read.setAttribute("BO", 'V');
                } else {
                    read.setAttribute("BO", 'I');
                    invalidMappingOrientation++;
                }

            }
        }
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


                String key = bsTemplateStrand + " " + read.getReadNegativeStrandFlag();
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
