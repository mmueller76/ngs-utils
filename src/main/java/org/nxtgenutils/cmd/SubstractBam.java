package org.nxtgenutils.cmd;

import net.sf.samtools.*;

import java.io.File;
import java.text.DecimalFormat;

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
 * Outputs a BAM file of all read pairs in BAM file A which
 * are flagged as unmapped in BAM file B. Input files have
 * to be read name sorted.
 * <p/>
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 20-Dec-2010
 * Time: 16:14:45
 */

/**
 * @author Michael Mueller
 */
public class SubstractBam {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(SubstractBam.class);

    private int readCount = 0;
    private int readPairCount = 0;
    private int readPairFilteredCount = 0;
    private int readMappingAOnly = 0;
    private int readMappingBOnly = 0;
    private int readMappingBoth = 0;
    private int readMappingNone = 0;

    private static DecimalFormat oneDecimalPlace = new DecimalFormat("#.#");

    /**
     * @param bamFileA read name sorted input BAM A
     * @param bamFileB read name sorted input BAM B
     */
    public SubstractBam(File bamFileA, File bamFileB, File bamFileOut) {

        logger.info("substracting BAM file B from A");
        logger.info("input BAM file A: " + bamFileA.getAbsolutePath());
        logger.info("input BAM file B: " + bamFileB.getAbsolutePath());
        logger.info("output BAM file : " + bamFileOut.getAbsolutePath());
        logger.info("the output BAM file will contain all read pairs in ");
        logger.info("input BAM file A which are not mapped in input BAM file B");

        //create SAM reader for forward and reverse strand
        SAMFileReader samFileReaderA = new SAMFileReader(bamFileA);
        SAMFileReader samFileReaderB = new SAMFileReader(bamFileB);

        SAMFileWriter bamWriter = new SAMFileWriterFactory().makeBAMWriter(samFileReaderA.getFileHeader(),
                true, bamFileOut);

        //set validation stringency
        samFileReaderA.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        samFileReaderB.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        //create SAM record iterator
        //if SAM files were not sorted with Picard tools but samtools
        //validator will throw exceptions because of wrong sort order
//        SAMRecordIterator iteratorForwardStrand = samFileReaderForward.iterator().assertSorted(SAMFileHeader.SortOrder.queryname);
//        SAMRecordIterator iteratorReverseStrand = samFileReaderReverse.iterator().assertSorted(SAMFileHeader.SortOrder.queryname);
        SAMRecordIterator iteratorA = samFileReaderA.iterator();
        SAMRecordIterator iteratorB = samFileReaderB.iterator();

//iterate synchronously over entries in
        //forward and reverse strand mapping
        while (iteratorA.hasNext() && iteratorB.hasNext()) {

            //read next forward and reverse strand
            // read pair mapping
            SAMRecord read1A = null;
            SAMRecord read1B = null;

            SAMRecord read2A = null;
            SAMRecord read2B = null;

            //count read pair
            readPairCount++;

            //count read
            readCount++;
            try {
                read1A = iteratorA.next();
            } catch (Exception e) {

                logger.error("exception while reading read1 in input file A read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }


            try {
                read1B = iteratorB.next();
            } catch (Exception e) {

                logger.error("exception while reading read1 in input file B read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }

            //count record
            readCount++;
            try {
                read2A = iteratorA.next();
            } catch (Exception e) {

                logger.error("exception while reading read2 in input file A read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }

            try {
                read2B = iteratorB.next();
            } catch (Exception e) {

                logger.error("exception while reading read2 in input file B read pair " + readPairCount);
                logger.error(e.getMessage());
                System.exit(1);

            }


            //make sure forward and reverse strand
            //mapping are not null
            if (read1A == null) {
                logger.error("read 1 input file A is NULL");
                System.exit(1);
            }
            if (read1B == null) {
                logger.error("read1 input file B is NULL");
                System.exit(1);
            }
            if (read2A == null) {
                logger.error("read2 input file A is NULL");
                System.exit(1);
            }
            if (read2B == null) {
                logger.error("read2 input file B is NULL");
                System.exit(1);
            }

            if(!read1A.getReadUnmappedFlag()&& read1B.getReadUnmappedFlag()){
                readMappingAOnly++;
            } else if(read1A.getReadUnmappedFlag()&& !read1B.getReadUnmappedFlag()){
                readMappingBOnly++;
            } else if (!read1A.getReadUnmappedFlag()&& !read1B.getReadUnmappedFlag()){
                readMappingBoth++;
            }else if(read1A.getReadUnmappedFlag() && read1B.getReadUnmappedFlag()){
                readMappingNone++;
            }

            if(!read2A.getReadUnmappedFlag()&& read2B.getReadUnmappedFlag()){
                readMappingAOnly++;
            } else if(read2A.getReadUnmappedFlag()&& !read2B.getReadUnmappedFlag()){
                readMappingBOnly++;
            } else if (!read2A.getReadUnmappedFlag()&& !read2B.getReadUnmappedFlag()){
                readMappingBoth++;
            }else if(read2A.getReadUnmappedFlag() && read2B.getReadUnmappedFlag()){
                readMappingNone++;
            }

            //filter reads
            if(read1B.getReadUnmappedFlag() && read2B.getReadUnmappedFlag()){

                bamWriter.addAlignment(read1A);
                bamWriter.addAlignment(read2A);

            } else {

                readPairFilteredCount++;

            }


            if (readCount > 0 && readCount % 1000000 == 0) {
                logger.info(readCount + " records processed");
            }

        }

        bamWriter.close();

        logger.info(readPairCount + " read pairs processed");
        logger.info(readPairFilteredCount + " (" + oneDecimalPlace.format((double)readPairFilteredCount/(double)readPairCount*100) + "%) read pairs filtered");
        logger.info((readPairCount - readPairFilteredCount) + " (" + oneDecimalPlace.format((double)(readPairCount - readPairFilteredCount)/(double)readPairCount*100) + "%) read pairs written to ouput BAM");
        logger.info(readMappingAOnly + " (" + oneDecimalPlace.format((double)readMappingAOnly/(double)readCount*100) + "%) reads mapped in BAM file A only");
        logger.info(readMappingBOnly + " (" + oneDecimalPlace.format((double)readMappingBOnly/(double)readCount*100) + "%) reads mapped in BAM file B only");
        logger.info(readMappingBoth + " (" + oneDecimalPlace.format((double)readMappingBoth/(double)readCount*100) + "%) reads mapped in BAM file A and B");
        logger.info(readMappingNone + " (" + oneDecimalPlace.format((double)readMappingNone/(double)readCount*100) + "%) reads mapped in neither BAM file");

    }

    public static void main(String[] args) {

//        String bamB = "/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.phix.fa.readname.sorted.head.bam";
//        String bamA = "/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.rn4.ct.fa.readname.sorted.head.bam";

        String bamA = args[0];
        String bamB = args[1];
        String bamOut = args[2];

        new SubstractBam(new File(bamA), new File(bamB), new File(bamOut));
                        
    }

}
