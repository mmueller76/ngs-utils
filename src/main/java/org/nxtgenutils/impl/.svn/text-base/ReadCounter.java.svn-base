package org.nxtgenutils.impl;

import java.io.*;
import java.util.Iterator;
import java.util.Set;
import java.util.HashSet;

import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;
import org.apache.log4j.Logger;
import org.nxtgenutils.io.BedRecord;
import org.nxtgenutils.io.impl.SimpleBedParser;

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
 * Date: 18-Oct-2011
 * Time: 10:19:10
 */

/**
 * @author Michael Mueller
 */
public class ReadCounter {

    private static Logger logger = Logger.getLogger(ReadCounter.class);

    public ReadCounter(File mappingInput, File regionInput, File countOut, int minMappingQuality, boolean mergeOverlappingReadPairs) {

        File mappingIndex = new File(mappingInput + ".bai");

        try {

            PrintWriter coverageOut = new PrintWriter(countOut);
//            BufferedReader regions = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(regionInput))));
            SimpleBedParser regions = new SimpleBedParser(regionInput);

            int intervalCount = 0;

            Iterator<BedRecord> regionsIterator = regions.iterator();
            while(regionsIterator.hasNext()){

                BedRecord region = regionsIterator.next();

                String regionChr = region.getChromosome();

                int regionStart = region.getStart()+1;
                int regionEnd = region.getEnd();
                String regionName = region.getName();

                intervalCount++;
                if (intervalCount % 1000 == 0) {
                    System.out.println(intervalCount + " regions processed");
                }

                Interval regionInterval = new Interval(regionChr, regionStart, regionEnd);

                SAMFileReader samFileReader = new SAMFileReader(mappingInput, mappingIndex);
                samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

                Iterator<SAMRecord> it = samFileReader.queryOverlapping(regionChr, regionStart, regionEnd);

                //retrieve reads in query name order to
                //facilitate detection of read pairs
                //mapping to the same region
//                if(mergeOverlappingReadPairs){
//                    it = samFileReader.queryOverlapping(regionChr, regionStart, regionEnd).assertSorted(SAMFileHeader.SortOrder.queryname);
//                }
//                //... else read order is not important
//                else {
//                    it = samFileReader.queryOverlapping(regionChr, regionStart, regionEnd).assertSorted(SAMFileHeader.SortOrder.queryname);
//                }

                int readCount = 0;
                Set<String> previousReadIds = new HashSet<String>();
                while (it.hasNext()) {

                    SAMRecord read = it.next();
                    String readId = read.getReadName();
                    if (read.getMappingQuality() >= minMappingQuality) {

                        for (AlignmentBlock block : read.getAlignmentBlocks()) {

                            int blockStart = block.getReferenceStart();
                            int blockEnd = block.getReferenceStart() + block.getLength() - 1;

                            Interval blockInterval = new Interval(regionChr, blockStart, blockEnd);
                            //if overlap found...
                            if (blockInterval.getIntersectionLength(regionInterval) > 0) {

                                //...if we do not care about overlapping
                                //read pairs
                                if(!mergeOverlappingReadPairs){
                                    readCount++;
                                }

                                //...if we do care about overlapping
                                //read pairs we check first if a read
                                //from the same fragment has already been seen
                                else if (!previousReadIds.contains(readId)){
                                    readCount++;
                                }

                                previousReadIds.add(readId);

                                break;
                                
                            }
                        }

                    }
                }

                coverageOut.println(regionChr + "\t" +
                        regionStart + "\t" +
                        regionEnd + "\t" +
                        regionName + "\t" +
                        readCount);
                coverageOut.flush();

                //close BAM file
                samFileReader.close();

            }

            coverageOut.close();

        } catch (IOException e) {
            logger.error(e);
        }


    }

}
