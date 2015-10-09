package org.nxtgenutils.bsseq;

import net.sf.samtools.*;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.util.Map;
import java.util.HashMap;

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
 * Date: 09-Nov-2011
 * Time: 14:58:15
 */

/**
 * @author Michael Mueller
 */
public class BSMappingFiltrator {

    private int readCount = 0;
    private int readFiltered = 0;
    private Map<String, Integer> filterResults = new HashMap<String, Integer>();

    private BSReadFilterSet bsReadFilterSet;
    private File inputMapping;
    private File outputMapping;

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(BSMappingFiltrator.class);

    public BSMappingFiltrator(File inputMapping, File outputMapping, File inputReadFilterSet) throws IllegalStateException {

        this.bsReadFilterSet = new FileBSReadFilterSet(inputReadFilterSet);
        this.inputMapping = inputMapping;
        this.outputMapping = outputMapping;

        //create stats output file
        String currentWorkingDirectory = System.getProperty("user.dir");
        File statsOutput = new File(currentWorkingDirectory + "/" + inputMapping.getName().replace(".bam", ".filtered.bam.stats"));

        //create BAM file reader
        SAMFileReader samFileReader = new SAMFileReader(inputMapping);

        SAMFileWriter bamFileWriter;
        SAMFileHeader samFileHeader = samFileReader.getFileHeader();
        addFilterPropertiesToSamFileHeader(samFileHeader, bsReadFilterSet);

        //write file
        if (outputMapping != null) {

            bamFileWriter = new SAMFileWriterFactory().makeBAMWriter(samFileHeader,
                    true, outputMapping);
        }
        //or write to STDOUT
        else {
            bamFileWriter = new SAMFileWriterFactory().makeSAMWriter(samFileHeader, true, System.out);
        }

        filterReads(samFileReader, bamFileWriter, bsReadFilterSet);

        samFileReader.close();
        bamFileWriter.close();

        writeStats(statsOutput);


    }

    private void writeStats(File statsOut) {

        try {

            PrintWriter pw = new PrintWriter(statsOut);

            pw.println("input BAM file: " + inputMapping.getAbsolutePath());
            if (outputMapping != null) {
                pw.println("output BAM file: " + outputMapping.getAbsolutePath());
            } else {
                pw.println("output BAM file: STDOUT");
            }

            pw.println(readCount + " reads processed");
            pw.println(readCount - readFiltered + " reads pass filter set");
            pw.println("applied filter set:");
            for (String filter : filterResults.keySet()) {
                int count = filterResults.get(filter);
                pw.println(filter + " -> " + count + " reads did not pass filter");

            }

            pw.close();

        } catch (FileNotFoundException e) {

            logger.error("Unable to write filter stats to " + statsOut.getAbsolutePath() + ".", e);

        }


    }

    private void filterReads(SAMFileReader samFileReader, SAMFileWriter samFileWriter, BSReadFilterSet readFilterSet) {

        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMRecordIterator iterator = samFileReader.iterator();


        for (BSReadFilter filter : readFilterSet)
            filterResults.put(filter.toString(), 0);

        //iterate over entries in mapping
        while (iterator.hasNext()) {

            SAMRecord record = iterator.next();
            readCount++;

            //do not check unmapped reads
            if (record.getReadUnmappedFlag()) {

                samFileWriter.addAlignment(record);

            } else if (bsReadFilterSet.passesFilterSet(record)) {

                samFileWriter.addAlignment(record);

            } else {

                Map<String, Boolean> currentFilterResult = bsReadFilterSet.getFilterResults();
                for (String filter : currentFilterResult.keySet()) {
                    boolean passed = currentFilterResult.get(filter);
                    if (!passed) {
                        int count = filterResults.get(filter);
                        count++;
                        filterResults.put(filter, count);
                    }
                }
                readFiltered++;

            }

        }


    }

    private void addFilterPropertiesToSamFileHeader(SAMFileHeader samFileHeader, BSReadFilterSet readFilterSet) {

        for (BSReadFilter filter : readFilterSet) {
            samFileHeader.addComment("BS read filter " + filter.toString());
        }

    }
}
