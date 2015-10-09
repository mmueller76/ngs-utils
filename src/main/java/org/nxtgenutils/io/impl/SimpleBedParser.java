package org.nxtgenutils.io.impl;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

import org.nxtgenutils.util.FileTypeDeterminator;
import org.nxtgenutils.Strand;
import org.nxtgenutils.io.BedRecord;

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
 * Date: 25-Oct-2011
 * Time: 17:06:18
 * To change this template use File | Settings | File Templates.
 */

/**
 * Minimal BED file parser implementation.
 *
 * @author Michael Mueller
 */
public class SimpleBedParser {

    private File bedFile;
    private String trackLine;
    private boolean hasTrackLine = false;
    
    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(SimpleBedParser.class);

    /**
     * Constructs a parser instance for a BED formatted input file.
     *
     * @param bedFile the BED file
     */
    public SimpleBedParser(File bedFile) {

        this.bedFile =  bedFile;
        this.trackLine = parseTrackLine();

    }

    /**
     * Returns the track line.
     *
     * @return the track line if present, null otherwise
     */
    public String getTrackLine() {
        return trackLine;
    }

    /**
     * Finds and returns the BED track line.
     *
     * @return the BED track line if present, null otherwise
     */
    private String parseTrackLine() {

        String retVal = null;

        try {

            //open input file
            BufferedReader br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(bedFile)));

            //find track line
            String line;
            if ((line = br.readLine()) != null && line.startsWith("track")) {
                retVal = line;
                hasTrackLine = true;
            }

            //close input file
            br.close();

        } catch (IOException e) {
            logger.error("Exception while parsing BED file track line.", e);
        }

        return retVal;
    }

    /**
     * Returns a BedRecord iterator.
     *
     * @return the BedRecord iterator
     */
    public Iterator<BedRecord> iterator() {

        BedRecordIterator retVal = null;

        try {
            retVal = new BedRecordIterator(bedFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return retVal;
    }

    /**
     * Iterator implementation to iterate over records in a BED formatted file
     * and create BedRecord objects.
     */
    class BedRecordIterator implements Iterator<BedRecord> {

        BufferedReader br;
        String nextLine = null;
        Map<Integer, String> sampleNames = new HashMap<Integer, String>();

        BedRecordIterator(File bedFile) throws IOException {

            br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(bedFile)));

            //skip track line
            if (hasTrackLine) {
                br.readLine();
            }

            nextLine = br.readLine();

        }

        @Override
        public boolean hasNext() {
            return nextLine != null;
        }

        @Override
        public BedRecord next() {

            BedRecord retVal;

            StringTokenizer tokenizer = new StringTokenizer(nextLine, "\t");

            String chr = tokenizer.nextToken();
            int start = Integer.parseInt(tokenizer.nextToken());
            int end = Integer.parseInt(tokenizer.nextToken());
            String name = tokenizer.nextToken();

            retVal = new BedRecordImpl(chr, start, end, name);

            int score = -1;
            if (tokenizer.hasMoreTokens()) {
                Integer.parseInt(tokenizer.nextToken());
                retVal.setScores(score);
            }

            Strand strand = Strand.UNKNOWN;
            if (tokenizer.hasMoreTokens()) {
                String strandString = tokenizer.nextToken();
                if (strandString.equals("+")) {
                    strand = Strand.FORWARD;
                } else if (strandString.equals("-")) {
                    strand = Strand.REVERSE;
                }
                retVal.setStrand(strand);
            }

            try {
                nextLine="";
                //go to next line skipping empty lines
                while(nextLine.equals("")){
                    nextLine = br.readLine();
                    //leave loop if no more line in file
                    if(nextLine == null){
                        break;
                    }
                }
            } catch (IOException e) {
                logger.error(e);
            }

            return retVal;

        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }


    }   

}
