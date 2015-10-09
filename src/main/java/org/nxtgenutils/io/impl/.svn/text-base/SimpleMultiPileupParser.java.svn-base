package org.nxtgenutils.io.impl;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.NoSuchElementException;

import org.nxtgenutils.util.FileTypeDeterminator;
import org.nxtgenutils.io.MultiPileup;
import org.nxtgenutils.io.MultiPileupParser;

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
 * Date: 26-Oct-2011
 * Time: 13:21:17
 */

/**
 * Implementation of MultiPileupParser interface.
 *
 * @author Michael Mueller
 */
public class SimpleMultiPileupParser implements MultiPileupParser {

    private File mPileupFile;

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(SimpleMultiPileupParser.class);

    /**
     * Constructs a parser instance for a multi-pileup file generated with samtools.
     *
     * @param multiPileUpFile the mpileup output file
     */
    public SimpleMultiPileupParser(File multiPileUpFile) {

        this.mPileupFile = multiPileUpFile;

    }

    /**
     * Returns a MultiPileup iterator.
     *
     * @return the MultiPileup iterator
     */
    public Iterator<MultiPileup> iterator() {

        PileupRecordIterator retVal = null;

        try {
            retVal = new PileupRecordIterator(mPileupFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return retVal;
    }

    /**
     * Iterator implementation to iterate over records in a multi-pileup file file
     * and create MultiPileup objects.
     */
    class PileupRecordIterator implements Iterator<MultiPileup> {

        BufferedReader br;
        String nextLine = null;

        PileupRecordIterator(File vcfFile) throws IOException {
            br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(vcfFile)));
            nextLine = br.readLine();
        }

        @Override
        public boolean hasNext() {
            return nextLine != null;
        }

        @Override
        public MultiPileup next() {

            if (nextLine == null) {
                throw new NoSuchElementException("End of pileup file reached.");
            }

            MultiPileup retVal;

            StringTokenizer tokenizer = new StringTokenizer(nextLine, "\t");

            String sequenceName = tokenizer.nextToken();
            int start = Integer.parseInt(tokenizer.nextToken());
            String referenceBase = tokenizer.nextToken();

            retVal = new MultiPileupImpl(sequenceName, start, referenceBase);
            logger.error(nextLine);
            while (tokenizer.hasMoreElements()) {

                int coverage = Integer.parseInt(tokenizer.nextToken());
                String bases = tokenizer.nextToken();
                String qualities = tokenizer.nextToken();

                retVal.addPileup(new PileupImpl(retVal, coverage, bases, qualities));

            }

            retVal.setRecordString(nextLine);

            try {
                nextLine = br.readLine();

                if (nextLine == null) {
                    br.close();
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
