package org.nxtgenutils.io.impl;

import org.apache.log4j.Logger;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

import org.nxtgenutils.util.FileTypeDeterminator;
import org.nxtgenutils.io.MethylationProfileParser;
import org.nxtgenutils.io.MethylationProfileRecord;
import org.nxtgenutils.Strand;

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
 * Date: 18-Nov-2011
 * Time: 12:24:51
 */

/**
 * @author Michael Mueller
 */
public class SimpleMethylationProfileParser implements MethylationProfileParser {

    private File methylationProfileFile;
    private String headerLine;
    private boolean hasHeaderLine = false;
    private List<String> sampleNames;

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(SimpleMethylationProfileParser.class);


    public SimpleMethylationProfileParser(File methylationProfile, List<String> sampleNames) {

        this.methylationProfileFile = methylationProfile;
        this.headerLine = parseHeaderLine();
        this.sampleNames=sampleNames;

    }

    public SimpleMethylationProfileParser(File methylationProfile) {

        this.methylationProfileFile = methylationProfile;
        this.headerLine = parseHeaderLine();
        if (hasHeaderLine) {
            sampleNames=extraxtSampleNamesFromHeaderLiner();
        } else {
            sampleNames=createSampleNames();
        }

    }

    public String getHeaderLine() {
        return headerLine;
    }

    private List<String> createSampleNames() {



        String firstLine = null;

        try {

            BufferedReader br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(methylationProfileFile)));
            firstLine = br.readLine();
            br.close();

        } catch (IOException e) {
            logger.error("Exception while determining sample count in methylation profile file.", e);
        }

        List<String> retVal = new ArrayList<String>();

        if (firstLine != null) {
            String[] columns = firstLine.split("\t");
            int sampleCount = (columns.length - 6) / 8;
            for (int s = 1; s <= sampleCount; s++) {
                retVal.add("" + s);
            }
        }

        return retVal;

    }

    private List<String> extraxtSampleNamesFromHeaderLiner() {

        List<String> retVal = new ArrayList<String>();

        String[] columns = headerLine.split("\t");
        int sampleCounter = 0;
        for (String column : columns) {

            if (column.endsWith("_smpl_context")) {
                String sampleName = column.replace("_smpl_context", "");
                retVal.add(sampleName);
            } else if (column.endsWith("smpl_context")) {
                sampleCounter++;
                retVal.add("" + sampleCounter);
            }

        }

        return retVal;

    }

    public List<String> getSampleNames() {
        return sampleNames;
    }

    public void setSampleNames(List<String> sampleNames) {
        this.sampleNames = sampleNames;
    }

    public int getSampleCounts() {
        return sampleNames.size();
    }

    private String parseHeaderLine() {

        String retVal = null;

        try {

            BufferedReader br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(methylationProfileFile)));

            String line;
            if ((line = br.readLine()) != null && line.startsWith("chr\tstart")) {
                retVal = line;
                hasHeaderLine = true;
            }

            br.close();

        } catch (IOException e) {
            logger.error("Exception while parsing VCF file header.", e);
        }

        return retVal;
    }

    public Iterator<MethylationProfileRecord> iterator() {
        MethylationProfileRecordIterator retVal = null;

        try {
            retVal = new MethylationProfileRecordIterator(methylationProfileFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return retVal;
    }

    class MethylationProfileRecordIterator implements Iterator<MethylationProfileRecord> {

        BufferedReader br;
        String nextLine = null;
        Map<Integer, String> sampleNames = new HashMap<Integer, String>();

        MethylationProfileRecordIterator(File methylationProfileFile) throws IOException {
            br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(methylationProfileFile)));

            //skip header line
            if (hasHeaderLine) {
                br.readLine();
            }

            nextLine = br.readLine();
        }

        @Override
        public boolean hasNext() {
            return nextLine != null;
        }

        @Override
        public MethylationProfileRecord next() {

            // 1 chr
        // 2 start
        // 3 end
        // 4 strand
        // 5 ref_context
        // 6 repeat_masked

        // 1 smpl_context
        // 2 c_count
        // 3 ct_count
        // 4 non_ct_count
        // 5 m%
        // 6 score
        // 7 snp
        // 8 indels

            MethylationProfileRecord retVal;

            StringTokenizer tokenizer = new StringTokenizer(nextLine, "\t");

            String chr = tokenizer.nextToken();
            int start = Integer.parseInt(tokenizer.nextToken());
            int end = Integer.parseInt(tokenizer.nextToken());
            String strandString = tokenizer.nextToken();
            Strand strand = Strand.FORWARD;
            if(strandString.equals("-")){
                strand = Strand.REVERSE;
            }
            String referenceContext = tokenizer.nextToken();
            String repeatMaskedString = tokenizer.nextToken();
            boolean repeatMasked = false;
            if(repeatMaskedString.equals("1")){
                repeatMasked = true;
            }

            retVal = new MethylationProfileRecordImpl(chr, start, end, strand, referenceContext, repeatMasked);

            int sampleIdx = 1;
            while (tokenizer.hasMoreTokens()) {

                String sampleName = sampleNames.get(sampleIdx);

                String sampleContext = tokenizer.nextToken();
                int cCount = Integer.parseInt(tokenizer.nextToken());
                int ctCount = Integer.parseInt(tokenizer.nextToken());
                int nonCTCount = Integer.parseInt(tokenizer.nextToken());
                //int methylationLevel = Integer.parseInt(tokenizer.nextToken());
                int score = Integer.parseInt(tokenizer.nextToken());
                int snpType = Integer.parseInt(tokenizer.nextToken());
                String indelString =tokenizer.nextToken();
                boolean indels = false;
                if(indelString.equals("1")){
                    indels = true;
                }

                retVal.addSampleRecord(new MethylationProfileSampleRecordImpl(retVal, sampleContext, cCount, ctCount, nonCTCount, sampleName, score, snpType, indels));

                sampleIdx++;

            }

            try {
                nextLine = br.readLine();
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
