package org.nxtgenutils.io.impl;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

import org.nxtgenutils.util.FileTypeDeterminator;
import org.nxtgenutils.io.VCFRecord;
import org.nxtgenutils.io.VCFParser;

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
 * Date: 06-Oct-2011
 * Time: 16:27:20
 */

/**
 * Minimal implementation of the VCFParser interface.
 *
 * @author Michael Mueller
 */
public class SimpleVCFParser implements VCFParser {

    private File vcfFile;
    private Set<String> headerLines;
    private double minorAlleleFrequencyCutoff = -1;

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(SimpleVCFParser.class);


    /**
     * Constructs a parser instance for a VCF file.
     *
     * @param vcfFile the mpileup output file
     */
    public SimpleVCFParser(File vcfFile) {

        this.vcfFile = vcfFile;
        this.headerLines = parseHeaderLines();

    }

    /**
     * Returns the header lines of the VCF file.
     *
     * @return a set of header lines
     */
    public Set<String> getHeaderLines() {
        return headerLines;
    }

    /**
     * Finds and extracts header lines from the VCF file.
     *
     * @return a set of header lines
     */
    private Set<String> parseHeaderLines() {
        Set<String> retVal = new LinkedHashSet<String>();

        try {
            BufferedReader br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(vcfFile)));

            String line;
            int lineCount = 0;
            while ((line = br.readLine()) != null && line.startsWith("#")) {
                lineCount++;
                retVal.add(line);
            }

            if (lineCount == 0) {
                throw new IOException("VCF file does not contain header information.");
            }

            br.close();

        } catch (IOException e) {
            logger.error("Exception while parsing VCF file header.", e);
        }

        return retVal;
    }

    /**
     * Returns a VCFRecord iterator.
     *
     * @return the BedRecord iterator
     */
    public Iterator<VCFRecord> iterator() {
        VCFRecordIterator retVal = null;

        try {
            retVal = new VCFRecordIterator(vcfFile);
        } catch (IOException e) {
            e.printStackTrace();
        }

        return retVal;
    }

    /**
     * Iterator implementation to iterate over records in a VCF formatted file
     * and create VCFRecord objects.
     */
    class VCFRecordIterator implements Iterator<VCFRecord> {

        BufferedReader br;
        String nextLine = null;
        Map<Integer, String> sampleNames = new HashMap<Integer, String>();

        VCFRecordIterator(File vcfFile) throws IOException {
            br = new BufferedReader(FileTypeDeterminator.determineFileType(new FileInputStream(vcfFile)));
            //br = new BufferedReader(new FileReader(vcfFile));
            boolean headerFound = false;
            while ((nextLine = br.readLine()) != null && !headerFound) {

                if (nextLine.startsWith("#CHROM")) {

                    headerFound = true;
                    StringTokenizer tokenizer = new StringTokenizer(nextLine, "\t");
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    tokenizer.nextToken();
                    int sampleIdx = 1;
                    while (tokenizer.hasMoreTokens()) {
                        String sampleName = tokenizer.nextToken();
                        sampleNames.put(sampleIdx, sampleName);
                        sampleIdx++;
                    }
                }
            }
            nextLine = br.readLine();
        }

        @Override
        public boolean hasNext() {
            return nextLine != null;
        }

        @Override
        public VCFRecord next() {

            VCFRecord retVal = new VCFRecordImpl();

            retVal.setVcfString(nextLine);

            StringTokenizer tokenizer = new StringTokenizer(nextLine, "\t");

            String chr = tokenizer.nextToken();
            int pos = Integer.parseInt(tokenizer.nextToken());
            tokenizer.nextToken();
            String refAllele = tokenizer.nextToken();
            String altAllele = tokenizer.nextToken();
            retVal.setChromsome(chr);
            retVal.setPosition(pos);
            retVal.setReferenceAllele(refAllele);
            retVal.setAlternativeAllele(altAllele);

            String qualityString = tokenizer.nextToken();
            String filterValue = tokenizer.nextToken();
            retVal.setFilterValue(filterValue);
            tokenizer.nextToken();
            tokenizer.nextToken();

            int sampleIdx = 1;
            while (tokenizer.hasMoreTokens()) {

                String sampleString = tokenizer.nextToken();
                StringTokenizer sampleTokenizer = new StringTokenizer(sampleString, ":");
                String genotype = sampleTokenizer.nextToken();

                if (!genotype.equals("./.")) {

                    String allelicDepth = sampleTokenizer.nextToken();
                    int allelicDepthReference = Integer.parseInt(allelicDepth.split(",")[0]);
                    int allelicDepthAlternative = Integer.parseInt(allelicDepth.split(",")[1]);

                    //apply genotype correction if heterozygous allele
                    //and minimum allele frequency cutoff specified
                    if (genotype.equals("1/1") && minorAlleleFrequencyCutoff != -1) {
                        genotype = correctGenotype(allelicDepthReference, allelicDepthAlternative);
                    }

                }

                String sampleName = sampleNames.get(sampleIdx);
                retVal.addSampleRecord(new VCFSampleRecordImpl(sampleName, genotype));
                sampleIdx++;

            }


            try {
                nextLine = br.readLine();
            } catch (IOException e) {
                logger.error(e);
            }

            return retVal;  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }


    }

    /**
     * Corrects heterozygous genotypes to homozygous reference/alternative allele
     * if the allele frequency of the minor allele is below a specified allele
     * frequence cutoff.
     *
     * @param allelicDepthReference   allele frequency of the reference allele
     * @param allelicDepthAlternative allele frequency of the alternative allele
     * @return genotype
     */
    private String correctGenotype(int allelicDepthReference, int allelicDepthAlternative) {

        String retVal = "0/1";
        int allelicDepthTotal = allelicDepthReference + allelicDepthAlternative;

        if (allelicDepthReference < allelicDepthAlternative && allelicDepthReference / allelicDepthTotal < minorAlleleFrequencyCutoff) {
            //correct SNP to homozygous alternative allele
            retVal = "1/1";
        } else if (allelicDepthReference > allelicDepthAlternative && allelicDepthAlternative / allelicDepthTotal < minorAlleleFrequencyCutoff) {
            //correct SNP to homozygous reference allele
            retVal = "0/0";
        }

        return retVal;

    }

    /**
     * Returns the minimum allele frequency at which a SNP called as heterozygous
     * in the VCF file will be returned as heterozyogous.
     *
     * @param cutoff the relative frequency
     */
    @Override
    public void setHeterozygousMinorAlleleFrequencyCutoff(double cutoff) {
        this.minorAlleleFrequencyCutoff = cutoff;
    }

    /**
     * Sets the minimum allele frequency at which a SNP called as heterozygous
     * in the VCF file will be returned as heterozyogous.
     *
     * @return the relative frequency
     */
    @Override
    public double getHeterozygousMinorAlleleFrequencyCutoff() {
        return minorAlleleFrequencyCutoff;
    }

}
