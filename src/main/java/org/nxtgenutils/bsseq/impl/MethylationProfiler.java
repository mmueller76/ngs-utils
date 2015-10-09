package org.nxtgenutils.bsseq.impl;

import org.nxtgenutils.io.impl.SimpleMultiPileupParser;
import org.nxtgenutils.io.impl.MethylationProfileRecordImpl;
import org.nxtgenutils.io.impl.MethylationProfileSampleRecordImpl;
import org.nxtgenutils.io.*;
import org.nxtgenutils.Strand;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPOutputStream;
import java.util.regex.Pattern;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import org.apache.log4j.Logger;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.commons.math.MathException;

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
 * Time: 11:44:35
 */

/**
 * @author Michael Mueller
 */
public class MethylationProfiler {

    //sequence name
    //position
    //reference base
    //strand
    //context
    //sample x coverage
    //sample x C
    //sample x T
    //sample x Other

    static private String profileColumnLabels = "chr\t" +
            "start\t" +
            "end\t" +
            "strand\t" +
            "ref_context\t" +
            "repeat_masked";

    static private String[] sampleColumnLabels = {"smpl_context",
            "c_count",
            "ct_count",
            "non_ct_count",
            "m%",
            "score",
            "snp",
            "indels"};

    static Pattern cg = Pattern.compile("CG");
    static Pattern tg = Pattern.compile("TG");
    static Pattern ag_g = Pattern.compile("[AG]G");
    static Pattern c_atc = Pattern.compile("C[ATC]");
    static Pattern ca = Pattern.compile("CA");
    static Pattern c_ct = Pattern.compile("C[CT]");
    static Pattern c_act = Pattern.compile("C[ACT]");
    static Pattern t_act = Pattern.compile("T[ACT]");
    static Pattern ag_act = Pattern.compile("AG[ACT]");
    static Pattern atg_atc = Pattern.compile("[ATG][ATC]");

    private int outputCounterForwardCpG = 0;
    private int outputCoutnerForwardNonCpG = 0;
    private int outputCounterReverseCpG = 0;
    private int outputCoutnerReverseNonCpG = 0;

    private Map<String, Map<String, Long>> nonCpGContextCCount;
    private Map<String, Map<String, Long>> nonCpGContextTCount;

    private List<String> sampleNames;

    private Set<String> nonCpGContexts = new HashSet<String>();

    private static String NON_CPG_CONTEXT = "non-CpG-context";

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(MethylationProfiler.class);


    /**
     * @param mPileupSample
     * @param outputFile
     * @param estimateBisulfiteConversionRateFrom
     *
     * @param sampleNames
     */
    public MethylationProfiler(File mPileupSample, File outputFile, String estimateBisulfiteConversionRateFrom, List<String> sampleNames) {
        this(mPileupSample, null, outputFile, estimateBisulfiteConversionRateFrom, sampleNames);
    }

    /**
     * @param mPileupSample
     * @param mPileupControl
     * @param outputFile
     * @param estimateBisulfiteConversionRateFrom
     *
     * @param sampleNames
     */
    public MethylationProfiler(File mPileupSample, File mPileupControl, File outputFile, String estimateBisulfiteConversionRateFrom, List<String> sampleNames) {

        if (estimateBisulfiteConversionRateFrom != null && estimateBisulfiteConversionRateFrom.equals("")) {
            estimateBisulfiteConversionRateFrom = NON_CPG_CONTEXT;
        }

        int sampleCount = getSampleCount(mPileupSample);

        if (sampleNames != null) {
            this.sampleNames = sampleNames;
        } else {
            this.sampleNames = generateSampleNames(sampleCount);
        }

        nonCpGContexts.add("CT");
        nonCpGContexts.add("CA");
        nonCpGContexts.add("CC");
        nonCpGContexts.add("total");

        initializeNonCpGCountMaps(this.sampleNames);

        Map<String, Double> bisulfiteConversionRate = null;
        if (estimateBisulfiteConversionRateFrom != null) {

            logger.info("Estimating bisulfite conversionr rate...");
            File mPileupConversionRate = mPileupSample;
            if (mPileupControl != null) {
                mPileupConversionRate = mPileupControl;
            }

            bisulfiteConversionRate = estimateBisulfiteConversionRate(mPileupConversionRate, estimateBisulfiteConversionRateFrom);

            NumberFormat formatter = new DecimalFormat("#.#");
            for (String name : this.sampleNames) {
                double rate = bisulfiteConversionRate.get(name);
                if (rate < 0) {
                    logger.info(name + ": -");
                } else {
                    logger.info(name + ": " + formatter.format(rate * 100.0) + "%");
                }
            }

            //print conversion rates to file
            String conversionRateOutputPath = outputFile.getAbsolutePath() + "." + NON_CPG_CONTEXT + ".conversion.rate.tsv";
            if (!estimateBisulfiteConversionRateFrom.equals(NON_CPG_CONTEXT)) {
                conversionRateOutputPath = outputFile.getAbsolutePath() + "." + estimateBisulfiteConversionRateFrom + ".conversion.rate.tsv";
            }

            logger.info("Writing bisulfite conversionr rate to " + conversionRateOutputPath + "...");
            writeConversionRate(bisulfiteConversionRate, conversionRateOutputPath);
            logger.info("--------------------------------------------------------");

        }

        MultiPileupParser pileupParser = new SimpleMultiPileupParser(mPileupSample);
        Iterator<MultiPileup> pileupIterator = pileupParser.iterator();

        GZIPOutputStream outputStream = null;

        try {
            outputStream = new GZIPOutputStream(new FileOutputStream(outputFile.getAbsolutePath() + ".gz"));
        } catch (IOException e) {
            logger.error("Exception while opening output file: " + outputFile.getAbsolutePath(), e);
            System.exit(1);
        }

        PrintWriter pw = new PrintWriter(outputStream);

        logger.info("Generating methylation profiles...");

        MultiPileup previousMultiPileup = null;
        int lineCounter = 0;
        while (pileupIterator.hasNext()) {

            lineCounter++;

            if (lineCounter % 1000000 == 0) {
                logger.info(lineCounter + " lines processed");
            }

            MultiPileup currentMultiPileup = pileupIterator.next();

            if (previousMultiPileup != null) {

                MethylationProfileRecord profileRecordForwardStrand = generateProfileRecord(previousMultiPileup, currentMultiPileup, Strand.FORWARD, this.sampleNames, estimateBisulfiteConversionRateFrom);
                MethylationProfileRecord profileRecordReverseStrand = generateProfileRecord(previousMultiPileup, currentMultiPileup, Strand.REVERSE, this.sampleNames, estimateBisulfiteConversionRateFrom);

                if (profileRecordForwardStrand != null) {

                    if (estimateBisulfiteConversionRateFrom != null) {
                        scoreProfileRecord(profileRecordForwardStrand, bisulfiteConversionRate);
                    }

                    pw.println(profileRecordForwardStrand.formatProfileRecord());
                    pw.flush();

                }

                if (profileRecordReverseStrand != null) {

                    if (estimateBisulfiteConversionRateFrom != null) {
                        scoreProfileRecord(profileRecordReverseStrand, bisulfiteConversionRate);
                    }

                    pw.println(profileRecordReverseStrand.formatProfileRecord());
                    pw.flush();

                }

            } else {

                pw.print(profileColumnLabels);
                for (int i = 0; i < sampleCount; i++) {
                    for(String label : sampleColumnLabels){
                        pw.print("\t" + this.sampleNames.get(i) + "_" + label);
                    }
                }
                pw.println();

            }

            previousMultiPileup = currentMultiPileup;


        }

        pw.close();

        try {
            outputStream.close();
        } catch (IOException e) {
            logger.error("Exception while closing output file: " + outputFile.getAbsolutePath(), e);
        }

        logger.info(lineCounter + " lines processed");
        logger.info("--------------------------------------------------------");
        logger.info("non-CpG conversion rate");

        printNonCpGConversionRate();

        String conversionRateOutputPath = outputFile.getAbsolutePath() + "." + NON_CPG_CONTEXT + ".conversion.rate.detailed.tsv";
        writeNonCpGConversionRate(conversionRateOutputPath);

    }


    private void writeConversionRate(Map<String, Double> bisulfiteConversionRate, String outputFilePath) {

        try {

            NumberFormat formatter = new DecimalFormat("#.#");

            PrintWriter conversionRateOutput = new PrintWriter(outputFilePath);
            for (String sample : this.sampleNames) {
                double rate = bisulfiteConversionRate.get(sample);
                if (rate < 0) {
                    conversionRateOutput.println(sample + "\t-1");
                    conversionRateOutput.flush();

                } else {
                    conversionRateOutput.println(sample + "\t" + formatter.format(rate * 100.0));
                    conversionRateOutput.flush();
                }
            }
            conversionRateOutput.close();

        } catch (FileNotFoundException e) {
            logger.error("Exception while writing conversion rate to " + outputFilePath + ".", e);
        }

    }

    /**
     *
     */
    private void printNonCpGConversionRate() {

        NumberFormat formatter = new DecimalFormat("#.#");
        logger.info("context\tsample\trate");

        for (String context : nonCpGContexts) {

            Map<String, Long> cCount = nonCpGContextCCount.get(context);
            Map<String, Long> tCount = nonCpGContextTCount.get(context);

            Map<String, Double> bisulfiteConversionRate = calculateBisfulfiteConversionRate(cCount, tCount);

            for (String sampleName : sampleNames) {

                double rate = bisulfiteConversionRate.get(sampleName);
                if (rate < 0) {
                    logger.info(context + "\t" + sampleName + "\t-1");
                } else {
                    logger.info(context + "\t" + sampleName + "\t" + formatter.format(rate * 100.0) + "%");
                }

            }

        }

    }

    private void writeNonCpGConversionRate(String outputFilePath) {

        try {

            NumberFormat formatter = new DecimalFormat("#.#");

            PrintWriter conversionRateOutput = new PrintWriter(outputFilePath);

            for (String context : nonCpGContexts) {

                Map<String, Long> cCount = nonCpGContextCCount.get(context);
                Map<String, Long> tCount = nonCpGContextTCount.get(context);

                Map<String, Double> bisulfiteConversionRate = calculateBisfulfiteConversionRate(cCount, tCount);

                for (String sampleName : sampleNames) {

                    double rate = bisulfiteConversionRate.get(sampleName);
                    if (rate < 0) {
                        conversionRateOutput.println(context + "\t" + sampleName + "\t-1");
                        conversionRateOutput.flush();
                        logger.info(context + "\t" + sampleName + "\t-1");
                    } else {
                        conversionRateOutput.println(context + "\t" + sampleName + "\t" + formatter.format(rate * 100.0));
                        conversionRateOutput.flush();
                        logger.info(context + "\t" + sampleName + "\t" + formatter.format(rate * 100.0) + "%");
                    }

                }

            }

            conversionRateOutput.close();

        } catch (FileNotFoundException e) {
            logger.error("Exception while writing conversion rate to " + outputFilePath + ".", e);
        }

    }

    /**
     * @param sampleNames
     */
    private void initializeNonCpGCountMaps(List<String> sampleNames) {

        nonCpGContextCCount = new HashMap<String, Map<String, Long>>();
        nonCpGContextTCount = new HashMap<String, Map<String, Long>>();

        for (String context : nonCpGContexts) {

            nonCpGContextCCount.put(context, new HashMap<String, Long>());
            nonCpGContextTCount.put(context, new HashMap<String, Long>());

            for (String sampleName : sampleNames) {

                nonCpGContextCCount.get(context).put(sampleName, (long) 0);
                nonCpGContextTCount.get(context).put(sampleName, (long) 0);

            }
        }

    }

    /**
     * @param profileRecord
     * @param conversionRate
     */
    private void scoreProfileRecord(MethylationProfileRecord profileRecord, Map<String, Double> conversionRate) {

        for (MethylationProfileSampleRecord sampleRecord : profileRecord.getSampleRecords()) {

            double rate = conversionRate.get(sampleRecord.getSampleName());
            int score = -1;
            boolean rateGreateZero = false;
            if (rate > 0) {
                rateGreateZero = true;
                int trials = sampleRecord.getCtCoverage();
                int successes = sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                BinomialDistribution bd = new BinomialDistributionImpl(trials, rate);
                double p = -1;
                try {

                    // BinomialDistributionImpl does return 0 of the p value is >10^-16.
                    // In this case we will assign a score of 161 which is one more then
                    // the maximum score that can be calculated from the p values returned
                    // by BinomialDistributionImpl and equivivalent to a p-value of 8.9 * 10^-17
                    p = bd.cumulativeProbability(successes);
                    if (p == 0.0) {
                        score = 161;
                    } else {
                        score = (int) Math.round(-10 * Math.log10(p));
                    }

                } catch (MathException e) {
                    logger.error("Exception while scoring methylation level.", e);
                }
                if (score == -1) {
                    logger.error(sampleRecord.getCCount() + "\t" + sampleRecord.getCtCoverage() + "\t" + rate + "\t" + score + "\t" + rateGreateZero + "\t" + successes + "\t" + trials + "\t" + p);
                }
            }
            sampleRecord.setScore(score);


        }

    }

    /**
     * Generates a list of sample names.
     *
     * @param sampleCount the number of samples in the input pileup file
     * @return a list of sample names
     */
    private List<String> generateSampleNames(int sampleCount) {

        List<String> retVal = new ArrayList<String>();
        for (int i = 0; i < sampleCount; i++) {
            retVal.add("s" + (i + 1));
        }
        return retVal;

    }

    /**
     * Returns the number of samples in a pileup file. For bisulfite sequencing
     * data there are two pileups, the reverse and forward strand pileup, for each
     * sample. Thus the number of samples equals the number of pileups devided by
     * two.
     *
     * @param mPileup the pileup file
     * @return the number of samples in the pileup file
     */
    private int getSampleCount(File mPileup) {

        int retVal = 0;

        MultiPileupParser pileupParser = new SimpleMultiPileupParser(mPileup);
        Iterator<MultiPileup> pileupIterator = pileupParser.iterator();
        if (pileupIterator.hasNext()) {
            MultiPileup currentMultiPileup = pileupIterator.next();
            retVal = currentMultiPileup.getPileupCount() / 2;
        }

        return retVal;
    }

    /**
     * @param sampleCCount
     * @param sampleTCount
     * @return
     */
    private Map<String, Double> calculateBisfulfiteConversionRate(Map<String, Long> sampleCCount, Map<String, Long> sampleTCount) {

        Map<String, Double> retVal = new HashMap<String, Double>();

        for (String name : sampleNames) {
            long c = sampleCCount.get(name);
            long t = sampleTCount.get(name);
            long ct = c + t;
            if (ct == 0) {
                retVal.put(name, -1.0);
            } else {
                double rate = (double) t / (double) ct;
                retVal.put(name, rate);
            }
        }

        return retVal;
    }

    /**
     * @param mPileup
     * @param estimateBisulfiteConversionRateFrom
     *
     * @return
     */
    private Map<String, Double> estimateBisulfiteConversionRate(File mPileup, String estimateBisulfiteConversionRateFrom) {

        Map<String, Double> retVal = new HashMap<String, Double>();

        if (estimateBisulfiteConversionRateFrom == null) {
            return retVal;
        }

        Map<String, Long> cCount = new HashMap<String, Long>();
        Map<String, Long> tCount = new HashMap<String, Long>();

        //initialise map keys
        for (String name : sampleNames) {
            cCount.put(name, (long) 0);
            tCount.put(name, (long) 0);
        }

        MultiPileupParser pileupParser = new SimpleMultiPileupParser(mPileup);
        Iterator<MultiPileup> pileupIterator = pileupParser.iterator();

        int lineCounter = 0;
        int usableLineCounter = 0;
        int controlSampleEntryCounter = 0;

        MultiPileup previousMultiPileup = null;
        while (pileupIterator.hasNext()) {

            lineCounter++;

            if (lineCounter % 1000000 == 0) {
                logger.info(lineCounter + " lines processed");
            }

            MultiPileup currentMultiPileup = pileupIterator.next();

            if (previousMultiPileup != null) {

                if (estimateBisulfiteConversionRateFrom.equals(NON_CPG_CONTEXT) ||
                        (!estimateBisulfiteConversionRateFrom.equals(NON_CPG_CONTEXT) && estimateBisulfiteConversionRateFrom.equals(currentMultiPileup.getSequenceName()))) {

                    controlSampleEntryCounter++;

                    MethylationProfileRecord profileRecordForwardStrand = generateProfileRecord(previousMultiPileup, currentMultiPileup, Strand.FORWARD, sampleNames, null);
                    MethylationProfileRecord profileRecordReverseStrand = generateProfileRecord(previousMultiPileup, currentMultiPileup, Strand.REVERSE, sampleNames, null);

                    if (updateCTCounts(profileRecordForwardStrand, cCount, tCount, estimateBisulfiteConversionRateFrom)) {
                        usableLineCounter++;
                    }

                    if (updateCTCounts(profileRecordReverseStrand, cCount, tCount, estimateBisulfiteConversionRateFrom)) {
                        usableLineCounter++;
                    }

                }

            }

            previousMultiPileup = currentMultiPileup;

        }

        logger.info(lineCounter + " lines processed");
        if (!estimateBisulfiteConversionRateFrom.equals(NON_CPG_CONTEXT) && controlSampleEntryCounter == 0) {
            logger.warn("BAM file did not contain any pileup records for control sample: " + estimateBisulfiteConversionRateFrom);
        }
        logger.info("conversion rate estimated from " + usableLineCounter + " cytosin positions");

        return calculateBisfulfiteConversionRate(cCount, tCount);

    }


    private boolean updateCTCounts(MethylationProfileRecord profileRecord, Map<String, Long> cCount, Map<String, Long> tCount, String estimateBisulfiteConversionRateFrom) {

        boolean retVal = false;

        if(estimateBisulfiteConversionRateFrom == null){
            return retVal;
        }

        if (profileRecord != null &&
                !profileRecord.isRepeatMasked()) {

            for (MethylationProfileSampleRecord sampleRecord : profileRecord.getSampleRecords()) {

                if ((estimateBisulfiteConversionRateFrom.equals(NON_CPG_CONTEXT) && c_atc.matcher(sampleRecord.getSampleContext()).matches()
                        && sampleRecord.getSnpType() == 0 && !sampleRecord.hasIndels())
                        ||
                        (!estimateBisulfiteConversionRateFrom.equals(NON_CPG_CONTEXT) && profileRecord.getSequenceName().equals(estimateBisulfiteConversionRateFrom))) {

                    retVal = true;

                    long currentCcount = cCount.get(sampleRecord.getSampleName());
                    currentCcount = currentCcount + sampleRecord.getCCount();
                    cCount.put(sampleRecord.getSampleName(), currentCcount);

                    long currentTcount = tCount.get(sampleRecord.getSampleName());
                    currentTcount = currentTcount + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                    tCount.put(sampleRecord.getSampleName(), currentTcount);

                }

            }

        }

        return retVal;

    }

    /**
     * @param previousMultiPileup
     * @param currentMultiPileup
     * @param strand
     * @param sampleNames
     * @return
     */
    private MethylationProfileRecord generateProfileRecord(MultiPileup previousMultiPileup, MultiPileup currentMultiPileup, Strand strand, List<String> sampleNames, String estimateBisulfiteConversionRateFrom) {

        MethylationProfileRecord retVal = null;

        //check if the reference or any of the samples
        //has a cytosin at this position
        boolean isCytosinDiNucleotide = isCytosinDiNucleotide(previousMultiPileup, currentMultiPileup, strand);

        if (isCytosinDiNucleotide) {

            String referenceContext = getReferenceContext(previousMultiPileup, currentMultiPileup, strand);

            boolean isRepeatMasked = (previousMultiPileup.isRepeatMasked() || currentMultiPileup.isRepeatMasked());

            String sequenceName = previousMultiPileup.getSequenceName();
            int positionStart = previousMultiPileup.getPosition();
            int positionEnd = currentMultiPileup.getPosition();

            retVal = new MethylationProfileRecordImpl(sequenceName, positionStart, positionEnd, strand, referenceContext, isRepeatMasked);

            //iterate over sample pileups

            Iterator<Pileup> previousPileups = previousMultiPileup.getPileups().iterator();
            Iterator<Pileup> currentPileups = currentMultiPileup.getPileups().iterator();

            int sampleCounter = 0;
            boolean hasCoverage = false;
            while (previousPileups.hasNext()) {

                String sampleName = sampleNames.get(sampleCounter);
                sampleCounter++;

                Pileup previousReverseStrandPileup = previousPileups.next();
                Pileup previousForwardStrandPileup = previousPileups.next();

                Pileup currentReverseStrandPileup = currentPileups.next();
                Pileup currentForwardStrandPileup = currentPileups.next();

                String sampleContext = getSampleContext(previousForwardStrandPileup,
                        currentForwardStrandPileup,
                        previousReverseStrandPileup,
                        currentReverseStrandPileup,
                        strand);

                //check for SNPs on the reverse strand
                int snpType = determineSNPType(referenceContext, sampleContext);

                boolean hasIndels = false;

                int cCount = 0;
                int tCount = 0;
                int totalCoverage = 0;

                if (strand == Strand.FORWARD) {

                    cCount = previousForwardStrandPileup.getBaseCallCount("C");
                    tCount = previousForwardStrandPileup.getBaseCallCount("T");
                    totalCoverage = previousForwardStrandPileup.getCoverage();

                    if (previousForwardStrandPileup.hasIndels()) {
                        hasIndels = true;
                    }

                }

                //in reverse strand alignment
                //Cs on the reverse strand equal to Gs in the forward reference strand
                //and Ts are equal to As on the forward reference strand
                else if (strand == Strand.REVERSE) {

                    cCount = currentReverseStrandPileup.getBaseCallCount("G");
                    tCount = currentReverseStrandPileup.getBaseCallCount("A");
                    totalCoverage = currentReverseStrandPileup.getCoverage();

                    if (previousReverseStrandPileup.hasIndels()) {
                        hasIndels = true;
                    }

                }

                if (totalCoverage > 0) {
                    hasCoverage = true;
                }

                int ctCoverage = cCount + tCount;
                int nonCtCoverage = totalCoverage - ctCoverage;


                MethylationProfileSampleRecordImpl sampleRecord = new MethylationProfileSampleRecordImpl(retVal, sampleContext, cCount, ctCoverage, nonCtCoverage);
                sampleRecord.setSnpType(snpType);
                sampleRecord.setHasIndels(hasIndels);
                sampleRecord.setSampleName(sampleName);
                retVal.addSampleRecord(sampleRecord);

                //update non-CpG C and T counts
                if (estimateBisulfiteConversionRateFrom != null && !retVal.getSequenceName().equals(estimateBisulfiteConversionRateFrom) &&
                        c_atc.matcher(sampleContext).matches() && sampleRecord.getSnpType() == 0 &&
                        !sampleRecord.hasIndels() && !retVal.isRepeatMasked()
                        ) {

                    long currentCcount = this.nonCpGContextCCount.get(sampleContext).get(sampleName);
                    currentCcount = currentCcount + cCount;
                    this.nonCpGContextCCount.get(sampleContext).put(sampleName, currentCcount);

                    long currentTotalCcount = this.nonCpGContextCCount.get("total").get(sampleName);
                    currentTotalCcount = currentTotalCcount + cCount;
                    this.nonCpGContextCCount.get("total").put(sampleName, currentTotalCcount);

                    long currentTcount = this.nonCpGContextTCount.get(sampleContext).get(sampleName);
                    currentTcount = currentTcount + tCount;
                    this.nonCpGContextTCount.get(sampleContext).put(sampleName, currentTcount);

                    long currentTotalTcount = this.nonCpGContextTCount.get("total").get(sampleName);
                    currentTotalTcount = currentTotalTcount + tCount;
                    this.nonCpGContextTCount.get("total").put(sampleName, currentTotalTcount);

                }

            }

            //if cytosin has no coverage in any
            //of the samples on the respective
            //strand
            if (!hasCoverage) {
                retVal = null;
            }

        }

        return retVal;

    }

    /**
     * @param previousMultiPileup
     * @param currentMultiPileup
     * @param strand
     * @return
     */
    private boolean isCytosinDiNucleotide(MultiPileup previousMultiPileup, MultiPileup currentMultiPileup, Strand strand) {

        boolean retVal = false;

        if (strand == Strand.FORWARD) {

            //check reference first
            if (previousMultiPileup.getReferenceBase().equalsIgnoreCase("C")) {
                retVal = true;
            }

            //and then sample
            else {

                Iterator<Pileup> previousPileups = previousMultiPileup.getPileups().iterator();
                while (previousPileups.hasNext()) {

                    //get reverse strand pileup
                    Pileup previousReverseStrandPileup = previousPileups.next();

                    //skip forward strand pileup
                    previousPileups.next();

                    //check on reverse strand if we are dealing with a cytosin
                    if (previousReverseStrandPileup.getConsensusBaseCall().contains("C")) {
                        retVal = true;
                    }

                }

            }

        } else if (strand == Strand.REVERSE) {

            //check reference first (G on the reference forward strand
            //means C on the reverse strand)
            if (currentMultiPileup.getReferenceBase().equalsIgnoreCase("G")) {
                retVal = true;
            }

            //and then sample
            else {                                                                             

                Iterator<Pileup> currentPileups = currentMultiPileup.getPileups().iterator();
                while (currentPileups.hasNext()) {

                    //skip reverse strand pileup
                    currentPileups.next();

                    //get forward strand pileup
                    Pileup currentForwardStrandPileup = currentPileups.next();

                    //check on forward strand if we are dealing with a C
                    //on the reverse strand, i.e. a G on the forward strand
                    if (currentForwardStrandPileup.getConsensusBaseCall().contains("G")) {
                        retVal = true;
                    }

                }

            }

        }

        return retVal;

    }

    /**
     * Returns the sample context of a cytosin di-nucleotide position, i.e. Cp[ACTGN]. Uses base calls
     * from the oposite genomic template strand to determine the genotype unless there is a G at the
     * respective position. In that case the oposite strand will have a cytosin position which is subject
     * to bisulfite conversion and cannot be used to reliably predict the genotype.
     *
     * @param previousForwardPileup the first forward strand pileup of the sample in the mpileup file
     * @param currentForwardPileup  the second forward strand pileup of the sample in the mpileup file
     * @param previousReversePileup the first reverse strand pileup of the sample in the mpileup file
     * @param currentReversePileup  the second reverse strand pileup of the sample in the mpileup file
     * @param strand                the strand for which to return the sample context
     * @return the CpN di-nucleotide context in the sample
     */
    private String getSampleContext(Pileup previousForwardPileup,
                                    Pileup currentForwardPileup,
                                    Pileup previousReversePileup,
                                    Pileup currentReversePileup,
                                    Strand strand) {

        StringBuffer retVal = new StringBuffer();

        String previousSampleBaseReverseStrand = previousReversePileup.getConsensusBaseCall();
        String currentSampleBaseReverseStrand = currentReversePileup.getConsensusBaseCall();

        String previousSampleBaseForwardStrand = previousForwardPileup.getConsensusBaseCall();
        String currentSampleBaseForwardStrand = currentForwardPileup.getConsensusBaseCall();

        //check cases where first or second position of di-nucleotide
        //are subject to methylation on either forward or reverse strand
        if (strand == Strand.FORWARD) {


            //if there is a G on the forward strand
            //the reverse strand will have a cytosin position
            //which might be methylated and can thus not be used
            //for SNP determination and we need to use the forward strand
            String previousSampleBase = previousSampleBaseReverseStrand;
            if (previousSampleBaseForwardStrand.contains("G")) {
                previousSampleBase = previousSampleBaseForwardStrand;
            }

            String currentSampleBase = currentSampleBaseReverseStrand;
            if (currentSampleBaseForwardStrand.contains("G")) {
                currentSampleBase = currentSampleBaseForwardStrand;
            }

            retVal.append(previousSampleBase).append(currentSampleBase);

        } else if (strand == Strand.REVERSE) {

            //if there is a G on the reverse strand (called as a C on the reference)
            //the forward strand will have a cytosin position
            //which might be methylated and can thus not be used
            //for SNP determination and we need to use reverse strand
            String previousSampleBase = previousSampleBaseForwardStrand;
            if (previousSampleBaseReverseStrand.contains("C")) {
                previousSampleBase = previousSampleBaseReverseStrand;
            }

            String currentSampleBase = currentSampleBaseForwardStrand;
            if (currentSampleBaseReverseStrand.contains("C")) {
                currentSampleBase = currentSampleBaseReverseStrand;
            }

            retVal.append(complement(currentSampleBase)).append(complement(previousSampleBase));

        } else {
            throw new IllegalArgumentException("Strand must be forward or reverse: " + strand);
        }

        return retVal.toString().toUpperCase();

    }

    /**
     * Returns the reference context of a cytosin di-nucleotide position, i.e. Cp[ACTGN].
     *
     * @param previousPileup the first multi-pileup record in the mpileup file
     * @param currentPileup  the second multi-pileup record in the mpileup file
     * @param strand         the stand for which the reference context is returned
     * @return the CpN di-nucleotide context
     */
    public String getReferenceContext(MultiPileup previousPileup, MultiPileup currentPileup, Strand strand) {

        StringBuffer retVal = new StringBuffer();

        String previousReferenceBase = previousPileup.getReferenceBase().toUpperCase();
        String currentReferenceBase = currentPileup.getReferenceBase().toUpperCase();

        if (strand == Strand.FORWARD) {
            retVal.append(previousReferenceBase).append(currentReferenceBase);
        } else if (strand == Strand.REVERSE) {
            retVal.append(complement(currentReferenceBase)).append(complement(previousReferenceBase));
        } else {                                                                               
            throw new IllegalArgumentException("Strand must be forward or reverse: " + strand);
        }

        return retVal.toString().toUpperCase();

    }


    /**
     * Pattern matching for SNP determination
     * <p/>
     * SNP
     * 0=no SNP
     * CpG context (CpG in reference)
     * 1=CG -> TG transition                    loss of cytosin -> false positive unmethylated
     * 2=CG -> [AG]G transversion               loss of cytosin -> no methylation signal
     * 3=CG -> C[ATC] transition/transversion   loss of CpG di-nucleotide -> potentially methylated site
     * 4=CG -> [ATG][ATC]                       loss of CpG di-nucleotide -> no methylation signal
     * 5=TG -> CG transition                    gain of CpG di-nucleotide -> potential CpG methylation site
     * 6=[AG]G -> CG transversion               gain of CpG di-nucleotide -> potential CpG methylation site
     * 7=CA -> CG  transition                   gain of CpG di-nucleotide -> potential CpG methylation site
     * 8=C[CT] -> CG transversion               gain of CpG di-nucleotide -> potential CpG methylation site
     * 9=[ATG][ATC] -> CG                       gain of CpG di-nucleotide -> potential CpG methylation site
     * non-CpG context (non-CpG in reference)
     * 10=C[ACT] -> T[ACT] transition            loss of cytosin -> false positive unmethylated
     * 11=C[ACT] -> [AG][ACT] transversion       loss of cytosin -> no methylation signal
     * 12=T[ACT] -> C[ACT] transition            gain of cytosin -> potential methylation site
     * 13=[AG][ACT] -> C[ACT] transversion       gain of cytosin -> potential methylation site
     * 14=other
     *
     * @param referenceContext
     * @param sampleContext
     * @return
     */
    private int determineSNPType(String referenceContext, String sampleContext) {


//        Pattern cg = Pattern.compile("CG");
//        Pattern tg = Pattern.compile("TG");
//        Pattern ag_g = Pattern.compile("[AG]G");
//        Pattern c_atc = Pattern.compile("C[ATC]");
//        Pattern ca = Pattern.compile("CA");
//        Pattern c_ct = Pattern.compile("C[CT]");
//        Pattern c_act = Pattern.compile("C[ACT]");
//        Pattern t_act = Pattern.compile("T[ACT]");
//        Pattern ag_act = Pattern.compile("AG[ACT]");

        if (referenceContext.equals(sampleContext)) {
            return 0;
        }
        //if one of the sample bases used for SNP determination is not called
        //we can return -1
        if (sampleContext.contains("N")) {
            return -1;
        }

        //if one of the sample base calls is heterozygous we return -1
        if (sampleContext.contains("/")) {
            return -1;
        }

        int retVal = 0;

        //CpG context (CpG in reference)

        //1 = CG -> TG transition
        if (cg.matcher(referenceContext).matches() && tg.matcher(sampleContext).matches()) {
            retVal = 1;
            logger.debug("1 = CG -> TG transition");
        }

        //2 = CG -> [AG]G transversion
        else if (cg.matcher(referenceContext).matches() && ag_g.matcher(sampleContext).matches()) {
            retVal = 2;
            logger.debug("2 = CG -> [AG]G transversion");
        }

        //3 = CG -> C[ATC] transition/transversion
        else if (cg.matcher(referenceContext).matches() && c_atc.matcher(sampleContext).matches()) {
            retVal = 3;
            logger.debug("3 = CG -> C[ATC] transition/transversion");
        }

        //4 = CG -> [ATG][ATC]
        else if (cg.matcher(referenceContext).matches() && atg_atc.matcher(sampleContext).matches()) {
            retVal = 4;
            logger.debug("4 = CG -> [ATG][ATC]");
        }

        //5 = TG -> CG transition
        else if (tg.matcher(referenceContext).matches() && cg.matcher(sampleContext).matches()) {
            retVal = 5;
            logger.debug("5 = TG -> CG transition");
        }

        //6 = [AG]G -> CG transversion
        else if (ag_g.matcher(referenceContext).matches() && cg.matcher(sampleContext).matches()) {
            retVal = 6;
            logger.debug("6 = [AG]G -> CG transversion");
        }

        //7 = CA -> CG  transition
        else if (ca.matcher(referenceContext).matches() && cg.matcher(sampleContext).matches()) {
            retVal = 7;
            logger.debug("7 = CA -> CG  transition");
        }

        //8 = C[CT] -> CG transversion
        else if (c_ct.matcher(referenceContext).matches() && cg.matcher(sampleContext).matches()) {
            retVal = 8;
            logger.debug("8 = C[CT] -> CG transversion");
        }


        //9 = [ATG][ATC] -> CG
        else if (atg_atc.matcher(referenceContext).matches() && cg.matcher(sampleContext).matches()) {
            retVal = 9;
            logger.debug("9 = [ATG][ATC] -> CG");
        }


        //non-CpG context (non-CpG in reference)

        //10 = C[ACT] -> T[ACT] transition
        else if (c_act.matcher(referenceContext).matches() && t_act.matcher(sampleContext).matches()) {
            retVal = 10;
            logger.debug("10 = C[ACT] -> T[ACT] transition");
        }

        //11 = C[ACT] -> [AG][ACT] transversion
        else if (c_act.matcher(referenceContext).matches() && ag_act.matcher(sampleContext).matches()) {
            retVal = 11;
            logger.debug("11 = C[ACT] -> [AG][ACT] transversion");
        }

        //12 = T[ACT] -> C[ACT] transition
        else if (t_act.matcher(referenceContext).matches() && c_act.matcher(sampleContext).matches()) {
            retVal = 12;
            logger.debug("12 = T[ACT] -> C[ACT] transition");
        }

        //13=[AG][ACT] -> C[ACT] transversion
        else if (ag_act.matcher(referenceContext).matches() && c_act.matcher(sampleContext).matches()) {
            retVal = 13;
            logger.debug("13=[AG][ACT] -> C[ACT] transversion");
        } else {
            retVal = 14;
        }

        return retVal;

    }

    /**
     * @param nucleotide
     * @return
     */
    private String complement(String nucleotide) {

        if (nucleotide.toUpperCase().equals("G")) {
            return "C";
        } else if (nucleotide.toUpperCase().equals("A")) {
            return "T";
        } else if (nucleotide.toUpperCase().equals("T")) {
            return "A";
        } else if (nucleotide.toUpperCase().equals("C")) {
            return "G";
        } else if (nucleotide.contains("/")) {

            String[] nucleotides = nucleotide.split("/");
            String firstComplement = complement(nucleotides[0]);
            String secondComplement = complement(nucleotides[1]);

            return firstComplement + "/" + secondComplement;

        } else {
            return nucleotide;
        }

    }

}
