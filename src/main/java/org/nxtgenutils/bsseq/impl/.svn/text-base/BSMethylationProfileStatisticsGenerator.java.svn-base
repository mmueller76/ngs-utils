package org.nxtgenutils.bsseq.impl;

import org.nxtgenutils.io.MethylationProfileParser;
import org.nxtgenutils.io.MethylationProfileRecord;
import org.nxtgenutils.io.MethylationProfileSampleRecord;
import org.nxtgenutils.io.impl.SimpleMethylationProfileParser;
import org.nxtgenutils.Strand;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.text.NumberFormat;
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
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 21-Nov-2011
 * Time: 13:28:28
 * To change this template use File | Settings | File Templates.
 */
public class BSMethylationProfileStatisticsGenerator {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(BSMethylationProfileStatisticsGenerator.class);


    private int coverageCutOff = 1000;

    private Map<String, int[]> cummulativeCpGCoverageBySample = new HashMap<String, int[]>();
    private Map<String, int[][]> cummulativeCpGCoverageAcrossSamples = new HashMap<String, int[][]>();
    private Map<String, int[]> cummulativeNonCpGCoverageBySample = new HashMap<String, int[]>();
    private Map<String, int[][]> cummulativeNonCpGCoverageAcrossSamples = new HashMap<String, int[][]>();
    private Map<String, Map<Integer, List<String>>> group2Coverage2SampleNames = new HashMap<String, Map<Integer, List<String>>>();

    private Map<String, Long> controlCCount = new HashMap<String, Long>();
    private Map<String, Long> controlTCount = new HashMap<String, Long>();

    private Map<String, Integer> controlConversionRatePositionCount = new HashMap<String, Integer>();
    private Map<String, Map<String, Integer>> nonCpGConversionRatePositionCount = new HashMap<String, Map<String, Integer>>();

    private Map<String, Map<String, Long>> nonCpGContextCCount = new HashMap<String, Map<String, Long>>();
    private Map<String, Map<String, Long>> nonCpGContextTCount = new HashMap<String, Map<String, Long>>();

    private Set<String> nonCpGContexts = new HashSet<String>();

    static Pattern c_atc = Pattern.compile("C[ATC]");

    public BSMethylationProfileStatisticsGenerator(File methylationProfileFile,
                                                   File sampleNamesFile,
                                                   File sampleGroupingFile,
                                                   String controleSequenceName
    ) {

        List<String> sampleNames = readSampleNames(sampleNamesFile);
        Map<String, List<String>> group2Sample = readSampleGrouping(sampleGroupingFile);
        Map<String, String> sample2Group = new LinkedHashMap<String, String>();
        for (String group : group2Sample.keySet()) {
            for (String sample : group2Sample.get(group)) {
                sample2Group.put(sample, group);
            }
        }

        String outputPrefix = methylationProfileFile.getAbsolutePath();

        initialiseMaps(group2Sample, controleSequenceName);

        MethylationProfileParser methylationProfileInput = new SimpleMethylationProfileParser(methylationProfileFile, sampleNames);
        methylationProfileInput.setSampleNames(sampleNames);

        Iterator<MethylationProfileRecord> profileRecordIterator = methylationProfileInput.iterator();

        try {

            while (profileRecordIterator.hasNext()) {

                MethylationProfileRecord record = profileRecordIterator.next();
                updateCounts(record, sample2Group);
                updateCTCounts(record, controleSequenceName);

            }

        } catch (Exception e) {

            logger.error("Exception while parsing methylation profile.", e);
            printStats(outputPrefix);

        }

        printStats(outputPrefix);

    }

    private void printStats(String outputPrefix) {


        //cummulative coverage by sample

        for (String sampleAndStrand : cummulativeCpGCoverageAcrossSamples.keySet()) {

            try {

                PrintWriter pw = new PrintWriter(outputPrefix + "cum_cpg_cov_by_sample." + sampleAndStrand + ".tsv");

                int[] coverage = cummulativeCpGCoverageBySample.get(sampleAndStrand);
                for (int i = 0; i < coverage.length; i++) {

                    String rowName = "" + (i + 1);
                    if (i == coverageCutOff) {
                        rowName = ">" + coverageCutOff;
                    }

                    //calculate cummulative coverage
                    int cummulativeFrequency = 0;
                    for (int i2 = i; i2 < coverage.length; i2++) {
                        cummulativeFrequency = cummulativeFrequency + coverage[i2];
                    }

                    pw.print(rowName + "\t" + cummulativeFrequency);
                    pw.flush();

                }

                pw.close();

            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }

        //cummulative non-CpG coverage
        for (String sampleAndStrand : cummulativeNonCpGCoverageBySample.keySet()) {

            try {

                PrintWriter pw = new PrintWriter(outputPrefix + "cum_non_cpg_cov_by_sample." + sampleAndStrand + ".tsv");

                int[] coverage = cummulativeNonCpGCoverageBySample.get(sampleAndStrand);
                for (int i = 0; i < coverage.length; i++) {

                    String rowName = "" + (i + 1);
                    if (i == coverageCutOff) {
                        rowName = ">" + coverageCutOff;
                    }

                    //calculate cummulative coverage
                    int cummulativeFrequency = 0;
                    for (int i2 = i; i2 < coverage.length; i2++) {
                        cummulativeFrequency = cummulativeFrequency + coverage[i2];
                    }

                    pw.print(rowName + "\t" + cummulativeFrequency);
                    pw.flush();

                }

                pw.close();

            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }

        //cummulative CpG coverage across samples
        for (String groupAndStrand : cummulativeCpGCoverageAcrossSamples.keySet()) {

            try {

                PrintWriter pw = new PrintWriter(outputPrefix + "cum_cpg_cov_across_samples." + groupAndStrand + ".tsv");

                int[][] coverage = cummulativeCpGCoverageAcrossSamples.get(groupAndStrand);

                //print column headers
                pw.print("cov");
                for (int c = 0; c < coverage[0].length; c++) {

                    pw.print("\t" + (c + 1) + "_sample");

                }

                pw.println();

                for (int i = 0; i < coverage.length; i++) {

                    String rowName = "" + (i + 1);
                    if (i == coverageCutOff) {
                        rowName = ">" + coverageCutOff;
                    }

                    pw.print(rowName);

                    for (int j = 0; j < coverage[i].length; j++) {

                        //calculate cummulative coverage
                        int cummulativeFrequency = 0;

                        //sum frequencies
                        for (int i2 = i; i2 < coverage.length; i2++) {
                            for (int j2 = i; j2 < coverage.length; j2++) {
                                cummulativeFrequency = cummulativeFrequency + coverage[i2][j2];
                            }
                        }
                        pw.print("\t" + cummulativeFrequency);

                    }

                    pw.println();
                    pw.flush();

                }

                pw.close();

            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }


        //cummulative non-CpG coverage across samples
        for (String groupAndStrand : cummulativeNonCpGCoverageAcrossSamples.keySet()) {

            try {

                PrintWriter pw = new PrintWriter(outputPrefix + "cum_non_cpg_cov_across_samples." + groupAndStrand + ".tsv");

                int[][] coverage = cummulativeNonCpGCoverageAcrossSamples.get(groupAndStrand);

                //print column headers
                pw.print("cov");
                for (int c = 0; c < coverage[0].length; c++) {

                    pw.print("\t" + (c + 1) + "_sample");

                }

                pw.println();

                for (int i = 0; i < coverage.length; i++) {

                    String rowName = "" + (i + 1);
                    if (i == coverageCutOff) {
                        rowName = ">" + coverageCutOff;
                    }

                    pw.print(rowName);

                    for (int j = 0; j < coverage[i].length; j++) {

                        //calculate cummulative coverage
                        int cummulativeFrequency = 0;

                        //sum frequencies
                        for (int i2 = i; i2 < coverage.length; i2++) {
                            for (int j2 = i; j2 < coverage.length; j2++) {
                                cummulativeFrequency = cummulativeFrequency + coverage[i2][j2];
                            }
                        }
                        pw.print("\t" + cummulativeFrequency);

                    }

                    pw.println();
                    pw.flush();

                }

                pw.close();

            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }

        //controle conversion rate
        NumberFormat formatter = new DecimalFormat("#.#");

        Map<String, Double> conversionRates = calculateBisfulfiteConversionRate(controlCCount, controlTCount);
        try {

            PrintWriter pw = new PrintWriter(outputPrefix + "controle_conversion_rates.tsv");


            pw.println("sample\trate\tpositions");
            for (String sampleName : conversionRates.keySet()) {

                int count = controlConversionRatePositionCount.get(sampleName);
                double rate = conversionRates.get(sampleName);
                pw.println(sampleName + "\t" + formatter.format(rate * 100.0) + "\t" + count);
                pw.flush();

            }

            pw.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


        try {

            PrintWriter pw = new PrintWriter(outputPrefix + "non_cpg_conversion_rates.tsv");
            pw.println("sample\tcontext\trate\tpositions");

            for (String context : nonCpGContextCCount.keySet()) {

                Map<String, Long> cCount = nonCpGContextCCount.get(context);
                Map<String, Long> tCount = nonCpGContextCCount.get(context);

                conversionRates = calculateBisfulfiteConversionRate(cCount, tCount);
                Map<String, Integer> positionCounts = nonCpGConversionRatePositionCount.get(context);

                for (String sampleName : conversionRates.keySet()) {

                    int count = positionCounts.get(sampleName);
                    double rate = conversionRates.get(sampleName);
                    pw.println(sampleName + "\t" + context + "\t" + formatter.format(rate * 100.0) + "\t" + count);
                    pw.flush();

                }

            }

            pw.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }


    }

    /**
     * @param sampleCCount
     * @param sampleTCount
     * @return
     */
    private Map<String, Double> calculateBisfulfiteConversionRate(Map<String, Long> sampleCCount, Map<String, Long> sampleTCount) {

        Map<String, Double> retVal = new HashMap<String, Double>();

        for (String name : sampleCCount.keySet()) {
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

    private void updateCounts(MethylationProfileRecord record, Map<String, String> sample2Group) {

        String referenceContext = record.getReferenceContext();
        Strand strand = record.getStrand();

        for (MethylationProfileSampleRecord sampleRecord : record.getSampleRecords()) {

            String sampleName = sampleRecord.getSampleName();
            String groupName = sample2Group.get(sampleName);
            int coverage = sampleRecord.getCtCoverage();
            if (coverage > 1000) {
                coverage = 1001;
            }

            if (!group2Coverage2SampleNames.get(groupName).containsKey(coverage)) {
                group2Coverage2SampleNames.get(groupName).put(coverage, new ArrayList<String>());
            }
            group2Coverage2SampleNames.get(groupName).get(coverage).add(sampleName);


            if (referenceContext.equals("CG")) {
                cummulativeCpGCoverageBySample.get(sampleName + "_" + strand)[coverage - 1]++;
            } else {
                cummulativeNonCpGCoverageBySample.get(sampleName + "_" + strand)[coverage - 1]++;
            }

        }

        for (String group : group2Coverage2SampleNames.keySet()) {

            for (int coverage : group2Coverage2SampleNames.get(group).keySet()) {

                int sampleCount = group2Coverage2SampleNames.get(group).get(coverage).size();
                if (coverage > 1000) {
                    coverage = 1001;
                }

                if (referenceContext.equals("CG")) {
                    cummulativeCpGCoverageAcrossSamples.get(group + "_" + strand)[coverage - 1][sampleCount]++;
                } else {
                    cummulativeNonCpGCoverageAcrossSamples.get(group + "_" + strand)[coverage - 1][sampleCount]++;
                }

            }

            group2Coverage2SampleNames.get(group).clear();

        }

    }


    private boolean updateCTCounts(MethylationProfileRecord profileRecord, String controlSequence) {

        boolean retVal = false;

        if (profileRecord == null ||
                profileRecord.isRepeatMasked()) {
            return false;
        }


        for (MethylationProfileSampleRecord sampleRecord : profileRecord.getSampleRecords()) {

            if (controlSequence != null && profileRecord.getSequenceName().equals(controlSequence)) {

                retVal = true;

                String sampleName = sampleRecord.getSampleName();

                long currentCcount = controlCCount.get(sampleName);
                currentCcount = currentCcount + sampleRecord.getCCount();
                controlCCount.put(sampleName, currentCcount);

                long currentTcount = controlTCount.get(sampleName);
                currentTcount = currentTcount + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                controlTCount.put(sampleName, currentTcount);

                int currentPositionCount = controlConversionRatePositionCount.get(sampleName);
                currentPositionCount++;
                controlConversionRatePositionCount.put(sampleName, currentPositionCount);

            } else if (sampleRecord.getSnpType() == 0 && !sampleRecord.hasIndels() &&
                    c_atc.matcher(sampleRecord.getSampleContext()).matches()) {

                retVal = true;

                String sampleName = sampleRecord.getSampleName();

                long currentCcount = nonCpGContextCCount.get(sampleRecord.getSampleContext()).get(sampleRecord.getSampleName());
                currentCcount = currentCcount + sampleRecord.getCCount();
                nonCpGContextCCount.get(sampleRecord.getSampleContext()).put(sampleRecord.getSampleName(), currentCcount);

                long currentTcount = nonCpGContextTCount.get(sampleRecord.getSampleContext()).get(sampleRecord.getSampleName());
                currentTcount = currentTcount + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                nonCpGContextTCount.get(sampleRecord.getSampleContext()).put(sampleRecord.getSampleName(), currentTcount);

                long currentCcountTotal = nonCpGContextCCount.get("total").get(sampleRecord.getSampleName());
                currentCcountTotal = currentCcountTotal + sampleRecord.getCCount();
                nonCpGContextCCount.get("total").put(sampleRecord.getSampleName(), currentCcountTotal);

                long currentTcountTotal = nonCpGContextTCount.get("total").get(sampleRecord.getSampleName());
                currentTcountTotal = currentTcountTotal + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                nonCpGContextTCount.get("total").put(sampleRecord.getSampleName(), currentTcountTotal);

                int currentPositionCount = nonCpGConversionRatePositionCount.get(sampleRecord.getSampleContext()).get(sampleName);
                currentPositionCount++;
                nonCpGConversionRatePositionCount.get(sampleRecord.getSampleContext()).put(sampleName, currentPositionCount);

                int currentPositionCountTotal = nonCpGConversionRatePositionCount.get("total").get(sampleName);
                currentPositionCountTotal++;
                nonCpGConversionRatePositionCount.get("total").put(sampleName, currentPositionCountTotal);

            }

        }


        return retVal;

    }

    private void initialiseMaps(Map<String, List<String>> sampleGrouping, String controleSequenceName) {

        //non-CpG context
        nonCpGContexts.add("CT");
        nonCpGContexts.add("CA");
        nonCpGContexts.add("CC");
        nonCpGContexts.add("total");

        for (String context : nonCpGContexts) {

            nonCpGContextCCount.put(context, new HashMap<String, Long>());
            nonCpGContextTCount.put(context, new HashMap<String, Long>());
            nonCpGConversionRatePositionCount.put(context, new HashMap<String, Integer>());
        }

        for (String group : sampleGrouping.keySet()) {

            int sampleCount = sampleGrouping.get(group).size();

            //cummulative CpG coverage across samples
            cummulativeCpGCoverageAcrossSamples.put(group + "_" + Strand.FORWARD, new int[coverageCutOff + 1][sampleCount]);
            cummulativeCpGCoverageAcrossSamples.put(group + "_" + Strand.REVERSE, new int[coverageCutOff + 1][sampleCount]);

            //cummulative non-CpG coverage across samples
            cummulativeNonCpGCoverageAcrossSamples.put(group + "_" + Strand.FORWARD, new int[coverageCutOff + 1][sampleCount]);
            cummulativeNonCpGCoverageAcrossSamples.put(group + "_" + Strand.REVERSE, new int[coverageCutOff + 1][sampleCount]);

            group2Coverage2SampleNames.put(group, new HashMap<Integer, List<String>>());

            for (String sample : sampleGrouping.get(group)) {

                //cummulative CpG coverage by sample
                cummulativeCpGCoverageBySample.put(sample + "_" + Strand.FORWARD, new int[coverageCutOff + 1]);
                cummulativeCpGCoverageBySample.put(sample + "_" + Strand.REVERSE, new int[coverageCutOff + 1]);

                //cummulative non-CpG coverage by sample
                cummulativeNonCpGCoverageBySample.put(sample + "_" + Strand.FORWARD, new int[coverageCutOff + 1]);
                cummulativeNonCpGCoverageBySample.put(sample + "_" + Strand.REVERSE, new int[coverageCutOff + 1]);

                //CT count
                //control sequence
                controlCCount.put(sample, (long) 0);
                controlTCount.put(sample, (long) 0);


                //non-CpG context
                for (String context : nonCpGContexts) {

                    nonCpGContextCCount.get(context).put(sample, (long) 0);
                    nonCpGContextTCount.get(context).put(sample, (long) 0);
                    nonCpGConversionRatePositionCount.get(context).put(sample, 0);
                }

                controlConversionRatePositionCount.put(sample, 0);

            }

        }

    }

//    private Map<String, BufferedWriter> intitialiseWriters(String outputPrefix){
//
//         Map<String, BufferedWriter>
//
//
//    }

    private List<String> readSampleNames(File sampleNames) {

        List<String> retVal = new ArrayList<String>();

        try {

            BufferedReader br = new BufferedReader(new FileReader(sampleNames));
            String line;
            while ((line = br.readLine()) != null) {
                line = line.replaceAll(" ", "");
                retVal.add(line);
            }
            br.close();

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        return retVal;

    }

    private Map<String, List<String>> readSampleGrouping(File sampleGrouping) {

        Map<String, List<String>> retVal = new LinkedHashMap<String, List<String>>();

        try {

            BufferedReader br = new BufferedReader(new FileReader(sampleGrouping));
            String line;

            while ((line = br.readLine()) != null) {
                line = line.replaceAll(" ", "");
                String[] sample2Group = line.split("\t");
                String group = sample2Group[1];
                String sample = sample2Group[0];

                if (!retVal.containsKey(group)) {
                    retVal.put(group, new ArrayList<String>());
                }

                retVal.get(group).add(sample);
            }
            br.close();

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        return retVal;

    }


    public static void main(String[] args) {

        String methylationProfilePath = "/home/mmuelle1/dev/java/nxtgen-utils/testdata/bn_shr.meth.profile.cg.chr1.profile.gz";
        String methylationProfileSampleNamesPath = "/home/mmuelle1/dev/java/nxtgen-utils/testdata/bn_shr.meth.profile.sampleNames.txt";
        String methylationProfileSampleGroupingPath = "/home/mmuelle1/dev/java/nxtgen-utils/testdata/bn_shr.meth.profile.sampleGrouping.tsv";
        String controlSequence = null;

        BSMethylationProfileStatisticsGenerator profileStats = new BSMethylationProfileStatisticsGenerator(new File(methylationProfilePath),
                new File(methylationProfileSampleNamesPath),
                new File(methylationProfileSampleGroupingPath),
                controlSequence);

    }

}
