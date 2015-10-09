package org.nxtgenutils.bsseq.impl;

import org.nxtgenutils.io.MethylationProfileParser;
import org.nxtgenutils.io.MethylationProfileRecord;
import org.nxtgenutils.io.MethylationProfileSampleRecord;
import org.nxtgenutils.io.impl.SimpleMethylationProfileParser;


import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.text.NumberFormat;
import java.text.DecimalFormat;

import org.apache.log4j.Logger;

/**
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 21-Nov-2011
 * Time: 21:36:42
 */

/**
 * @author Michael Mueller
 */
public class BSMethylationProfileConversionRateEstimator {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(BSMethylationProfileConversionRateEstimator.class);

    private static Pattern c_atc = Pattern.compile("C[ATC]");
    private static Set<String> contexts = new HashSet<String>();

    static {
        contexts.add("CT");
        contexts.add("CA");
        contexts.add("CC");
        contexts.add("CG");
        contexts.add("total");
        contexts.add("non_cpg_total");

    }

    private Map<String, Map<String, Long>> context2Sample2CCount = new HashMap<String, Map<String, Long>>();
    private Map<String, Map<String, Long>> context2Sample2TCount = new HashMap<String, Map<String, Long>>();
    private Map<String, Map<String, Long>> context2Sample2PositionCount = new HashMap<String, Map<String, Long>>();


    public BSMethylationProfileConversionRateEstimator(File methylationProfileFile,
                                                       File sampleNamesFile
    ) {

        List<String> sampleNames = readSampleNames(sampleNamesFile);
        String outputPrefix = methylationProfileFile.getAbsolutePath();

        initialiseMaps(sampleNames);

        MethylationProfileParser methylationProfileInput = new SimpleMethylationProfileParser(methylationProfileFile, sampleNames);

        Iterator<MethylationProfileRecord> profileRecordIterator = methylationProfileInput.iterator();

        try {

            while (profileRecordIterator.hasNext()) {

                MethylationProfileRecord record = profileRecordIterator.next();
                updateCTCounts(record);

            }

        } catch (Exception e) {

            logger.error("Exception while parsing methylation profile.", e);
            printStats(outputPrefix);

        }

        printStats(outputPrefix);

    }

    private void printStats(String outputPrefix) {

        try {

            NumberFormat formatter = new DecimalFormat("#.#");

            PrintWriter pw = new PrintWriter(outputPrefix + "non_cpg_conversion_rates.tsv");
            pw.println("sample\tcontext\trate\tpositions");

            for (String context : context2Sample2CCount.keySet()) {

                Map<String, Long> cCount = context2Sample2CCount.get(context);
                Map<String, Long> tCount = context2Sample2CCount.get(context);

                Map<String, Double> conversionRates = calculateBisfulfiteConversionRate(cCount, tCount);
                Map<String, Long> positionCounts = context2Sample2PositionCount.get(context);

                for (String sampleName : conversionRates.keySet()) {

                    long count = positionCounts.get(sampleName);
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


    private void initialiseMaps(List<String> sampleNames) {

        for (String context : contexts) {
            context2Sample2CCount.put(context, new HashMap<String, Long>());
            context2Sample2TCount.put(context, new HashMap<String, Long>());
            context2Sample2PositionCount.put(context, new HashMap<String, Long>());
        }

        for (String sample : sampleNames) {

            for (String context : contexts) {

                context2Sample2CCount.get(context).put(sample, (long) 0);
                context2Sample2TCount.get(context).put(sample, (long) 0);
                context2Sample2PositionCount.get(context).put(sample, (long) 0);

            }

        }

    }

    private boolean updateCTCounts(MethylationProfileRecord profileRecord) {

        boolean retVal = false;

        if (profileRecord == null ||
                profileRecord.isRepeatMasked()) {
            return false;
        }

        for (MethylationProfileSampleRecord sampleRecord : profileRecord.getSampleRecords()) {

            if (sampleRecord.getSnpType() == 0 && !sampleRecord.hasIndels()) {

                retVal = true;

                String sampleName = sampleRecord.getSampleName();

                long currentCcount = context2Sample2CCount.get(sampleRecord.getSampleContext()).get(sampleName);
                currentCcount = currentCcount + sampleRecord.getCCount();
                context2Sample2CCount.get(sampleRecord.getSampleContext()).put(sampleName, currentCcount);

                long currentTcount = context2Sample2TCount.get(sampleRecord.getSampleContext()).get(sampleName);
                currentTcount = currentTcount + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                context2Sample2TCount.get(sampleRecord.getSampleContext()).put(sampleName, currentTcount);

                long currentPositionCount = context2Sample2PositionCount.get(sampleRecord.getSampleContext()).get(sampleName);
                currentPositionCount++;
                context2Sample2PositionCount.get(sampleRecord.getSampleContext()).put(sampleName, currentPositionCount);


                long currentCcountTotal = context2Sample2CCount.get("total").get(sampleName);
                currentCcountTotal = currentCcountTotal + sampleRecord.getCCount();
                context2Sample2CCount.get("total").put(sampleName, currentCcountTotal);

                long currentTcountTotal = context2Sample2TCount.get("total").get(sampleName);
                currentTcountTotal = currentTcountTotal + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                context2Sample2TCount.get("total").put(sampleName, currentTcountTotal);

                long currentPositionCountTotal = context2Sample2PositionCount.get("total").get(sampleName);
                currentPositionCountTotal++;
                context2Sample2PositionCount.get("total").put(sampleName, currentPositionCountTotal);


                if (c_atc.matcher(sampleRecord.getSampleContext()).matches()) {

                    long currentCcountNonCpG = context2Sample2CCount.get("non_cpg_total").get(sampleName);
                    currentCcountNonCpG = currentCcountNonCpG + sampleRecord.getCCount();
                    context2Sample2CCount.get("non_cpg_total").put(sampleName, currentCcountNonCpG);

                    long currentTcountNonCpG = context2Sample2TCount.get("non_cpg_total").get(sampleName);
                    currentTcountNonCpG = currentTcountNonCpG + sampleRecord.getCtCoverage() - sampleRecord.getCCount();
                    context2Sample2TCount.get("non_cpg_total").put(sampleName, currentTcountNonCpG);

                    long currentPositionCountNonCpG = context2Sample2PositionCount.get("non_cpg_total").get(sampleName);
                    currentPositionCountNonCpG++;
                    context2Sample2PositionCount.get("non_cpg_total").put(sampleName, currentPositionCountNonCpG);

                }
            }

        }


        return retVal;

    }

}
