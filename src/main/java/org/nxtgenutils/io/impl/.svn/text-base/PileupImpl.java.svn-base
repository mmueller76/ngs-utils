package org.nxtgenutils.io.impl;

import org.nxtgenutils.io.MultiPileup;
import org.nxtgenutils.io.Pileup;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.*;

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
 * Date: 26-Oct-2011
 * Time: 13:47:40
 */

/**
 * @author Michael Mueller
 */
public class PileupImpl implements Pileup {

    public MultiPileup parent;
    public int coverage;
    public String bases;
    public String basesWithoutIndels;
    public String qualities;
    private double heterozygousAlleleFrequencyCutOff = 0.4;
    private List<String> indels = null;

    public static Pattern patternReferenceMatch = Pattern.compile("[\\.|,]");
    public static Pattern patternIndel = Pattern.compile("([\\+-])([0-9]+)([ACGTNacgtn]+)");

    public PileupImpl(MultiPileup parent, int coverage, String bases, String qualities) {
        this.parent = parent;
        this.coverage = coverage;
        this.bases = bases;
        this.qualities = qualities;
    }

    public int getCoverage() {
        return coverage;
    }

    public String getBases() {
        return bases;
    }

    public int getReferenceMatchCount() {

        Matcher m = patternReferenceMatch.matcher(bases);
        int retVal = 0;
        while (m.find()) {
            retVal++;
        }
        return retVal;

    }

    public int getBaseCallCount(String base) {

        Pattern p;
        if (base.equalsIgnoreCase(parent.getReferenceBase())) {
            p = patternReferenceMatch;
        } else {
            p = Pattern.compile("[" + base.toLowerCase() + base.toUpperCase() + "]");
        }

        String basesNoIndels = getBasesWithoutIndels();
        Matcher m = p.matcher(basesNoIndels);
        int retVal = 0;
        while (m.find()) {
            retVal++;
        }
        return retVal;

    }

    public String getConsensusBaseCall() {
        return this.getConsensusBaseCall(heterozygousAlleleFrequencyCutOff);
    }

    public String getConsensusBaseCall(double alleleFrequencyCutoff) {

        String retVal = "N";
        int maxCount = 0;

        String[] bases = {"A", "T", "C", "G"};

        Map<String, Double> frequencies = new TreeMap<String, Double>();

        for (String base : bases) {
            int baseCount = this.getBaseCallCount(base);
            double frequency = (double) baseCount / (double) this.getCoverage();

            if (baseCount > maxCount) {
                maxCount = baseCount;
                retVal = base;
            }

            frequencies.put(base, frequency);

        }

        //check if there is a base call with
        //a frequency above the allele frequency
        //cutoff for heterozyous calls
        for (String base : frequencies.keySet()) {

            double frequency = frequencies.get(base);
            if (frequency >= alleleFrequencyCutoff && !base.equals(retVal)) {

                retVal = retVal + "/" + base;
                break;

            }

        }


//        if(!hetAllele.toString().equals(retVal)){
//            retVal = hetAllele.toString();
//        }

        return retVal;

    }

    public MultiPileup getPileupRecord(){
        return parent;
    }


    public String getQualities() {
        return qualities;
    }

    private List<String> parseIndels(){

        List<String> retVal = new ArrayList<String>();
        Matcher m = patternIndel.matcher(bases);

        while (m.find()) {
            String sign = m.group(1);
            int indelSize = Integer.parseInt(m.group(2));
            String indelBases = m.group(3);
            indelBases = indelBases.substring(0, indelSize);

            retVal.add(sign + indelSize + indelBases);

        }

        return retVal;
    }

    public boolean hasIndels(){
        if(indels == null){
            indels = parseIndels();
        }
        return indels.size() > 0;
    }

    public String getBasesWithoutIndels() {

        if(basesWithoutIndels != null){
            return basesWithoutIndels;
        }

        basesWithoutIndels = bases;
        if(indels == null){
            indels = parseIndels();
        }

        for (String indel : indels) {
            basesWithoutIndels = basesWithoutIndels.replace(indel, "");
        }
        
        return basesWithoutIndels;

    }

    @Override
    public String toString() {
        return "Pileup{" +
                "coverage=" + coverage +
                ", bases='" + bases + '\'' +
                ", qualities='" + qualities + '\'' +
                '}';
    }

}
