package org.nxtgenutils.bsseq;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecordUtil;
import net.sf.samtools.SAMValidationError;

import java.util.HashMap;
import java.util.TreeMap;
import java.util.Map;
import java.io.File;
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
 * Date: 18-Feb-2011
 * Time: 10:06:15
 */

/**
 * @author Michael Mueller
 */
public abstract class AbstractBSMappingProcessor implements BSMappingProcessor {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(AbstractBSMappingProcessor.class);

    /**
     * how many reads had CIGAR string but were flagged as unmapped
     */
    protected int readMappedFlagFixedCount = 0;
    protected int readCigarStringFixedCount = 0;
    protected int mappingForwardStrandOnly = 0;
    protected int mappingReverseStrandOnly = 0;
    protected int mapingBothStrands = 0;
    protected int unmapped = 0;
    protected int properlyPairedForwardStrandOnly = 0;
    protected int properlyPairedReverseStrandOnly = 0;
    protected int properlyPairedBothStrands = 0;
    protected int notProperlyPaired = 0;
    protected int properlyPairedBothStrandsAndRepeat = 0;
    protected Map<Integer, Integer> nonBSMismatchDistributionForwardStrand = new TreeMap<Integer, Integer>();
    protected Map<Integer, Integer> nonBSMismatchDistributionReverseStrand = new TreeMap<Integer, Integer>();
    protected Map<Integer, Integer> bsInvalidMismatchDistributionForwardStrand = new TreeMap<Integer, Integer>();
    protected Map<Integer, Integer> bsInvalidMismatchDistributionReverseStrand = new TreeMap<Integer, Integer>();
    protected Map<String, Map<Integer, Integer>> bsConversionRate = new TreeMap<String, Map<Integer, Integer>>();
    protected static DecimalFormat oneDecimalPlace = new DecimalFormat("#.#");
    protected static DecimalFormat twoDecimalPlaces = new DecimalFormat("#.##");
    protected int readCount = 0;
    protected int readPairCount = 0;
    protected int invalidMappingOrientation = 0;
    protected Map<String, ReferenceSequence> referenceSequence;

    protected void loadReferenceSequences(ReferenceSequenceFile referenceSequenceFile) {

        logger.info("loading reference sequence");
        referenceSequence = new HashMap<String, ReferenceSequence>();
        ReferenceSequence rs;
        while ((rs = referenceSequenceFile.nextSequence()) != null) {


            String name = rs.getName();

            logger.info(name);
            referenceSequence.put(name, rs);


        }

//        rs = referenceSequenceFile.getSequence("chr1");
//        String name = rs.getName();
//        logger.info(name);
//        referenceSequence.put(name, rs);

    }

    protected abstract void printStatistics(File outputPath);

    protected abstract void appendBSMappingOrientationTag(SAMRecord read);

    /**
     * Returns the number of converted and unconverted Cs in a
     * bisulphite converted read with respect to the mapping
     * position.
     *
     * @param read            the read
     * @param referenceBases  the reference sequence the read is mapped to
     * @param referenceOffset the offset if reference bases are a sub-sequence of the reference
     * @param cpgContext      if true conversion in CpG context is assessed; outside CpG context otherwise
     * @return integer array of length 2, the first value holds the number of unconverted Cs and
     *         the second value the number of converted Cs
     */
    protected int[] getBSConversion(final SAMRecord read, final byte[] referenceBases, final int referenceOffset, final boolean cpgContext) {


        int[] retVal = new int[2];

        final byte[] readBases = read.getReadBases();

        boolean reverseStrandTemplate = false;
        if (read.getAttribute("BT").toString().equals("R")) {
            reverseStrandTemplate = true;
        }

        for (final AlignmentBlock block : read.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - 1 - referenceOffset;
            final int length = block.getLength();


            for (int i = 0; i < length; i++) {

                try {

                    byte readBase = readBases[readBlockStart + i];
                    byte referenceBase = referenceBases[referenceBlockStart + i];
                    byte previousReferenceBase = (byte) 'N';
                    if (!((referenceBlockStart + i - 1) < 0)) {
                        previousReferenceBase = referenceBases[referenceBlockStart + i - 1];
                    }

                    byte nextReferenceBase = (byte) 'N';
                    if (!((referenceBlockStart + i - 1) > referenceBases.length - 1)) {
                        nextReferenceBase = referenceBases[referenceBlockStart + i + 1];
                    }

                    if (!reverseStrandTemplate) {

                        //if unconverted C in read from forward strand template in CpG context
                        if (basesEqual(readBase, C) &&
                                basesEqual(referenceBase, C) &&
                                ((cpgContext && basesEqual(nextReferenceBase, G))
                                        || (!cpgContext && !basesEqual(nextReferenceBase, G)))) {

                            retVal[0]++;

                            //if converted C in read from forward strand template in CpG context
                        } else if (basesEqual(readBase, T) &&
                                basesEqual(referenceBase, C) &&
                                ((cpgContext && basesEqual(nextReferenceBase, G))
                                        || (!cpgContext && !basesEqual(nextReferenceBase, G)))) {

                            retVal[1]++;

                        }

                    } else {

                        //if unconverted C in read from reverse strand template in CpG context
                        if (basesEqual(readBase, G) &&
                                basesEqual(referenceBase, G) &&
                                ((cpgContext && basesEqual(previousReferenceBase, C))
                                     || (!cpgContext && !basesEqual(previousReferenceBase, C)))) {

                            retVal[0]++;

                            //if converted C in read from reverse strand template in CpG context
                        } else if (basesEqual(readBase, A) &&
                                basesEqual(referenceBase, G) &&
                                ((cpgContext && basesEqual(previousReferenceBase, C))
                                    || (!cpgContext && !basesEqual(previousReferenceBase, C)))) {

                            retVal[1]++;

                        }


                    }

                } catch (ArrayIndexOutOfBoundsException e) {

                    logger.error("Warning: exception while determining BS conversion rate for for read " + read.getReadName() + ": " + e.getMessage());

                }

            }


        }

        return retVal;

    }

    /**
     *
     * @param read
     */
    protected abstract void appendBSConversionTags(SAMRecord read);

    protected void appendBSTemplateStrandTag(SAMRecord record, BSTemplateStrand bsMappingStrand) {

        if (bsMappingStrand == BSTemplateStrand.FORWARD) {
            record.setAttribute("BT", 'F');
        } else if (bsMappingStrand == BSTemplateStrand.REVERSE) {
            record.setAttribute("BT", 'R');
        }

    }

    protected BSMappingStrandSpecificity appendMappingStrandSpecificityTag(SAMRecord recordfowardStrandMapping, SAMRecord recordReverseStrandMapping) {
        //check if read maps only to one strand
        BSMappingStrandSpecificity retVal = BSMappingStrandSpecificity.UNMAPPED;

        if (!recordfowardStrandMapping.getReadUnmappedFlag() && recordReverseStrandMapping.getReadUnmappedFlag()) {

            retVal = BSMappingStrandSpecificity.FORWARD;
            recordfowardStrandMapping.setAttribute("BS", 'F');
            recordReverseStrandMapping.setAttribute("BS", 'F');
            mappingForwardStrandOnly++;

        } else if (recordfowardStrandMapping.getReadUnmappedFlag() && !recordReverseStrandMapping.getReadUnmappedFlag()) {

            retVal = BSMappingStrandSpecificity.REVERSE;
            recordfowardStrandMapping.setAttribute("BS", 'R');
            recordReverseStrandMapping.setAttribute("BS", 'R');
            mappingReverseStrandOnly++;

        } else if (!recordfowardStrandMapping.getReadUnmappedFlag() && !recordReverseStrandMapping.getReadUnmappedFlag()) {

            retVal = BSMappingStrandSpecificity.BOTH;
            recordfowardStrandMapping.setAttribute("BS", 'B');
            recordReverseStrandMapping.setAttribute("BS", 'B');
            mapingBothStrands++;

        } else if (recordfowardStrandMapping.getReadUnmappedFlag() && recordReverseStrandMapping.getReadUnmappedFlag()) {

            retVal = BSMappingStrandSpecificity.UNMAPPED;
            recordfowardStrandMapping.setAttribute("BS", 'N');
            recordReverseStrandMapping.setAttribute("BS", 'N');
            unmapped++;

        }

        return retVal;

    }

    /**
     * @param read
     */
    protected void appendMismatchTags(SAMRecord read) {

        if (!read.getReadUnmappedFlag()) {

            String referenceSequenceName = read.getReferenceName().replace(".ct", "").replace(".ga", "");
            ReferenceSequence rs = referenceSequence.get(referenceSequenceName);

            if(rs == null){
                logger.error("Error while processing read " + read.getReadName() + ". Reference sequence with name " + referenceSequenceName + " not in reference sequence file.");
            }

            int totalMismatches = countMismatches(read, rs.getBases(), 0, false);
            int bsInvalidMismatches = countBsInvalidMismatches(read, rs.getBases(), 0);
            int nonBsMismatches = countMismatches(read, rs.getBases(), 0, true) - bsInvalidMismatches;
            int bsValidMismatches = totalMismatches - nonBsMismatches - bsInvalidMismatches;

            read.setAttribute("BN", nonBsMismatches);
            read.setAttribute("BI", bsInvalidMismatches);
            read.setAttribute("BV", bsValidMismatches);


            BSTemplateStrand bsTemplateStrand = BSTemplateStrand.FORWARD;
            if (read.getAttribute("BT").toString().equals("R")) {
                bsTemplateStrand = BSTemplateStrand.REVERSE;
            }

            if (bsTemplateStrand == BSTemplateStrand.FORWARD) {

                if (!nonBSMismatchDistributionForwardStrand.containsKey(nonBsMismatches)) {
                    nonBSMismatchDistributionForwardStrand.put(nonBsMismatches, 0);
                }
                int freq = nonBSMismatchDistributionForwardStrand.get(nonBsMismatches);
                freq++;
                nonBSMismatchDistributionForwardStrand.put(nonBsMismatches, freq);

                if (!bsInvalidMismatchDistributionForwardStrand.containsKey(bsInvalidMismatches)) {
                    bsInvalidMismatchDistributionForwardStrand.put(bsInvalidMismatches, 0);
                }
                freq = bsInvalidMismatchDistributionForwardStrand.get(bsInvalidMismatches);
                freq++;
                bsInvalidMismatchDistributionForwardStrand.put(bsInvalidMismatches, freq);

            } else {

                if (!nonBSMismatchDistributionReverseStrand.containsKey(nonBsMismatches)) {
                    nonBSMismatchDistributionReverseStrand.put(nonBsMismatches, 0);
                }

                int freq = nonBSMismatchDistributionReverseStrand.get(nonBsMismatches);
                freq++;
                nonBSMismatchDistributionReverseStrand.put(nonBsMismatches, freq);

                if (!bsInvalidMismatchDistributionReverseStrand.containsKey(bsInvalidMismatches)) {
                    bsInvalidMismatchDistributionReverseStrand.put(bsInvalidMismatches, 0);
                }
                freq = bsInvalidMismatchDistributionReverseStrand.get(bsInvalidMismatches);
                freq++;
                bsInvalidMismatchDistributionReverseStrand.put(bsInvalidMismatches, freq);

            }

        }
    }

    /**
     * Calculates the number of mismatches between the read and the reference sequence provided.
     */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases) {
        return countMismatches(read, referenceBases, 0, false);
    }

    /**
     * Calculates the number of mismatches between the read and the reference sequence provided.
     */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases, final int referenceOffset) {
        return countMismatches(read, referenceBases, referenceOffset, false);
    }

    /**
     * Calculates the number of mismatches between the read and the reference sequence provided.
     *
     * @param referenceBases    Array of ASCII bytes that covers at least the the portion of the reference sequence
     *                          to which read is aligned from getReferenceStart to getReferenceEnd.
     * @param referenceOffset   0-based offset of the first element of referenceBases relative to the start
     *                          of that reference sequence.
     * @param bisulfiteSequence If this is true, it is assumed that the reads were bisulfite treated
     *                          and C->T on the positive strand and G->A on the negative strand will not be counted
     *                          as mismatches.
     */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases, final int referenceOffset,
                                      final boolean bisulfiteSequence) {
        int mismatches = 0;

        final byte[] readBases = read.getReadBases();

        boolean reverseStrandTemplate = false;
        if (read.getAttribute("BT").toString().equals("R")) {
            reverseStrandTemplate = true;
        }

        for (final AlignmentBlock block : read.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - 1 - referenceOffset;
            final int length = block.getLength();

            for (int i = 0; i < length; ++i) {

                try {


                    if (!bisulfiteSequence) {

                        if (!basesEqual(readBases[readBlockStart + i], referenceBases[referenceBlockStart + i])) {
                            ++mismatches;
                        }

                    } else {
                        if (!bisulfiteBasesEqual(reverseStrandTemplate, readBases[readBlockStart + i],
                                referenceBases[referenceBlockStart + i])) {
                            ++mismatches;
                        }
                    }

                } catch (ArrayIndexOutOfBoundsException e) {

                    logger.error("Warning: exception while counting mismatches for read " + read.getReadName() + ": " + e.getMessage());

                }
            }
        }
        return mismatches;
    }

    public static int countBsInvalidMismatches(final SAMRecord read, final byte[] referenceBases, final int referenceOffset) {
        int mismatches = 0;

        final byte[] readBases = read.getReadBases();

        boolean reverseStrandTemplate = false;
        if (read.getAttribute("BT").toString().equals("R")) {
            reverseStrandTemplate = true;
        }

        for (final AlignmentBlock block : read.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - 1 - referenceOffset;
            final int length = block.getLength();


            for (int i = 0; i < length; ++i) {

                try {
                    if (bisulfiteInvalid(reverseStrandTemplate, readBases[readBlockStart + i],
                            referenceBases[referenceBlockStart + i])) {
                        ++mismatches;
                    }
                } catch (ArrayIndexOutOfBoundsException e) {

                    logger.error("Warning: exception while counting BS invalid mismatches for read " + read.getReadName() + ": " + e.getMessage());

                }

            }


        }
        return mismatches;
    }

    /**
     * Attempts to efficiently compare two bases stored as bytes for equality.
     */
    public static boolean basesEqual(byte lhs, byte rhs) {
        if (lhs == rhs) return true;
        else {
            if (lhs > 90) lhs -= 32;
            if (rhs > 90) rhs -= 32;
        }

        return lhs == rhs;
    }

    /**
     * Returns true if the bases are equal OR if the mismatch cannot be accounted for by
     * bisfulite treatment.  C->T on the positive strand and G->A on the negative strand
     * do not count as mismatches
     */
    public static boolean bisulfiteBasesEqual(boolean reverseTemplateStrand, byte read, byte reference) {


        if (basesEqual(read, reference)) {
            return true;
        }

        if (reverseTemplateStrand) {
            if (basesEqual(reference, (byte) 'G') && basesEqual(read, (byte) 'A')) {
                return true;
            }
        } else {
            if (basesEqual(reference, (byte) 'C') && basesEqual(read, (byte) 'T')) {
                return true;
            }
        }

        return false;

    }

    public static boolean bisulfiteInvalid(boolean reverseTemplateStrand, byte read, byte reference) {


        if (basesEqual(read, reference)) {
            return false;
        }

        if (reverseTemplateStrand) {
            if (basesEqual(reference, (byte) 'A') && basesEqual(read, (byte) 'G')) {
                return true;
            }
        } else {
            if (basesEqual(reference, (byte) 'T') && basesEqual(read, (byte) 'C')) {
                return true;
            }
        }

        //System.out.println("bisulfiteBasesEqual: read = " + byteToChar(read) + " reference = " + byteToChar(reference) + " -> " + retVal);

        return false;
    }

    public static char byteToChar(final byte b) {

        switch (b) {
            case t:
                return 't';
            case g:
                return 'g';
            case c:
                return 'c';
            case a:
                return 'a';
            case T:
                return 'T';
            case G:
                return 'G';
            case C:
                return 'C';
            case A:
                return 'A';
            default:
                return (char) b;
        }
    }

    public static String byteArrayToString(byte[] sequence) {

        StringBuffer retVal = new StringBuffer();

        for (byte b : sequence) {

            retVal.append(byteToChar(b));

        }

        return retVal.toString();

    }

    public static boolean readIsEqual(SAMRecord record1, SAMRecord record2) {

        if (record1 == null || record2 == null) {
            return false;
        }

        boolean nameIsEqual = record1.getReadName().equals(record2.getReadName());
        String read1String = "";
        String read2String = "";

        if (record1.getReadNegativeStrandFlag()) {
            read1String = getReverseComplementReadString(record1);
        } else {
            read1String = record1.getReadString();
        }

        if (record2.getReadNegativeStrandFlag()) {
            read2String = getReverseComplementReadString(record2);
        } else {
            read2String = record2.getReadString();
        }

        boolean sequenceIsEqual = read1String.equals(read2String);

        return nameIsEqual && sequenceIsEqual;

    }

    public static String getReverseComplementReadString(SAMRecord record) {

        String retVal;

        //reverse sequence...
        SAMRecordUtil.reverseComplement(record);
        //get sequence string
        retVal = record.getReadString();
        //...and back again
        SAMRecordUtil.reverseComplement(record);

        return retVal;
    }

    protected void fixRecord(SAMRecord record) {


        if (record.isValid() == null) {
            return;
        }

        for (SAMValidationError error : record.isValid()) {

            //fix read mapped flag
            //System.out.println(error.getMessage());
            if (error.getMessage().equalsIgnoreCase("CIGAR should have zero elements for unmapped read.")) {
                //System.out.println(record.getMateReferenceName() + " " + record.getAlignmentStart());
                if (record.getMateReferenceName() != null &&
                        record.getAlignmentStart() != 0 &&
                        (!record.getMateReferenceName().equals(record.getReferenceName()) &&
                                record.getMateAlignmentStart() != record.getAlignmentStart())) {
                    record.setReadUnmappedFlag(false);
                    readMappedFlagFixedCount++;
                } else {
                    record.setCigarString("*");
                    readCigarStringFixedCount++;
                }

            }

        }

    }
}
