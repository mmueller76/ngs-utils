package org.nxtgenutils.impl;

import net.sf.picard.sam.FixMateInformation;
import net.sf.picard.util.CigarUtil;
import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.nxtgenutils.PrimerLocation;
import org.nxtgenutils.Strand;
import org.nxtgenutils.io.BedRecord;
import org.nxtgenutils.io.impl.SimpleBedParser;

import java.io.File;
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
 * Created with IntelliJ IDEA.
 * User: mmuelle1
 * Date: 8/8/13
 * Time: 11:09 AM
 * To change this template use File | Settings | File Templates.
 */
public class PrimerClipper {


    private File mappingInput;
    private File mappingOutput;
    private File primerCoordinateBed;
    private SAMFileHeader mappingOutputHeader;
    private int primerCoordinatesOffSet = 2;

    private static Logger logger = Logger.getLogger(PrimerClipper.class);

    public PrimerClipper(File mappingInput, File mappingOutput, File primerCoordinateBed){
        this(mappingInput, mappingOutput, primerCoordinateBed, 2);
    }

    public PrimerClipper(File mappingInput, File mappingOutput, File primerCoordinateBed, int primerCoordinatesOffSet) throws IllegalStateException {

        this.mappingInput = mappingInput;
        this.mappingOutput = mappingOutput;
        this.primerCoordinateBed = primerCoordinateBed;
        this.primerCoordinatesOffSet = primerCoordinatesOffSet;

        //create SAM reader & set validation stringency
        SAMFileReader samFileReader = new SAMFileReader(mappingInput);
        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        mappingOutputHeader = samFileReader.getFileHeader();
        mappingOutputHeader.addComment("primer/probe sequences soft clipped");
        mappingOutputHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        samFileReader.close();

        SimpleBedParser bedParser = new SimpleBedParser(this.primerCoordinateBed);

        Iterator<BedRecord> bedIterator = bedParser.iterator();

//        Map<String, Map<BedRecord, PrimerLocation>> primerCoordByChr = new LinkedHashMap<String, Map<BedRecord, PrimerLocation>>();

        Map<String, Set<BedRecord>> primerCoordByChr = new LinkedHashMap<String, Set<BedRecord>>();

        //determine primer amplicon location
        //start or end of amplicon relative to the
        // forward strand reference sequence
//        logger.info("Determining amplicon primer location relative to reference sequence forward strand...");
//
        int count = 0;
        while (bedIterator.hasNext()) {

            count++;
            BedRecord primerCoord = bedIterator.next();
            String primerChr = primerCoord.getChromosome();
            if (!primerCoordByChr.containsKey(primerChr)) {
                primerCoordByChr.put(primerChr, new LinkedHashSet<BedRecord>());
            }

            primerCoordByChr.get(primerChr).add(primerCoord);

        }




        logger.info(count + " primer coordinates in input file.");


        logger.info("Soft clipping primer sequence...");

        //create SAM reader & set validation stringency
        samFileReader = new SAMFileReader(mappingInput);
        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        Iterator<SAMRecord> readIterator = samFileReader.iterator();

        SAMFileWriter bamFileWriter = new SAMFileWriterFactory().makeBAMWriter(mappingOutputHeader,
                true, this.mappingOutput);

        int readProcessedCount = 0;
        int readClippedCount = 0;

        while (readIterator.hasNext()) {

            readProcessedCount++;

            SAMRecord read = readIterator.next();

            if (softClipPrimerSequence(read, primerCoordByChr)) {
                readClippedCount++;
            }

            bamFileWriter.addAlignment(read);

            if (readProcessedCount % 1000000 == 0) {
                logger.info(readProcessedCount + " records processed, " + readClippedCount + " records clipped.");
            }

        }

        logger.info(readProcessedCount + " records processed, " + readClippedCount + " records clipped.");
        logger.info("Fixing mate information and re-sorting clipped reads...");

        bamFileWriter.close();

        String[] fixMateArgs = new String[3];
        fixMateArgs[0] = "INPUT=" + mappingOutput.getAbsolutePath();
        fixMateArgs[1] = "SORT_ORDER=coordinate";
        fixMateArgs[2] = "VALIDATION_STRINGENCY=SILENT";

        FixMateInformation.main(fixMateArgs);

        logger.info("done");

    }

    private PrimerLocation determinePrimerLocation(BedRecord primerCoordinates) {

        PrimerLocation retVal = PrimerLocation.UNKNOWN;

        String primerChr = primerCoordinates.getChromosome();
        int primerStart = primerCoordinates.getStart() + 1;
        int primerEnd = primerCoordinates.getEnd();

        //create SAM reader & set validation stringency
        SAMFileReader samFileReader = new SAMFileReader(mappingInput);
        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        Iterator<SAMRecord> readIterator = samFileReader.queryOverlapping(primerChr, primerStart, primerEnd);

        int readOverlapAmplStartCount = 0;
        int readOverlapAmplEndCount = 0;

        while (readIterator.hasNext()) {


            SAMRecord read = readIterator.next();
            int alignStart = read.getUnclippedStart();
            int alignEnd = read.getUnclippedEnd();


            if ((alignStart >= primerStart - primerCoordinatesOffSet && alignStart <= primerStart + primerCoordinatesOffSet) &&
                    primerEnd < alignEnd) {

                readOverlapAmplStartCount++;

            }

            if ((alignEnd >= primerEnd - primerCoordinatesOffSet && alignEnd <= primerEnd + primerCoordinatesOffSet) &&
                    primerStart < alignEnd) {

                readOverlapAmplEndCount++;

            }


        }

        samFileReader.close();

        //if target region start primer
        if (readOverlapAmplStartCount > readOverlapAmplEndCount) {
            retVal = PrimerLocation.AMPLICON5PRIME;
        }

        //if target region start primer
        if (readOverlapAmplStartCount < readOverlapAmplEndCount) {
            retVal = PrimerLocation.AMPLICON3PRIME;
        }

        return retVal;

    }

    private boolean softClipPrimerSequence(SAMRecord read, Map<String, Set<BedRecord>> primerCoordinatesByChromosome) {

        String readChr = read.getReferenceName();

        int alignStart = read.getUnclippedStart();
        int alignEnd = read.getUnclippedEnd();

        if (primerCoordinatesByChromosome.containsKey(readChr)) {

            Set<BedRecord> primerCoordinates = primerCoordinatesByChromosome.get(readChr);

            for (BedRecord primer : primerCoordinates) {

                String primerName = primer.getName();
                int primerStart = primer.getStart();
                int primerEnd = primer.getEnd();
                Strand primerStrand = primer.getStrand();

                int clipFrom = 0;
//                if (primerAmpliconLocation.equals(PrimerLocation.AMPLICON5PRIME) &&
//                        (alignStart >= primerStart - primerCoordinatesOffSet && alignStart <= primerStart + primerCoordinatesOffSet) &&
//                        primerEnd < alignEnd) {

                if (primerStrand.equals(Strand.FORWARD) &&
                        !read.getReadNegativeStrandFlag() &&
                            (alignStart >= primerStart - primerCoordinatesOffSet && alignStart <= primerStart + primerCoordinatesOffSet) &&
                            primerEnd < alignEnd) {

                    clipFrom = primerEnd - alignStart + 1;

//                    System.out.println("alignStart="+alignStart);
//                    System.out.println("alignEnd="+alignEnd);
//                    System.out.println("primerName="+primerName);
//                    System.out.println("primerStart="+primerStart);
//                    System.out.println("primerEnd="+primerEnd);
//                    System.out.println("clipFrom="+clipFrom);
//                    System.out.println("primerStrand="+primerStrand);

                }

                if (primerStrand.equals(Strand.REVERSE) &&
                        read.getReadNegativeStrandFlag() &&
                        (alignEnd >= primerEnd - primerCoordinatesOffSet && alignEnd <= primerEnd + primerCoordinatesOffSet) &&
                        primerStart < alignEnd) {

//                    if (primerAmpliconLocation.equals(PrimerLocation.AMPLICON3PRIME) &&
//                        (alignEnd >= primerEnd - primerCoordinatesOffSet && alignEnd <= primerEnd + primerCoordinatesOffSet) &&
//                        primerStart < alignEnd) {

                    clipFrom = alignEnd - primerStart + 1;

//                    System.out.println("alignStart="+alignStart);
//                    System.out.println("alignEnd="+alignEnd);
//                    System.out.println("primerName="+primerName);
//                    System.out.println("primerStart="+primerStart);
//                    System.out.println("primerEnd="+primerEnd);
//                    System.out.println("clipFrom="+clipFrom);
//                    System.out.println("primerStrand="+primerStrand);

                }


                //if overlap found,
                //clip read and return
                //true
                if (clipFrom > 0) {

                    int invertedClipFrom = read.getReadLength() - clipFrom + 1;

                    //if the read maps on the forward strand, we have to invert
                    //the cigar string to clip from the 5' end
                    if (!read.getReadNegativeStrandFlag()) {

                        List<CigarElement> invertedCigar = invertCigarElements(read.getCigar().getCigarElements());

                        List<CigarElement> clippedInvertedCigar = CigarUtil.softClipEndOfRead(invertedClipFrom, invertedCigar);

                        List<CigarElement> clippedCigar = invertCigarElements(clippedInvertedCigar);

                        //String cigarString = CigarUtil.cigarStringFromArray(CigarUtil.cigarArrayFromElements(clippedCigar));
                        read.setCigar(new Cigar(clippedCigar));

                        //adjust alignment start
                        int clippedAlignmentStart = alignStart + clipFrom;
                        read.setAlignmentStart(clippedAlignmentStart);


                    } else {

                        List<CigarElement> clippedCigar = CigarUtil.softClipEndOfRead(invertedClipFrom, read.getCigar().getCigarElements());

                        //String cigarString = CigarUtil.cigarStringFromArray(CigarUtil.cigarArrayFromElements(clippedCigar));
                        read.setCigar(new Cigar(clippedCigar));


                    }

                    read.setAttribute("PC", clipFrom);

                    return true;

                }

            }

        }

        return false;

    }

    public List<CigarElement> invertCigarElements(List<CigarElement> cigar) {

        List<CigarElement> retVal = new ArrayList<CigarElement>();
        for (int i = cigar.size() - 1; i >= 0; i--) {
            retVal.add(cigar.get(i));
        }

        return retVal;

    }

    public static void main(String[] args) {

        //new ReadPairClipper(new File("/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.rn4.ct.fa.sorted.nophix.processed.head.bam"));
        //new ReadPairClipper(new File("/home/mmuelle1/dev/java/bs-utils/testdata/100805_SN172_0185_B20473ABXX_L3_1_sequence.fq.vs.rn4.ga.fa.sampe.uncon.read.readname.sorted.bam"));
//        new PrimerClipper(new File("/home/mmuelle1//dev/java/nxtgen-utils/testdata/5154.testamplicon.bam"),
//                new File("/home/mmuelle1//dev/java/nxtgen-utils/testdata/5154.testamplicon.clipped.bam"),
//                new File("/home/mmuelle1/dev/java/nxtgen-utils/testdata/testamplicon_primer.bed"));

//        new PrimerClipper(new File("/home/mmuelle1/Documents/aitman_ph/5154.bam"),
//                new File("/home/mmuelle1/Documents/aitman_ph/5154.clipped.bam"),
//                new File("/home/mmuelle1/Documents/aitman_ph/truseq_probe_coords.bed"));

//        new PrimerClipper(new File("/home/mmuelle1/dev/java/nxtgen-utils/testdata/1027.realigned.recalibrated.2013-08-16.test.bam"),
//                new File("/home/mmuelle1/dev/java/nxtgen-utils/testdata/1027.realigned.recalibrated.2013-08-16.test.clipped.bam"),
//                new File("/home/mmuelle1/dev/java/nxtgen-utils/testdata/aitman_fh_primer_coords.AA_11v_1085_final.bed"));

        new PrimerClipper(new File("/home/mmuelle1/Documents/aitman_fh/1027.chunk_1.realigned.unclipped.2013-08-22.igv.bam"),
                new File("/home/mmuelle1/Documents/aitman_fh/1027.chunk_1.realigned.unclipped.2013-08-22.igv.clipped.bam"),
                new File("/home/mmuelle1/dev/java/nxtgen-utils/testdata/aitman_fh_primer_coords.AA_11v_1085_final.bed"));


    }

}
