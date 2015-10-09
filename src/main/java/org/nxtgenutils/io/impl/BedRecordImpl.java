package org.nxtgenutils.io.impl;

import org.nxtgenutils.Strand;
import org.nxtgenutils.io.BedRecord;

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
 * Date: 25-Oct-2011
 * Time: 17:07:44
 */

/**
 * Implementation of the BedRecord interface.
 *
 * @author Michael Mueller
 */
public class BedRecordImpl implements BedRecord {

    private String chromosome;
    private int start;
    private int end;
    private String name;
    private int scores;
    private Strand strand;

    /**
     * Constructs a new BedRecord with the mandatory fields.
     *
     * @param chromosome the chromosome name
     * @param start the feature start (0 based)
     * @param end the feature end (1 based)
     * @param name the feature name
     */
    public BedRecordImpl(String chromosome, int start, int end, String name) {
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.name = name;
    }

    /**
     * Returns chromosome.
     *
     * @return the chromosome name
     */
    @Override
    public String getChromosome() {
        return chromosome;
    }

    /**
     * Sets chromosome.
     *
     * @param chromosome the chromosome name
     */
    @Override
    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    /**
     * Returns feature start (0 based).
     *
     * @return the feature start
     */
    @Override
    public int getStart() {
        return start;
    }

    /**
     * Sets feature start (0 based).
     *
     * @param start the feature start
     */
    @Override
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * Returns feature end (1 based).
     *
     * @return the feature end
     */
    @Override
    public int getEnd() {
        return end;
    }

    /**
     * Sets feature end (1 based).
     *
     * @param end the feature end
     */
    @Override
    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * Returns feature name.
     *
     * @return the feature name
     */
    @Override
    public String getName() {
        return name;
    }

    /**
     * Sets feature name.
     *
     * @param name the feature name
     */
    @Override
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Returns feature score.
     *
     * @return the score
     */
    @Override
    public int getScores() {
        return scores;
    }

    /**
     * Sets feature score.
     *
     * @param scores the score
     */
    @Override
    public void setScores(int scores) {
        this.scores = scores;
    }

    /**
     * Returns feature strand.
     *
     * @return the strand
     */
    @Override
    public Strand getStrand() {
        return strand;
    }

    /**
     * Sets feature strand.
     *
     * @param strand the strand
     */
    @Override
    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BedRecordImpl)) return false;

        BedRecordImpl bedRecord = (BedRecordImpl) o;

        if (end != bedRecord.end) return false;
        if (scores != bedRecord.scores) return false;
        if (start != bedRecord.start) return false;
        if (chromosome != null ? !chromosome.equals(bedRecord.chromosome) : bedRecord.chromosome != null) return false;
        if (name != null ? !name.equals(bedRecord.name) : bedRecord.name != null) return false;
        if (strand != bedRecord.strand) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = chromosome != null ? chromosome.hashCode() : 0;
        result = 31 * result + start;
        result = 31 * result + end;
        result = 31 * result + (name != null ? name.hashCode() : 0);
        result = 31 * result + scores;
        result = 31 * result + (strand != null ? strand.hashCode() : 0);
        return result;
    }
}
