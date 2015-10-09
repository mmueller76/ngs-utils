package org.nxtgenutils.io.impl;

import org.nxtgenutils.io.MultiPileup;
import org.nxtgenutils.io.Pileup;

import java.util.Set;
import java.util.LinkedHashSet;

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
 * Time: 13:23:18
 * To change this template use File | Settings | File Templates.
 */

/**
 * @author Michael Mueller
 */
public class MultiPileupImpl implements MultiPileup {

    private String sequenceName;
    private int position;
    private String referenceBase;
    private Set<Pileup> pileups = new LinkedHashSet<Pileup>();
    private String recordString;

    public MultiPileupImpl(String sequenceName, int position, String referenceBase) {
        this.sequenceName = sequenceName;
        this.position = position;
        this.referenceBase = referenceBase;
    }

    public void addPileup(Pileup pileup){
        this.pileups.add(pileup);
    }

    public int getPileupCount(){
        return this.pileups.size();
    }

    public Set<Pileup> getPileups(){
        return this.pileups;
    }

    public String getSequenceName() {
        return sequenceName;
    }

    public int getPosition() {
        return position;
    }

    public String getReferenceBase() {
        return referenceBase;
    }

    public String getRecordString() {
        return recordString;
    }

    public void setRecordString(String recordString) {
        this.recordString = recordString;
    }

    public boolean isRepeatMasked(){
        //check if bases are in lower case
        return this.getReferenceBase().matches("[a-z]");
    }

    @Override
    public String toString() {
        return "MultiPileup{" +
                "sequenceName='" + sequenceName + '\'' +
                ", position=" + position +
                ", referenceBase='" + referenceBase + '\'' +
                ", pileups=" + pileups +
                '}';
    }
    
}
