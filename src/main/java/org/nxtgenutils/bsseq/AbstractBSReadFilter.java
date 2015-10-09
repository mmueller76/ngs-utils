package org.nxtgenutils.bsseq;

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
 * Date: 24-Feb-2011
 * Time: 17:00:00
 */

/**
 * @author Michael Mueller
 * @param <V>
 */
public abstract class AbstractBSReadFilter<V> implements BSReadFilter<V> {

    V filterValue;

    @Override
    public V getFilterValue() {
        return filterValue;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void setFilterValue(V value) {
        this.filterValue=value;
    }
    
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof AbstractBSReadFilter)) return false;

        AbstractBSReadFilter that = (AbstractBSReadFilter) o;

        if (!filterValue.equals(that.filterValue)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return filterValue.hashCode();
    }

}
