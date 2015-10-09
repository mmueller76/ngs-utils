package org.nxtgenutils.bsseq;

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
 * Date: 24-Feb-2011
 * Time: 13:40:35
 */

/**
 * @author Michael Mueller
 */
public class BSReadFilterFactory {

    /**
     * the log4j logger
     */
    private static Logger logger = Logger.getLogger(FileBSReadFilterSet.class);


    private static BSReadFilterFactory ourInstance = new BSReadFilterFactory();

    public static BSReadFilterFactory getInstance() {
        return ourInstance;
    }

    private BSReadFilterFactory() {
    }

    BSReadFilter createBSReadFilter(String filterClassName, Object value){

        BSReadFilter retVal = null;

        try {
            
            retVal = (BSReadFilter)Class.forName("org.nxtgenutils.bsseq.impl." + filterClassName).newInstance();

            if(value instanceof String){
                retVal.setFilterValueFromString((String)value);
            } else {
                retVal.setFilterValue(value);
            }

        } catch (ClassNotFoundException e) {
            logger.error("Exception while creating filter of class " + filterClassName + " and value " + value + ".", e);
        } catch (IllegalAccessException e) {
            logger.error("Exception while creating filter of class " + filterClassName + " and value " + value + ".", e);
        } catch (InstantiationException e) {
            logger.error("Exception while creating filter of class " + filterClassName + " and value " + value + ".", e);
        }

        return retVal;

    }

}
