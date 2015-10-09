package org.nxtgenutils.bsseq;

import org.apache.log4j.Logger;

import java.util.Properties;
import java.util.Iterator;
import java.io.FileInputStream;
import java.io.File;
import java.io.IOException;

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
 * Date: 03-Feb-2011
 * Time: 13:06:56
 */

/**
 * @author Michael Mueller
 */
public class FileBSReadFilterSet extends AbstractBSReadFilterSet {

    /**
      * the log4j logger
      */
    private static Logger logger = Logger.getLogger(FileBSReadFilterSet.class);


    private Properties filterProperties;


    public FileBSReadFilterSet(String filePath){

       this(new File(filePath));

    }

    public FileBSReadFilterSet(File file){

        filterProperties = new Properties();
        try {
            filterProperties.load(new FileInputStream(file));
        } catch (IOException e) {
            logger.error("Unable to read filter properties from file " + file.getAbsolutePath() + "." + e.getMessage());            
        }
        readFilterProperties();

    }

    private void readFilterProperties() {

        Iterator propertyIterator = filterProperties.keySet().iterator();

        while (propertyIterator.hasNext()) {

            String filterClassName =  (String)propertyIterator.next();
            String value = filterProperties.getProperty(filterClassName);

            this.add(BSReadFilterFactory.getInstance().createBSReadFilter(filterClassName, value));

        }

    }

    public void writeFilterProperties() {

        for(BSReadFilter filter : this){

            filterProperties.setProperty(filter.getClass().toString().replace("class org.nxtgenutils.bsseq.impl.",""), filter.getFilterValue().toString());

        }

    }

    public static void main(String[] args) {

        String filterPropertiesFile = "/home/mmuelle1/dev/java/bs-utils/testdata/bs_read_filter.properties";
        BSReadFilterSet rfs = new FileBSReadFilterSet(filterPropertiesFile);

        for(BSReadFilter filter : rfs){

            System.out.println(filter.toString());
            

        }

    }


}
