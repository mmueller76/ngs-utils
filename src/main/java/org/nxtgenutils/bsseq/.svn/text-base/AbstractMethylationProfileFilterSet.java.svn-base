package org.nxtgenutils.bsseq;

import java.util.HashSet;
import java.util.Map;
import java.util.HashMap;

import org.nxtgenutils.io.MethylationProfileRecord;

/**
 * Created by IntelliJ IDEA.
 * User: mmuelle1
 * Date: 18-Nov-2011
 * Time: 17:14:37
 * To change this template use File | Settings | File Templates.
 */
public class AbstractMethylationProfileFilterSet extends HashSet<AbstractMethylationProfileFilter> implements MethylationProfileFilterSet  {

    private Map<String, Boolean> filterResults = new HashMap<String, Boolean>();

    public AbstractMethylationProfileFilterSet() {
    }

    @Override
    public boolean passesFilterSet(MethylationProfileRecord record) {

        boolean retVal = true;

        for(MethylationProfileFilter filter : this){

            boolean passesFilter = filter.passesFilter(record);
            filterResults.put(filter.toString(), passesFilter);
            if(!passesFilter){
                retVal = false;
            }

        }

        return retVal;
    }

    @Override
    public Map<String, Boolean> getFilterResults() {
        return filterResults;
    }

}