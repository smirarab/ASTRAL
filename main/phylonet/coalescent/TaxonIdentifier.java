package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class TaxonIdentifier {
    private HashMap<String, Integer> nameToId = new HashMap<String, Integer>();
    private List<String> idToName = new ArrayList<String>();
    private int taxonCount = 0;
    private boolean locked = false;

    public String[] getAllTaxonNames(){
        return idToName.toArray(new String[]{});
    }

    public void lock() {
        this.locked = true;
    }
    public Integer taxonId(String name) {
        Integer a = nameToId.get(name);
        if (a ==  null){ 
            if (locked) {
                throw new RuntimeException("Error: "+name+" was not seen in main input trees.");
            }
            nameToId.put(name,taxonCount);
            idToName.add(name);
            a = taxonCount;
            taxonCount ++;
        }
        return a;
    }

    public String getTaxonName(Integer id) {
        return idToName.get(id);
    }
    
    public int taxonCount(){
        return taxonCount;
    }
}