package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.sti.STITreeCluster;

public class TaxonIdentifier {
    protected HashMap<String, Integer> nameToId = new HashMap<String, Integer>();
    protected List<String> idToName = new ArrayList<String>();
    protected int taxonCount = 0;
    protected boolean locked = false;
    
    public String[] getAllTaxonNames(){
        return idToName.toArray(new String[]{});
    }

    //public STITreeCluster newCluster() {
    //	return new STITreeCluster(this);
    //}
    
    public void lock() {
        this.locked = true;
    }
    /**
     * Returns the ID corresponding to the taxon name , if the taxon name is new it adds it unless the taxonidentifier is locked
     * @param name
     * @return
     */
    public Integer taxonId(String name) {
    	if ("".equals(name)) {
    		throw new RuntimeException("Empty name observed; likely, an input tree has an error");
    	}
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
    
	public  STITreeCluster getClusterForNodeName(String nodeName) {
		STITreeCluster cluster = Factory.instance.newCluster(this);
		Integer taxonID = this.taxonId(nodeName);
		cluster.addLeaf(taxonID);
		return cluster;
	}
}