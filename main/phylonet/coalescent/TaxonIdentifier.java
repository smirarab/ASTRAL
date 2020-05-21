package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class TaxonIdentifier {
    private HashMap<String, Integer> nameToId = new HashMap<String, Integer>();
    private List<String> idToName = new ArrayList<String>();
    private int taxonCount = 0;
    private boolean locked = false;
    
	/**
	 * Hash values for each taxon
	 */
	public long[] hash1, hash2;
	


    public String[] getAllTaxonNames(){
        return idToName.toArray(new String[]{});
    }

    public STITreeCluster newCluster() {
    	return new STITreeCluster(this);
    }
    
    public STITreeCluster newCluster(BitSet bitset) {
    	return new STITreeCluster(this,bitset);
    }
    
    public void lock() {
        this.locked = true;
	hash1 = new long[taxonCount]; 
	hash2 = new long[taxonCount];
	Random rnd = GlobalMaps.random;
	for (int i = 0; i < taxonCount; i++) {
		hash1[i] = rnd.nextLong() ^ rnd.nextLong() ^ rnd.nextLong(); 
		hash2[i] = ( rnd.nextLong() ^ rnd.nextLong() ) +  ( rnd.nextLong() ^ rnd.nextLong() );
		//System.err.println(idToName.get(i) +" "+hash1[i] +" "+ hash2[i]);
	}
    }
    
    /**
     * Returns the ID corresponding to the taxon name , if the taxon name is new it adds it unless the taxonidentifier is locked
     * @param name
     * @return
     */
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
    
	public  STITreeCluster getClusterForNodeName(String nodeName) {
		STITreeCluster cluster = this.newCluster();;
		Integer taxonID = this.taxonId(nodeName);
		cluster.addLeaf(taxonID);
		return cluster;
	}
}
