package phylonet.coalescent;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import phylonet.tree.model.Tree;
import phylonet.tree.util.Trees;

/**
 * Maps gene names to species names (as strings)
 * 
 * @author smirarab
 *
 */
public class TaxonNameMap {
    private Map<String, String> taxonMap = null;
    private String pattern = null;
    private String rep = null;
    private SpeciesMapper speciesIdMapper;
    //Map<String, BitSet> speciesBitSet = null;

    public TaxonNameMap (Map<String, String> taxonMap) {
        this.taxonMap = taxonMap;
        this.initializeSpeciesMapper();
    }
    public TaxonNameMap (String pattern, String rep) {
        this.pattern = pattern;
        this.rep = rep;
        throw new RuntimeException("Not implemented yet");
    }
    public TaxonNameMap () {
        this.initializeSpeciesMapper();
    }

    private void initializeSpeciesMapper() {
    	this.speciesIdMapper = newSpeciesMapper();
        if (this.taxonMap != null) {
            for (Entry<String, String> entry: this.taxonMap.entrySet()) {
                speciesIdMapper.setSpeciesIdForTaxon(entry.getKey(), entry.getValue());
            }
        } else {
            for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++) {
                speciesIdMapper.setSpeciesIdForTaxon(i,GlobalMaps.taxonIdentifier.getTaxonName(i));
            }
        }
    }
	protected SpeciesMapper newSpeciesMapper() {
		return Factory.instance.newSpeciesMapper(); //GlobalMaps.taxonIdentifier.taxonCount());
	}

    /**
     * For a given gene name, give the species name
     * @param geneName
     * @return
     */
    public String getTaxonName(String geneName) {
        if (geneName == null || "".equals(geneName)) {
            throw new RuntimeException("Empty name?");
        }
        if (this.pattern != null) {
            String s = geneName.replaceAll(pattern,rep);
            //System.err.println("Mapped " + geneName + " to " + s);
            return s;
        } else if (this.taxonMap != null) {
            return taxonMap.get(geneName);
        } else {
            return geneName;
        }
    }

    /**
     * Returns SpeciesIDMapper instance associated with this mapping. 
     * @return
     */
    public SpeciesMapper getSpeciesIdMapper() {
        return this.speciesIdMapper;
    }
    
    public void checkMapping(List<Tree> trees) {
        if (this.taxonMap != null) {
            Map<String,String> taxonMap = GlobalMaps.taxonNameMap.taxonMap;
            String error = Trees.checkMapping(trees, taxonMap);
            if (error != null) {
                throw new RuntimeException("Gene trees have a leaf named "
                        + error
                        + " that hasn't been defined in the mapping file");
            }
        } 

    }

    /*		public BitSet getTaxonBitSet(String speciesName) {
		    if (speciesName == null || "".equals(speciesName)) {
                throw new RuntimeException("Empty name?");
            }
	        return this.speciesBitSet.get(speciesName);
	    }*/
}