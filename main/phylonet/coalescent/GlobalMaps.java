package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;


public class GlobalMaps{
	
	public static class TaxonNameMap {
		Map<String, String> taxonMap = null;
		String pattern = null;
		String rep = null;
		SpeciesMapper speciesIdMapper;
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
		    
		}
		
		public void initializeSpeciesMapper() {
		    speciesIdMapper = new SpeciesMapper(GlobalMaps.taxonIdentifier.taxonCount());
		    if (this.taxonMap != null) {
    		    for (Entry<String, String> entry: this.taxonMap.entrySet()) {
    		        speciesIdMapper.setSpeciesIdForTaxon(entry.getKey(), entry.getValue());
    		    }
		    } else {
		        for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++) {
		            speciesIdMapper.setSpeciesIdForTaxon(i,i);
                }
		    }
		}
		
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
		
		public SpeciesMapper getSpeciesIdMapper() {
		    return this.speciesIdMapper;
		}
		
/*		public BitSet getTaxonBitSet(String speciesName) {
		    if (speciesName == null || "".equals(speciesName)) {
                throw new RuntimeException("Empty name?");
            }
	        return this.speciesBitSet.get(speciesName);
	    }*/
	}


    
	public static class TaxonIdentifier {
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
		
/*		public List<String> getTaxonList() {
			return idToName;
		}*/
	}
	
	public static class SpeciesMapper {
	    
	    private int [] taxonIdToSpeciesId;
	    TaxonIdentifier speciesNameIdMap;
	    
	    public SpeciesMapper(int taxonCount) {
	        this.taxonIdToSpeciesId = new int[taxonCount];
	        this.speciesNameIdMap = new TaxonIdentifier();
	    }
	    
	    public int getSpeciesIdForTaxon(int id) {
	        return this.taxonIdToSpeciesId[id];
	    }
	    
	    protected void setSpeciesIdForTaxon(int taxonId, int speciesId) {
	        this.taxonIdToSpeciesId[taxonId] = speciesId;
            //System.err.println("Mapped taxon "+taxonId +" to species "+speciesId);
	    }
	    
	    protected void setSpeciesIdForTaxon(int taxonId, String speciesName) {
            this.setSpeciesIdForTaxon(taxonId, this.speciesNameIdMap.taxonId(speciesName));
        }
	    
       protected void setSpeciesIdForTaxon(String taxonName, String speciesName) {
           this.setSpeciesIdForTaxon(GlobalMaps.taxonIdentifier.taxonId(taxonName),
                   this.speciesNameIdMap.taxonId(speciesName));
       }
       
       public int getSpeciesCount() {
           return this.speciesNameIdMap.taxonCount();
       }
       
       public String getSpeciesName(int stId) {
           return this.speciesNameIdMap.getTaxonName(stId);
       }
	    
	}
	
	public static TaxonIdentifier taxonIdentifier = new TaxonIdentifier();
	public static TaxonNameMap taxonNameMap;
}