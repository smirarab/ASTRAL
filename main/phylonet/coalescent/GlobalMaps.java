package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import phylonet.util.BitSet;


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
	    private List<Set<Integer>> speciesIdtoTaxonId;
	    TaxonIdentifier speciesNameIdMap;
	    
	    public SpeciesMapper(int taxonCount) {
	        this.taxonIdToSpeciesId = new int[taxonCount];
	        this.speciesNameIdMap = new TaxonIdentifier();
	        this.speciesIdtoTaxonId = new ArrayList<Set<Integer>>();
	    }
	    
	    public TaxonIdentifier getSTTaxonIdentifier() {
	        return this.speciesNameIdMap;
	    }
	    
	    public int getSpeciesIdForTaxon(int id) {
	        return this.taxonIdToSpeciesId[id];
	    }
	    
	    private void setSpeciesIdForTaxon(int taxonId, int speciesId) {
	        this.taxonIdToSpeciesId[taxonId] = speciesId;
	        for (int i = this.speciesIdtoTaxonId.size(); i <= speciesId; i++) {
                this.speciesIdtoTaxonId.add(new TreeSet<Integer>());
            }
	        this.speciesIdtoTaxonId.get(speciesId).add(taxonId);
            //System.err.println("Mapped taxon "+taxonId +" to species "+speciesId);
	    }
	    
	    protected void setSpeciesIdForTaxon(int taxonId, String speciesName) {
            this.setSpeciesIdForTaxon(taxonId, this.speciesNameIdMap.taxonId(speciesName));
        }
	    
       protected void setSpeciesIdForTaxon(String taxonName, String speciesName) {
           this.setSpeciesIdForTaxon(GlobalMaps.taxonIdentifier.taxonId(taxonName),
                   this.speciesNameIdMap.taxonId(speciesName));
       }
       
       public Set<Integer> getTaxonsForSpecies(Integer species){
           return this.speciesIdtoTaxonId.get(species);
       }
       
       public int getSpeciesCount() {
           return this.speciesNameIdMap.taxonCount();
       }
       
       public String getSpeciesName(int stId) {
           return this.speciesNameIdMap.getTaxonName(stId);
       }
       
       public Integer speciesId(String name) {
           return this.speciesNameIdMap.taxonId(name);
       }
       
       protected BitSet geneBitsetToSTBitSt(BitSet bs) {
           BitSet stbs = new BitSet(this.getSpeciesCount());
           for (int i = bs.nextSetBit(0); i >=0 ; i = bs.nextSetBit(i+1)) {
               stbs.set(this.getSpeciesIdForTaxon(i));
           }
           return stbs;
       }
       
       public void addMissingIndividuals(BitSet geneBS) {
           BitSet stBS = this.geneBitsetToSTBitSt(geneBS);
           for (int i = stBS.nextSetBit(0); i >=0 ; i = stBS.nextSetBit(i+1)) {
               Set<Integer> taxonsForSpecies = this.getTaxonsForSpecies(i);
               for (Iterator<Integer> it = taxonsForSpecies.iterator(); it.hasNext();) {
                Integer gi = it.next();
                geneBS.set(gi);               
            }
           }
       }

    public boolean isSingleSP(BitSet bs) {
        int i = bs.nextSetBit(0);
        int prevId =  this.getSpeciesIdForTaxon(i);
        for (; i >=0 ; i = bs.nextSetBit(i+1)) {
            int stID = this.getSpeciesIdForTaxon(i);
            if (prevId != stID) {
                return false;
            }
            prevId = stID;
        }
        return true;
    }
	    
	}
	
	public static TaxonIdentifier taxonIdentifier = new TaxonIdentifier();
	public static TaxonNameMap taxonNameMap;
}