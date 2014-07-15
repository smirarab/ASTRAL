package phylonet.coalescent;

import java.util.Map;
import java.util.Map.Entry;

public class TaxonNameMap {
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
		    this.initializeSpeciesMapper();
		}
		
		public void initializeSpeciesMapper() {
		    speciesIdMapper = new SpeciesMapper(GlobalMaps.taxonIdentifier.taxonCount());
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