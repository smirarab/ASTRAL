package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class GlobalMaps{
	
	public static class TaxonNameMap {
		Map<String, String> taxonMap;
		String pattern = null;
		String rep = null;
		public TaxonNameMap (Map<String, String> taxonMap) {
			this.taxonMap = taxonMap;
		}
		public TaxonNameMap (String pattern, String rep) {
			this.pattern = pattern;
			this.rep = rep;
		}
		private String getTaxonName(String geneName) {
			if (geneName == null || "".equals(geneName)) {
				throw new RuntimeException("Empty name?");
			}
			if (pattern != null) {
				String s = geneName.replaceAll(pattern,rep);
				//System.err.println("Mapped " + geneName + " to " + s);
				return s;
			} else {
				return taxonMap.get(geneName);
			}
		}	
		
	}
	
    protected static String getSpeciesName(String geneName) {
        if (GlobalMaps.taxonNameMap != null) {
            return  GlobalMaps.taxonNameMap.getTaxonName(geneName);
        }
        return geneName;
    }
    
	public static class TaxonIdentifier {
		private HashMap<String, Integer> nameToId = new HashMap<String, Integer>();
		private List<String> idToName = new ArrayList<String>();
		private int taxonCount = 0;
		
		public String[] getAllTaxonNames(){
		    return idToName.toArray(new String[]{});
		}
		
		public Integer taxonId(String name) {
			Integer a = nameToId.get(name);
			if (a ==  null){
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
	
	public static TaxonIdentifier taxonIdentifier = new TaxonIdentifier();
	public static TaxonNameMap taxonNameMap;
}