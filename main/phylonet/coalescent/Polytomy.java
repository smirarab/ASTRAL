package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;

public class Polytomy extends AbstractPartition {
	
	STITreeCluster[] clusters;
	private int _hash = 0;
	
	public Polytomy(STITreeCluster[] cs){
		clusters = new STITreeCluster[cs.length];
		long[] ns = new long[cs.length];
		for (int i = 0; i < cs.length; i++){
			cs[i].updateHash();
			clusters[i] = new STITreeCluster(cs[i]);
			ns[i] = cs[i].hash1;
		}
		for (int i = 0; i < cs.length; i++){
			for (int j = i + 1; j < cs.length; j++){
				if (ns[i] < ns[j]){
					long tn = ns[i];
					ns[i] = ns[j];
					ns[j] = tn;
					STITreeCluster tc = clusters[i];
					clusters[i] = clusters[j];
					clusters[j] = tc;
				}
			}
		}
	}
	
	public STITreeCluster[] getClusters(){
		return clusters;
	}
	
	@Override
	public boolean equals(Object obj) {
		if ((obj instanceof Polytomy) == false) return false;
		Polytomy p = (Polytomy) obj; 
		
		if (this == obj) return true;
		for (int i = 0; i < clusters.length; i++){
			if (clusters[i].equals(p.clusters[i]) == false) return false;
		}
		return true;					
	}
	@Override
	public int hashCode() {
		if (_hash == 0) {
			for (int i = 0; i < clusters.length; i++){
				_hash = _hash * 33 + clusters[i].hashCode();
			}
		}
		return _hash;
	}
	@Override
	public String toString() {
		String s = "Polytomy:";
		for (int i = 0; i < clusters.length; i++){
			s += "|" + clusters[i].toString();
		}
		return s;
	}
}
