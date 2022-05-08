package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeClusterMP;

public class PolytomyMP extends Polytomy<STITreeClusterMP> {
	
	
	@Override
	public void initalize(STITreeClusterMP[] cs){
		clusters = new STITreeClusterMP[cs.length];
		long[] ns = new long[cs.length];
		for (int i = 0; i < cs.length; i++){
			cs[i].updateHash();
			clusters[i] = cs[i].clone();
			ns[i] = cs[i].hash1;
		}
		arrange(cs, ns);
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
