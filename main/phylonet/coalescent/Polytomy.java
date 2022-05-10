package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;

public class Polytomy<T extends STITreeCluster> extends AbstractPartition<T> {

	protected T[] clusters;
	protected int _hash = 0;

	public void initalize(T[] cs){
		clusters = (T[]) new STITreeCluster[cs.length];
		long[] ns = new long[cs.length];
		for (int i = 0; i < cs.length; i++){
			cs[i].updateHash();
			clusters[i] = (T) new STITreeCluster(cs[i]);
			ns[i] = cs[i].partionId();
		}
		arrange(cs, ns);
	}

	public boolean equals(Object obj) {
		if ((obj instanceof Polytomy) == false) return false;
		Polytomy p = (Polytomy) obj; 

		if (this == obj) return true;
		for (int i = 0; i < clusters.length; i++){
			if (clusters[i].equals(p.clusters[i]) == false) return false;
		}
		return true;					
	}

	protected void arrange(T[] cs, long[] ns) {
		for (int i = 0; i < cs.length; i++){
			for (int j = i + 1; j < cs.length; j++){
				if (ns[i] < ns[j]){
					long tn = ns[i];
					ns[i] = ns[j];
					ns[j] = tn;
					T tc = clusters[i];
					clusters[i] = clusters[j];
					clusters[j] = tc;
				}
			}
		}
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
