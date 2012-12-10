package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class STBipartition {
	
	STITreeCluster cluster1;
	STITreeCluster cluster2;		
	STITreeCluster c;
	private int _hash = 0;
	
	public STBipartition(STITreeCluster c1, STITreeCluster c2, STITreeCluster cluster) {
		if (c1.getBitSet().nextSetBit(0) > c2.getBitSet().nextSetBit(0)) {
			cluster1 = c1;
			cluster2 = c2;
		} else {
			cluster1 = c2;
			cluster2 = c1;				
		}
		c = cluster;				
	}
	@Override
	public boolean equals(Object obj) {
		STBipartition stb2 = (STBipartition) obj; 
		
		return this == obj ||
				((stb2.cluster1.equals(this.cluster1) && stb2.cluster2.equals(this.cluster2)));					
	}
	@Override
	public int hashCode() {
		if (_hash == 0) {
			_hash = cluster1.hashCode() * cluster2.hashCode();
		}
		return _hash;
	}
	@Override
	public String toString() {		
		return cluster1.toString()+"|"+cluster2.toString();
	}
	public boolean isDominatedBy(STBipartition dominant) {
/*		if (! dominant.c.containsCluster(this.c)) {
			return false;
		}*/
		//cnt++;
		return (dominant.cluster1.containsCluster(this.cluster1) && dominant.cluster2.containsCluster(this.cluster2)) ||
				(dominant.cluster2.containsCluster(this.cluster1) && dominant.cluster1.containsCluster(this.cluster2));
	}
	
	public boolean isDominatedBy(ClusterCollection left_contained,ClusterCollection rigth_contained, Vertex lv, Vertex rv) {
		/*Vertex thisc = contained.getVertexForCluster(this.c);
		if (! contained.contains(thisc)) {
			return false;
		}*/
		//cnt++;
//		Vertex lv = left_contained.getVertexForCluster(this.cluster1);
	//	Vertex rv = rigth_contained.getVertexForCluster(this.cluster2);
		return (left_contained.contains(lv) && rigth_contained.contains(rv)) ||
				(rigth_contained.contains(lv) && left_contained.contains(rv));
	}
	
	public STBipartition getInducedSTB(STITreeCluster cluster) {
		STITreeCluster lf = new STITreeCluster(c.getTaxa());
		lf.setCluster((BitSet) this.cluster1.getBitSet().clone());
		lf.getBitSet().and(cluster.getBitSet());
		
		STITreeCluster rf = new STITreeCluster(c.getTaxa());
		rf.setCluster((BitSet) this.cluster2.getBitSet().clone());
		rf.getBitSet().and(cluster.getBitSet());

		STITreeCluster cf = new STITreeCluster(c.getTaxa());
		cf.setCluster((BitSet) this.c.getBitSet().clone());
		cf.getBitSet().and(cluster.getBitSet());
		
		return new STBipartition(lf, rf, cf);
	}
	
}
