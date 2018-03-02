package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class HashClusterCollection extends AbstractClusterCollection {

	HashMap<Long, Vertex> map = null;
	
	public HashMap<Long, Vertex> getMap() {
		return map;
	}

	public void setMap(HashMap<Long, Vertex> map) {
		this.map = map;
	}

	public HashClusterCollection(int n) {
		this.initialize(n);
	}
	
	@Override
	public AbstractClusterCollection newInstance(int size) {
		HashClusterCollection newInst = new HashClusterCollection(size);
		if (this.getMap() != null) {
			newInst.setMap(this.getMap());
		}
		return newInst;
		
	}
	
	public void preComputeHashValues() {
		System.err.println("Computing hash values ...");
		Random rnd = GlobalMaps.random;
		int n = GlobalMaps.taxonIdentifier.taxonCount();
		long[] hash1 = new long[n], hash2 = new long[n];
		boolean succeed = false;
		while (!succeed) {
			map = new HashMap<Long, Vertex>();
			succeed = true;
			for (int i = 0; i < n; i++) {
				hash1[i] = rnd.nextLong();
				hash2[i] = rnd.nextLong();
			}
			for (int i = 1; i <= n; i++) {
				for (Vertex v: this.clusters.get(i)) {
					long h1 = 0, h2 = 0;
					BitSet b = v.getCluster().getBitSet();
					
					for (int k = b.nextSetBit(0); k >= 0; k = b.nextSetBit(k + 1)) {
							h1 += hash1[k];
							h2 ^= hash2[k];
					}
					if (map.containsKey(h1)) {
						succeed = false;
						break;
					}
					map.put(h1, v);
					v.hash1 = h1;
					v.hash2 = h2;
					v.clusterSize = i;
				}
				if (!succeed) break;
			}
		}
		System.err.println("Done computing hash values ...");
	}


	@Override
	public getClusterResolutionsLoop getClusterResolutionLoop(int i, Vertex vert, int clusterSize) {
		return new hashClusterResolutionsLoop(this.clusters, i, vert, clusterSize);
	}
	
	public class hashClusterResolutionsLoop extends getClusterResolutionsLoop{

		public hashClusterResolutionsLoop(ArrayList<Set<Vertex>> cluster, int i, Vertex v, int clusterSize) {
			super(cluster, i, v, clusterSize);
		}
		public ArrayList<VertexPair> call() {
			ArrayList<VertexPair> ret = new ArrayList<VertexPair>();

			if (this.clusters.get(i) == null || this.clusters.get(i).size() == 0 ||
				this.clusters.get(clusterSize - i) == null || this.clusters.get(clusterSize - i).size() == 0) {
				return null;
			}
			
			Set<Vertex> small = (this.clusters.get(i).size() < this.clusters.get(clusterSize - i).size()) ? this.clusters.get(i) : this.clusters.get(clusterSize - i);

			for (Vertex v1 : small) {
				Vertex v2 = map.get(v.hash1 - v1.hash1);
				if (v2 != null) {
					if (v.hash2 == (v1.hash2 ^ v2.hash2)
						&& v.clusterSize == v1.clusterSize + v2.clusterSize) {
						if (v1.clusterSize != v2.clusterSize || v1.hash1 < v2.hash1) {
							VertexPair bi = new VertexPair(v1, v2, v);
							ret.add(bi);
						}
					}
				}
			}
			return ret;
		}
	}

	

}
