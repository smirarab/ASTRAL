package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class HashClusterCollection extends AbstractClusterCollection {

	HashMap<Long, Vertex> h1ToVertexMap = null;
	
	public HashMap<Long, Vertex> getH1ToVertexMap() {
		return h1ToVertexMap;
	}

	public void setH1ToVertexMap(HashMap<Long, Vertex> map) {
		this.h1ToVertexMap = map;
	}

	public HashClusterCollection(int n) {
		this.initialize(n);
	}
	
	@Override
	public AbstractClusterCollection newInstance(int size) {
		HashClusterCollection newInst = new HashClusterCollection(size);
		if (this.getH1ToVertexMap() != null) {
			newInst.setH1ToVertexMap(this.getH1ToVertexMap());
		}
		return newInst;
		
	}
	
	public void preComputeHashValues() {
		long t = System.currentTimeMillis();
		Random rnd = GlobalMaps.random;
		int n = GlobalMaps.taxonIdentifier.taxonCount();
		long[] hash1 = new long[n], hash2 = new long[n];
		boolean succeed = false;
		while (!succeed) {
			h1ToVertexMap = new HashMap<Long, Vertex>();
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
							h2 += hash2[k];
					}
					if (h1ToVertexMap.containsKey(h1)) {
						succeed = false;
						break;
					}
					h1ToVertexMap.put(h1, v);
					v.hash1 = h1;
					v.hash2 = h2;
					v.clusterSize = i;
				}
				if (!succeed) break;
			}
		}
		System.err.println("Computing hash values took "+((System.currentTimeMillis()-t)/1000)+" seconds");
	}


	@Override
	public getClusterResolutionsLoop getClusterResolutionLoop(int i, Vertex vert, int clusterSize) {
		return new hashClusterResolutionsLoop( i, vert, clusterSize);
	}
	public class hashClusterResolutionsLoop extends getClusterResolutionsLoop{

		public hashClusterResolutionsLoop(int i, Vertex v, int clusterSize) {
			super( i, v, clusterSize);
		}
		public ArrayList<VertexPair> call() {
			ArrayList<VertexPair> ret = new ArrayList<VertexPair>();

			if (clusters.get(i) == null || clusters.get(i).size() == 0 ||
				clusters.get(clusterSize - i) == null || clusters.get(clusterSize - i).size() == 0) {
				return null;
			}
			Set<Vertex> small = (clusters.get(i).size() < clusters.get(clusterSize - i).size()) ? clusters.get(i) : clusters.get(clusterSize - i);

			for (Vertex v1 : small) {
				Vertex v2 = h1ToVertexMap.get(v.hash1 - v1.hash1);
				if (v2 != null) {
					if (v.hash2 == (v1.hash2 + v2.hash2)
						&& v.clusterSize == v1.clusterSize + v2.clusterSize) {
						if (v1.clusterSize != v2.clusterSize || v1.hash1 < v2.hash1) { // To avoid a pair of the same size twice
							VertexPair bi = new VertexPair(v1, v2, v);
							ret.add(bi);
						}
					}
				}
			}
			return ret;
		}
	}
	

	
	public IClusterCollection getContainedClusters(Vertex v) {
		STITreeCluster cluster = v.getCluster();
		int size = cluster.getClusterSize();
		AbstractClusterCollection ret = newInstance(size);
		ret.clusters = this.clusters;
		ret.topV = v;
		ret.totalcount = -1;
		return ret;
	}
	

}
