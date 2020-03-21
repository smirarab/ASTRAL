package phylonet.coalescent;

import java.util.HashMap;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex; 

public class HashClusterCollection extends AbstractClusterCollection {

	protected HashMap<Long, Vertex> h1ToVertexMap = null;

	public HashClusterCollection(int n) {
		this.initialize(n);
	}
	
	@Override
	public AbstractClusterCollection newInstance(int size) {
		HashClusterCollection newInst = new HashClusterCollection(size);
		if (this.h1ToVertexMap != null) {
			newInst.h1ToVertexMap = (this.h1ToVertexMap);
		}
		return newInst;
		
	}
	
	public void preComputeHashValues() {
		long t = System.currentTimeMillis();
		int n = GlobalMaps.taxonIdentifier.taxonCount();
		boolean succeed = false;
		while (!succeed) {
			h1ToVertexMap = new HashMap<Long, Vertex>();
			succeed = true;
			for (int i = 1; i <= n; i++) {
				for (Vertex v: this.clusters.get(i)) {
					v.getCluster().updateHash();
					long h1 = v.getCluster().hash1;
					if (h1ToVertexMap.containsKey(h1) || h1 == 0 || h1 == (1L << 63)) { // make sure h1 != -h1
						succeed = false;
						break;
					}
					h1ToVertexMap.put(h1, v);
					v.clusterSize = i;
				}
				if (!succeed) throw new RuntimeException("Bad Random Bits, bad luck. Please rerun the program.");
			}
		}
		System.err.println("Computing hash values took "+((System.currentTimeMillis()-t)/1000)+" seconds");
	}
	

	
	public IClusterCollection getContainedClusters(Vertex v) {
		STITreeCluster cluster = v.getCluster();
		int size = cluster.getClusterSize();
		AbstractClusterCollection ret = newInstance(size);
		ret.clusters = this.clusters;
		ret.topV = v;
		return ret;
	}
	

}
