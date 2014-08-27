package phylonet.coalescent;

import java.util.HashMap;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;


public class WQClusterCollection extends AbstractClusterCollection{

	public WQClusterCollection(int len,  HashMap<STITreeCluster, Vertex> globalVertexCash) {
		initialize(len, globalVertexCash);
		//geneTreeSTBByCluster = new HashMap<STITreeCluster, HashSet<STBipartition>>();
	}

	@Override
	public AbstractClusterCollection newInstance(int size) {
		return new WQClusterCollection(size, this.globalVertexCash);
	}
		

}
