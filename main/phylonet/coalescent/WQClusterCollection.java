package phylonet.coalescent;

import java.util.HashMap;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;


public class WQClusterCollection extends BasicClusterCollection{

	public WQClusterCollection(int len,  HashMap<STITreeCluster, Vertex> globalVertexCash) {
		initialize(len, globalVertexCash);
		//geneTreeSTBByCluster = new HashMap<STITreeCluster, HashSet<STBipartition>>();
	}
		

}
