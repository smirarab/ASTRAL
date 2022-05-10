package phylonet.coalescent;

import java.util.HashSet;

import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQClusterCollection extends AbstractClusterCollection{

	public WQClusterCollection(int len) {
		initialize(len);
	}

	@Override
	public AbstractClusterCollection newInstance(int size) {
		return new WQClusterCollection(size);
	}

	public void printDiff(WQClusterCollection other) {
		for (int i = 0; i< this.clusters.size(); i++) {
			HashSet<Vertex> temp = new HashSet<Vertex>(this.clusters.get(i));
			temp.retainAll(other.clusters.get(i));
			System.err.println(temp);
		}
	}

}
