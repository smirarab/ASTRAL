package phylonet.coalescent;

import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster.Vertex;

public interface IClusterCollection {

	Vertex getTopVertex();

	int getClusterCount();

	boolean addCluster(Vertex nv, int size);

	boolean contains(Vertex reverse);

	IClusterCollection getContainedClusters(Vertex v);

	Iterable<VertexPair> getClusterResolutions();

	Iterable<Set<Vertex>> getSubClusters();
	
	Set<Vertex> getSubClusters(int size);
	
	class VertexPair {
		public final Vertex cluster1;
		public final Vertex cluster2;
		public final Vertex both;
		public VertexPair(Vertex c1, Vertex c2, Vertex b) {
			cluster1 = c1;
			cluster2 = c2;
			both = b;
		}
	}

}
