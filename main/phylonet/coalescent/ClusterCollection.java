package phylonet.coalescent;

import java.util.Collection;
import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public interface ClusterCollection {

	Vertex getTopVertex();

	int getClusterCount();

	boolean addCluster(Vertex nv, int size);

	boolean contains(Vertex reverse);

	ClusterCollection getContainedClusters(STITreeCluster _cluster);

	Vertex getVertexForCluster(STITreeCluster cluster1);

	Collection<STBipartition> getClusterResolutions();

	Iterable<Set<Vertex>> getSubClusters();

	void addGeneTreeSTB(STBipartition stb, int size);

	Iterable<STBipartition> getContainedGeneTreeSTBs();

}
