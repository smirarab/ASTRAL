package phylonet.coalescent;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import phylonet.coalescent.DuplicationWeightCounter.STBipartition;
import phylonet.coalescent.MGDInference_DP.Vertex;
import phylonet.tree.model.sti.STITreeCluster;
public interface ClusterCollection {

	Vertex getTopCluster();

	int getClusterCount();

	boolean addCluster(Vertex nv, int size);

	boolean contains(Vertex reverse);

	ClusterCollection getContainedClusters(STITreeCluster _cluster);

	Vertex getVertexForCluster(STITreeCluster cluster1);

	Collection<STBipartition> getClusterResolutions(STITreeCluster cluster);


	Iterator<Set<Vertex>> getSubClusters();

}
