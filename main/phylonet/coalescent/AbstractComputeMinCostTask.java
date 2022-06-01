package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;

/**
 * This class implements the dynamic programming
 * @author smirarab
 *
 * @param <T>
 */
public abstract class AbstractComputeMinCostTask<T> {

	protected AbstractInference<T> inference;
	protected SpeciesMapper spm;

	IClusterCollection clusters;
	long target = 0;
	IClusterCollection containedVertecies;

	public AbstractComputeMinCostTask(AbstractInference<T> inference, 
			IClusterCollection clusters) {
		this.inference = inference;
		this.clusters = clusters;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
	}


	protected Long compute() {
		try {
			return (long) computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}

	protected abstract AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
			IClusterCollection clusters);

	protected AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
			IClusterCollection clusters, long target){
		AbstractComputeMinCostTask<T> task = newMinCostTask(v, clusters);
		task.target = target;
		return task;
	}

	protected abstract double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom);

	protected abstract Long defaultWeightForFullClusters();

	protected abstract long calculateClusterLevelCost();

	protected abstract long scoreBaseCase(boolean rooted, List<Tree> trees);

	protected abstract T STB2T(VertexPair stb);

	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
	protected abstract long computeMinCost() throws CannotResolveException;


	protected void addComplementaryClusters(int clusterSize, Vertex v) {
		Iterator<Set<Vertex>> it = containedVertecies.getSubClusters().iterator();
		while (it.hasNext()) {
			Collection<Vertex> subClusters = new ArrayList<Vertex>(it.next());
			int i = -1;
			for (Vertex x : subClusters) {
				i = i > 0 ? i : x.getCluster().getClusterSize();
				int complementarySize = clusterSize - i;
				containedVertecies.addCluster(
						getCompleteryVertx(x, v.getCluster()),
						complementarySize);
			}
			if (i < clusterSize * inference.getCD()) {
				return;
			}

		}
	}

	/***
	 * Used in the exact version
	 * @param cluster
	 * @param containedVertecies
	 */
	void addAllPossibleSubClusters(STITreeCluster cluster,
			IClusterCollection containedVertecies) {
		int size = cluster.getClusterSize();
		for (int i = cluster.getBitSet().nextSetBit(0); i >= 0; i = cluster
				.getBitSet().nextSetBit(i + 1)) {
			STITreeCluster c = Factory.instance.newCluster(cluster);
			c.getBitSet().clear(i);

			Vertex nv = c.newVertex();
			containedVertecies.addCluster(nv, size - 1);

			addAllPossibleSubClusters(c, containedVertecies);
		}
	}

	public Vertex getCompleteryVertx(Vertex x, STITreeCluster refCluster) {
		STITreeCluster c = x.getCluster();

		STITreeCluster revcluster = Factory.instance.newCluster(refCluster);
		revcluster.getBitSet().xor(c.getBitSet());
		Vertex reverse = revcluster.newVertex();
		return reverse;
	}

}
