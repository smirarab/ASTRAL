package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractComputeMinCostTask<T> {

	protected AbstractInference<T> inference;
	protected Vertex v;
	protected SpeciesMapper spm;

	final byte getDoneState = 1;
	final byte getOtherDoneState = 3;

	public AbstractComputeMinCostTask(AbstractInference<T> inference, Vertex v) {
		this.inference = inference;
		this.v = v;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
	}

	protected Double compute() {
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}
	
	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
	abstract double computeMinCost() throws CannotResolveException;

	public Long getWeight(T t) {
		return inference.weightCalculator.getWeight(t);
	}


	abstract Long defaultWeightForFullClusters();

	protected abstract AbstractComputeMinCostTask<T> newMinCostTask(Vertex v);

	protected abstract double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom);

	protected abstract long calculateClusterLevelCost();

	protected abstract long scoreBaseCase(boolean rooted, List<Tree> trees);

	protected abstract T STB2T(VertexPair stb);

	/***
	 * Used in the exact version
	 * @param cluster
	 * @param containedVertecies
	 */
	void addAllPossibleSubClusters(STITreeCluster cluster, IClusterCollection containedVertecies) {
		int size = cluster.getClusterSize();
		for (int i = cluster.getBitSet().nextSetBit(0); i >= 0; i = cluster
				.getBitSet().nextSetBit(i + 1)) {
			STITreeCluster c = new STITreeCluster(cluster);
			c.getBitSet().clear(i);
	
			Vertex nv = c.new Vertex();
			containedVertecies.addCluster(nv, size - 1);
	
			addAllPossibleSubClusters(c, containedVertecies);
		}
	}
	



}
