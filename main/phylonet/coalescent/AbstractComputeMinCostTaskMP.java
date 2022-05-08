package phylonet.coalescent;

import java.util.List;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeClusterMP.VertexMP;

public abstract class AbstractComputeMinCostTaskMP<T> extends AbstractComputeMinCostTask<T> {

	VertexMP v; 
	public AbstractComputeMinCostTaskMP(AbstractInference<T> inference, VertexMP v) {
		super(inference, inference.dataCollection.clusters);
		this.v = v;
	}

	protected Long compute() {
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}
	
	protected abstract AbstractComputeMinCostTaskMP<T> newMinCostTask(VertexMP v);
	
	@Override
	protected AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
			IClusterCollection clusters){
		return newMinCostTask((VertexMP)v);
	}

	@Override
	protected double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom) {	
		return Wdom;
	}
	

	@Override
	protected long scoreBaseCase(boolean rooted, List<Tree> trees) {	
		return 0l;
	}

	protected abstract T STB2T(VertexPair stb);
	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
	abstract long computeMinCost() throws CannotResolveException;


	public Long getWeight(T t) {
		return inference.weightCalculator.getWeight(t);
	}

}
