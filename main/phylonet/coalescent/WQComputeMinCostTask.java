package phylonet.coalescent;

import java.util.List;

import phylonet.coalescent.ClusterCollection;
import phylonet.coalescent.ComputeMinCostTask;
import phylonet.coalescent.DPInference;
import phylonet.coalescent.Counter.CalculateWeightTask;
import phylonet.coalescent.DLWeightCounter.DPWeightTask;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQComputeMinCostTask extends ComputeMinCostTask<Tripartition>{

	WQWeightCounter wqCounter;
	
	public WQComputeMinCostTask(DPInference<Tripartition> inference, Vertex v,
			ClusterCollection clusters) {
		super(inference, v, clusters);
		this.wqCounter = (WQWeightCounter)inference.counter;
	}
	
	protected double adjustWeight(int clusterLevelCost, Vertex smallV,
			Vertex bigv, Integer Wdom) {	
		return Wdom;
	}
	
	@Override
	protected int scoreBaseCase(boolean rooted, List<Tree> trees) {	
		return 0;
	}

	@Override
	protected ComputeMinCostTask<Tripartition> newMinCostTask(
			 Vertex v, ClusterCollection clusters) {
		return new WQComputeMinCostTask(inference, v, clusters);
	}
	
	@Override
	protected int calculateClusterLevelCost() {
		return 0;
	}
	
	@Override
	protected void initializeWeightTask(CalculateWeightTask weigthWork) {
	}

	@Override
	protected Tripartition STB2T(STBipartition stb) {
		return new Tripartition(stb.cluster1, stb.cluster2, stb.c.complementaryCluster());
	}

}
