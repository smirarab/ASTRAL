package phylonet.coalescent;


import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInference extends Inference<Tripartition> {
	
	public WQInference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution, boolean duploss) {
		super(rooted, extrarooted, trees, extraTrees, exactSolution);
	}

	

	public void scoreGeneTree(STITree st) {
		throw new RuntimeException("Not implemented yet");
	}


	@Override
	int getTotalCost(Vertex all) {
		return (int) (all._max_score/4);
	}


	@Override
	ComputeMinCostTask<Tripartition> newComputeMinCostTask(Inference<Tripartition> dlInference,
			Vertex all, ClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	ClusterCollection newClusterCollection() {
		return new WQClusterCollection(stTaxa.length);
	}
	
	Counter<Tripartition> newCounter(ClusterCollection clusters) {
		return new WQWeightCounter(gtTaxa, stTaxa, rooted, (WQClusterCollection)clusters);
	}

}
