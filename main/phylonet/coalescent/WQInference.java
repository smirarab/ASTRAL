package phylonet.coalescent;


import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInference extends Inference<Tripartition> {
	
	int forceAlg = -1;
	
	public WQInference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution, boolean duploss, int alg) {
		super(rooted, extrarooted, trees, extraTrees, exactSolution);
		this.forceAlg = alg;
	}

	public void scoreGeneTree(Tree st) {
		throw new RuntimeException("Not implemented yet");
	}


	@Override
	Long getTotalCost(Vertex all) {
		return (long) (all._max_score/4);
	}


	@Override
	ComputeMinCostTask<Tripartition> newComputeMinCostTask(Inference<Tripartition> dlInference,
			Vertex all, ClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	ClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount(), new HashMap<STITreeCluster, STITreeCluster.Vertex>());
	}
	
	DataCollection<Tripartition> newCounter(ClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this.forceAlg);
	}



	@Override
	WeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}

}
