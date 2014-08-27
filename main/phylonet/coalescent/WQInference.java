package phylonet.coalescent;


import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInference extends AbstractInference<Tripartition> {
	
	int forceAlg = -1;
	
	public WQInference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution, boolean duploss, int alg, boolean addExtra) {
		super(rooted, extrarooted, trees, extraTrees, exactSolution, addExtra);
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
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(AbstractInference<Tripartition> dlInference,
			Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount(), new HashMap<STITreeCluster, STITreeCluster.Vertex>());
	}
	
	AbstractDataCollection<Tripartition> newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this.forceAlg);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}

}
