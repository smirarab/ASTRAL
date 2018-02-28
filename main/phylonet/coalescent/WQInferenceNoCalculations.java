package phylonet.coalescent;
import java.util.List;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInferenceNoCalculations extends AbstractInferenceNoCalculations<Tripartition> {

	int forceAlg = -1;
	long maxpossible;
	public WQInferenceNoCalculations(Options inOptions, List<Tree> trees, List<Tree> extraTrees) {
		super(inOptions, trees, extraTrees);

		this.forceAlg = inOptions.getAlg();
	}


	public WQInferenceNoCalculations(AbstractInference semiDeepCopy) {
		super(semiDeepCopy);
	}




	boolean skipNode (TNode node) {
		return 	node.isLeaf() || node.isRoot() || node.getChildCount() > 2 || 
				(node.getParent().getChildCount() >3) ||
				(node.getParent().getChildCount() >2 && !node.getParent().isRoot())||
				((node.getParent().isRoot() && node.getParent().getChildCount() == 2));
	}

	class NodeData {
		Double mainfreq, alt1freqs, alt2freqs;
		Long quartcount;
		Integer effn ;
		Quadrapartition [] quads;
		STBipartition[] bipartitions;

	}



	@Override
	Long getTotalCost(Vertex all) {
		System.err.println("Normalized score (portion of input quartet trees satisfied): " + 
				all._max_score/4./this.maxpossible);
		return (long) (all._max_score/4l);
	}


	@Override

	AbstractComputeMinCostTaskNoCalculations<Tripartition> newComputeMinCostTaskNoCalculations(AbstractInferenceNoCalculations<Tripartition> dlInference,
			Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTaskNoCalculations( (WQInferenceNoCalculations) dlInference, all, clusters);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}

	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this);
	}


	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this, super.queue2);
	}
	AbstractWeightCalculatorNoCalculations<Tripartition> newWeightCalculatorNoCalculations() {
		return new WQWeightCalculatorNoCalculations(this, super.queue1);
	}


	@Override
	public double scoreGeneTree(Tree scorest, boolean initialize) {
		// TODO Auto-generated method stub
		return 0;
	}


	@Override
	public double scoreSpeciesTreeWithGTLabels(Tree scorest, boolean initialize) {
		// TODO Auto-generated method stub
		return 0;
	}


	@Override
	void initializeWeightCalculator() {
		// TODO Auto-generated method stub

	}


	@Override
	void setupMisc() {
		// TODO Auto-generated method stub

	}


	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(AbstractInference<Tripartition> dlInference,
			Vertex all) {
		// TODO Auto-generated method stub
		return null;
	}

}
