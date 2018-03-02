package phylonet.coalescent;
import java.util.List;

import phylonet.coalescent.BipartitionWeightCalculator.Quadrapartition;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInferenceNoCalculations extends AbstractInferenceProducer<Tripartition> {
	
	public WQInferenceNoCalculations(Options inOptions, List<Tree> trees, List<Tree> extraTrees) {
		super(inOptions, trees, extraTrees);
	}

	public WQInferenceNoCalculations(AbstractInference semiDeepCopy) {
		super(semiDeepCopy);
	}

	@Override
	Long getTotalCost(Vertex all) {
		return 0l;
	}

	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(
			AbstractInference<Tripartition> inference, Vertex all) {
		return new WQComputeMinCostTaskProducer((AbstractInferenceProducer<Tripartition>) inference, all);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}

	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this);
	}

	AbstractWeightCalculatorProducer<Tripartition> newWeightCalculator() {
		return new WQWeightCalculatorProducer(this, super.queue1);
	}

	@Override
	public double scoreGeneTree(Tree scorest, boolean initialize) {
		return 0;
	}

	@Override
	public double scoreSpeciesTreeWithGTLabels(Tree scorest, boolean initialize) {
		return 0;
	}

	@Override
	void initializeWeightCalculator() {
		this.weightCalculator.initializeWeightContainer(0);

	}

	@Override
	void setupMisc() {

	}

}
