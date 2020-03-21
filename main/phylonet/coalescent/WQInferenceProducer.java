package phylonet.coalescent;
import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInferenceProducer extends AbstractInferenceProducer<Tripartition> {
	
	public WQInferenceProducer(Options inOptions, List<Tree> trees, List<Tree> extraTrees,List<Tree> toRemoveExtraTrees) {
		super(inOptions, trees, extraTrees, toRemoveExtraTrees);
	}

	public WQInferenceProducer(AbstractInference semiDeepCopy) {
		super(semiDeepCopy);
	}

	@Override
	Long getTotalCost(Vertex all) {
		return 0l;
	}

	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(
			AbstractInference<Tripartition> inference, Vertex all) {
		return new ComputeMinCostTaskProducer((AbstractInferenceProducer<Tripartition>) inference, all);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}

	WQDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this);
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

	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return null;
	}

}
