package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.Tree;

class WQWeightCalculatorProducer extends AbstractWeightCalculatorProducer<Tripartition> {

	public WQWeightCalculatorProducer(AbstractInference<Tripartition> inference, LinkedBlockingQueue<Tripartition> queue1) {
		super(queue1);	
	}

	@Override
	public void preCalculateWeights(List trees, List extraTrees) {
		
	}

	@Override
	protected Long calculateWeight(Tripartition t, AbstractComputeMinCostTask<Tripartition> minCostTask) {
		return 0l;
	}

}
