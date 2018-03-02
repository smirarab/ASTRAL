package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.Tree;

class WQWeightCalculatorNoCalculations extends AbstractWeightCalculatorProducer<Tripartition> {
	AbstractInference<Tripartition> inference;

	public WQWeightCalculatorNoCalculations(AbstractInference<Tripartition> inference, LinkedBlockingQueue<Tripartition> queue1) {
		super(queue1);
		this.inference =  inference;	
	}

	@Override
	public void preCalculateWeights(List trees, List extraTrees) {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected Long calculateWeight(Tripartition t, AbstractComputeMinCostTask<Tripartition> minCostTask) {
		return 0l;
	}


	

}
