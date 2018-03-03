package phylonet.coalescent;

import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
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
	protected Long[] calculateWeight(Tripartition[] t) {
		Long[] ret = new Long[t.length];
		for (int i = 0; i < ret.length; i++) {
			ret[i] = 0l;
		}
		return ret;
		
	}

}
