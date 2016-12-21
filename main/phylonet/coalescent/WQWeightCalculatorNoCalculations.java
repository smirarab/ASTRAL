package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;

import phylonet.tree.model.Tree;

class WQWeightCalculatorNoCalculations extends AbstractWeightCalculatorNoCalculations<Tripartition> {
	public static boolean HAS_NOT = true;
	public static boolean WRITE_OR_DEBUG = false;
	AbstractInference<Tripartition> inference;
	private WQDataCollection dataCollection;

	public WQWeightCalculatorNoCalculations(AbstractInference<Tripartition> inference, ConcurrentLinkedQueue<Tripartition> queue1) {
		super(false, queue1);
		this.dataCollection = (WQDataCollection) inference.dataCollection;
		if(inference instanceof AbstractInferenceNoCalculations) {
			this.inference = (WQInferenceNoCalculations) inference;
		}
		else
			this.inference = (WQInference) inference;	
		}

	@Override
	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {
		// TODO Auto-generated method stub
		
	}

	

}
