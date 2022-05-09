package phylonet.coalescent;


import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;

/**
 * Knows how to compute the score of a given tripartition
 * @author smirarab
 *
 */
class WQWeightCalculatorMP extends WQWeightCalculator{

	private TraversalWeightCalculator tmpalgorithm;
	
	private LinkedBlockingQueue<Long> queue;
	private boolean threadingOff = false;

	public WQWeightCalculatorMP(AbstractInference<Tripartition> inference, LinkedBlockingQueue<Long> queue2) {
		super(inference);
		this.queue = queue2;
		this.lastTime = System.currentTimeMillis();
		this.algorithm = new CondensedTraversalWeightCalculator();
		tmpalgorithm = new TraversalWeightCalculator();
		//tmpalgorithm.setupGeneTrees((WQInference) inference);
	}
	
	@Override
	public int getCalculatedWeightCount() {
		return this.callcounter;
	}
	
	@Override
	public Long calculateWeight(Tripartition t) {
		return this.calculateWeight(new Tripartition[] {(Tripartition) t})[0];
	}

	@Override
	public Long getWeight(Tripartition t) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight != null) {
			return weight;
		}

		if(isThreadingOff()) { // After the main DP, for computing final score. 		
			weight =  calculateWeight(t);
			saveWeight(t, weight);
			return weight;
		}

		try {
			weight = queue.take();
		}
		catch(Exception e) {
			throw new RuntimeException(e);
		}
		if(weight == TurnTaskToScores.THEEND) {// After the main DP, for computing final score, we switch to local calculator
			setThreadingOff(true);
			weight =  calculateWeight(t);
		}
		saveWeight(t, weight);

		return weight;

	}

	/**
	 * one of ASTRAL-III way of calculating weights
	 * Should be memory efficient
	 * @author chaoszhang
	 *
	 */
	class CondensedTraversalWeightCalculator extends WeightCalculatorAlgorithm {
		Polytree polytree;

		@Override
		Long[] calculateWeight(Tripartition[] trip) {
			return polytree.WQWeightByTraversal(trip);
		}
		
		@Override
		Long calculateWeight(Tripartition t) {
			return polytree.WQWeightByTraversal(t);
		}

		@Override
		void setupGeneTrees(WQInference inference) {
			polytree = new Polytree(inference.trees, (WQDataCollectionMP) dataCollection);
		}
	}

	/**
	 * Each algorithm will have its own data structure for gene trees
	 * @param wqInference
	 */
	@Override
	public void setupGeneTrees(WQInference wqInference) {
		tmpalgorithm.setupGeneTrees((WQInferenceConsumerMP) wqInference);
		this.algorithm.setupGeneTrees(wqInference);
	}

	//TODO: this is algorithm-specific should not be exposed. Fix.
	public int[] geneTreesAsInts() {
		return (tmpalgorithm).geneTreesAsInts;
	}
	@Override
	public Long[] calculateWeight(Tripartition[] t) {
		return this.algorithm.calculateWeight(t);
	}

	//TODO: this is algorithm-specific should not be exposed. Fix. 
	public int maxHeight() {
		return ((TraversalWeightCalculator)tmpalgorithm).maxHeight;
	}


	public boolean isThreadingOff() {
		return threadingOff;
	}

	public void setThreadingOff(boolean done) {
		this.threadingOff = done;

	}


}
