package phylonet.coalescent;


import java.util.concurrent.LinkedBlockingQueue;

/**
 * Knows how to compute the score of a given tripartition
 * @author smirarab
 *
 */
class WQWeightCalculatorMP extends WQWeightCalculator{

	private LinkedBlockingQueue<Long> queue;
	private boolean threadingOff = false;

	public WQWeightCalculatorMP(AbstractInference<Tripartition> inference, LinkedBlockingQueue<Long> queue2) {
		super(inference);
		this.queue = queue2;
		this.lastTime = System.currentTimeMillis();
		this.algorithm = new CondensedTraversalWeightCalculator();
		//tmpalgorithm.setupGeneTrees((WQInference) inference);
	}

	@Override
	public int getCalculatedWeightCount() {
		return this.callcounter;
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


	public boolean isThreadingOff() {
		return threadingOff;
	}

	public void setThreadingOff(boolean done) {
		this.threadingOff = done;

	}


}
