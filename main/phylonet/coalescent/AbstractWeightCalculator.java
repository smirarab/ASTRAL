package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculator<T> {
	
	private static final boolean TESTRUN = false;
	private int callcounter = 0;
	protected HashMap<T, Long> weights;

	public AbstractWeightCalculator() {
	}
	
	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Long>(size);
	}
	
	public int getCalculatedWeightCount() {
		System.err.println("Weights requested "+this.callcounter +" times");
		return weights.size();
	}
	
	public Long getCalculatedWeight(T bi) {
//		if (!weights.containsKey(bi)) {
//			// weights.put(bi,calculateMissingWeight(bi));
//			return null;
//		}
		return weights.get(bi);
	}
	
	public Long getWeight(T t, AbstractComputeMinCostTask<T> minCostTask) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight == null) {
//			if (clusterSize > 9 && v._max_score > (2-0.02*clusterSize)*(lscore + rscore))
//			continue;
			ICalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
			prepareWeightTask(weigthWork, minCostTask);
			// MP_VERSION: smallWork.fork();
			weight = TESTRUN ? 0 : weigthWork.calculateWeight();
			weights.put(t, weight);
		}
		return weight;
	}
	
	protected abstract void prepareWeightTask(ICalculateWeightTask<T> weigthWork,AbstractComputeMinCostTask<T> task);
	
	public abstract ICalculateWeightTask<T> getWeightCalculateTask(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}