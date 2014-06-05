package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class WeightCalculator<T> {
	
	protected HashMap<T, Long> weights;

	public WeightCalculator() {
	}
	
	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Long>(size);
	}
	
	public int getCalculatedWeightCount() {
		return weights.size();
	}
	
	public Long getCalculatedWeight(T bi) {
//		if (!weights.containsKey(bi)) {
//			// weights.put(bi,calculateMissingWeight(bi));
//			return null;
//		}
		return weights.get(bi);
	}
	
	public Long getWeight(T t, ComputeMinCostTask<T> minCostTask) {
		Long weight = getCalculatedWeight(t);
		
		if (weight == null) {
//						if (clusterSize > 9 && v._max_score > (2-0.02*clusterSize)*(lscore + rscore))
//							continue;
			CalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
			prepareWeightTask(weigthWork, minCostTask);
			// MP_VERSION: smallWork.fork();
			weight = weigthWork.calculateWeight();
		}
		return weight;
	}
	
	protected abstract void prepareWeightTask(CalculateWeightTask<T> weigthWork,ComputeMinCostTask<T> task);
	
	public abstract CalculateWeightTask<T> getWeightCalculateTask(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}