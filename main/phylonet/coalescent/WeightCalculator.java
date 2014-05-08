package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class WeightCalculator<T> {
	
	protected HashMap<T, Integer> weights;

	public WeightCalculator() {
	}
	
	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Integer>(size);
	}
	
	public int getCalculatedWeightCount() {
		return weights.size();
	}
	
	public Integer getCalculatedWeight(T bi) {
//		if (!weights.containsKey(bi)) {
//			// weights.put(bi,calculateMissingWeight(bi));
//			return null;
//		}
		return weights.get(bi);
	}
	
	public abstract CalculateWeightTask<T> getWeightCalculateTask(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}