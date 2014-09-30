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
		//return this.callcounter;	
		System.err.println("Weights requested "+this.callcounter +" times");
		return weights.size();
	}
	
	public Long getCalculatedWeight(T t) {
		return weights.get(t);
	}
	
	public Long getWeight(T t, AbstractComputeMinCostTask<T> minCostTask) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight == null) {
			ICalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
			prepareWeightTask(weigthWork, minCostTask);
			// MP_VERSION: smallWork.fork();
			weight = TESTRUN ? 0 : weigthWork.calculateWeight();
			int count;
			if (! TESTRUN ) {
				weights.put(t, weight);
				count = weights.size();
			} else
				count = this.callcounter;
			if (count % 100000 == 0)
				System.err.println("Calculated "+ count +" weights; time:" + System.currentTimeMillis());
/*			if (weights.size() == 75318) {
				System.err.println("here");
			}*/
		} else {
			//System.err.println("Found " + t );
		}
		return weight;
	}
	
	protected abstract void prepareWeightTask(ICalculateWeightTask<T> weigthWork,AbstractComputeMinCostTask<T> task);
	
	public abstract ICalculateWeightTask<T> getWeightCalculateTask(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}