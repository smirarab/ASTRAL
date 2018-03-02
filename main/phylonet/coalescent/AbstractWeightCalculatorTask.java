package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculatorTask<T> {

	protected int callcounter = 0;
	protected HashMap<T, Long> weights;
	protected long lastTime;
	int counter;

	public AbstractWeightCalculatorTask() {
		super();
	}

	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Long>(size);
	}

	public int getCalculatedWeightCount() {
			return this.callcounter;
	}

	public Long getCalculatedWeight(T t) {
		return weights.get(t);
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
	
	protected abstract Long calculateWeight(T t, AbstractComputeMinCostTask<T> minCostTask);

	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);

	public abstract  Long getWeight(T t, AbstractComputeMinCostTask<T> minCostTask);
}