package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculator<T> {

	protected int callcounter = 0;
	protected HashMap<T, Long> weights;
	protected long lastTime;
	int counter;

	public AbstractWeightCalculator() {
		super();
	}

	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Long>(size);
	}

	public int getCalculatedWeightCount() {
			return this.callcounter;
	}

	Long getCalculatedWeight(T t) {
		return weights.get(t);
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
	public abstract Long getWeight2(T t);

	protected abstract Long[] calculateWeight(T[] t);

	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);

	public abstract  Long getWeight(T t);
}