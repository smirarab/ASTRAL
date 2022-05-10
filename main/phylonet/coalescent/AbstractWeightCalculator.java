package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculator<T> implements Cloneable {
	
	protected static final boolean TESTRUN = false;
	protected int callcounter = 0;
	protected HashMap<T, Long> weights;
	boolean save;
	public long lastTime;

	public AbstractWeightCalculator(boolean save) {
		this.save = save;
		this.lastTime = System.currentTimeMillis();
	}
	
	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Long>(size);
	}
	
	public int getCalculatedWeightCount() {
		if (!save)
			return this.callcounter;
		else 
			return weights.size();
	}
	
	public Long getCalculatedWeight(T t) {
		return weights.get(t);
	}
	
	public Long getWeight(T t, AbstractComputeMinCostTask<T> minCostTask) {
		return getWeight(t);
	}
		
	public Long getWeight(T t) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight != null) {
			return weight;
		}

		weight = TESTRUN ? 0 : calculateWeight(t);
		saveWeight(t, weight);

		return weight;
	}
	
	protected void saveWeight(T t, Long weight) {
		int count;
		if (save && !TESTRUN) {
			weights.put(t, weight);
			count = weights.size();
		} else {
			count = this.callcounter;
		}
		if (count % 100000 == 0) {
			Logging.log("Calculated "+ count +" weights; time (seconds): " + (System.currentTimeMillis() - lastTime)/1000);
			lastTime = System.currentTimeMillis();
		}

	}

	/*public Long[] calculateWeight(T[] ts) {
		Long [] ret = new Long[ts.length];
		int i = 0;
		for (T t: ts)
			ret[i++] = this.calculateWeight(t);
		return ret;
	}*/
	
	public abstract Long[] calculateWeight(T[] ts);
	
	protected abstract Long calculateWeight(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
	
		@Override
		public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}
		
	public abstract void setupGeneTrees(WQInference wqInference);
	
}