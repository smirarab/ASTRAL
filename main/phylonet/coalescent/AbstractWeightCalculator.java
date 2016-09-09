package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculator<T> implements Cloneable {
	
	private static final boolean TESTRUN = false;
	private int callcounter = 0;
	protected HashMap<T, Long> weights;
	boolean save;
	long lastTime;

	public AbstractWeightCalculator(boolean save) {
		this.save = save;
		this.lastTime = System.currentTimeMillis();
	}
	
	public void initializeWeightContainer(int size) {
		weights = new HashMap<T, Long>(size);
	}
	
	public int getCalculatedWeightCount() {
		//return this.callcounter;	
		if (!save)
			return this.callcounter;
		else 
			return weights.size();
	}
	
	public Long getCalculatedWeight(T t) {
		return weights.get(t);
	}
	
		
	public Long getWeight(T t, AbstractComputeMinCostTask<T> minCostTask) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight == null) {
			//ICalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
			//prepareWeightTask(weigthWork, minCostTask);
			// MP_VERSION: smallWork.fork();
			weight = TESTRUN ? 0 : calculateWeight(t, minCostTask);
			int count;
			if (save && !TESTRUN ) {
				weights.put(t, weight);
				count = weights.size();
			} else {
				count = this.callcounter;
			}
			if (count % 100000 == 0) {
				System.err.println("Calculated "+ count +" weights; time (seconds): " + (System.currentTimeMillis() - lastTime)/1000);
				lastTime = System.currentTimeMillis();
			}
/*			if (weights.size() == 75318) {
				System.err.println("here");
			}*/
		} else {
			//System.err.println("Found " + t );
		}
		return weight;
	}
	
	abstract Long calculateWeight(T t, AbstractComputeMinCostTask<T> minCostTask);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
	
		@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
	
}