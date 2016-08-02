package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculatorNoCalculations<T> {
	private int callcounter = 0;
	protected HashMap<T, Long> weights;
	boolean save;
	long lastTime;
	ConcurrentLinkedQueue<ICalculateWeightTask<T>> queue;
	
	public AbstractWeightCalculatorNoCalculations(boolean save, ConcurrentLinkedQueue<ICalculateWeightTask<T>> queue1) {
		this.save = save;
		this.queue = queue1;
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
			ICalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
			prepareWeightTask(weigthWork, minCostTask);
			// MP_VERSION: smallWork.fork();
			queue.add(weigthWork);
			
			
/*			if (weights.size() == 75318) {
				System.err.println("here");
			}*/
			//System.out.println(t.toString());
		return 0L;
	}
	protected abstract void prepareWeightTask(ICalculateWeightTask<T> weigthWork,AbstractComputeMinCostTask<T> task);
	
	public abstract ICalculateWeightTask<T> getWeightCalculateTask(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}
