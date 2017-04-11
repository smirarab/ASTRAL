package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculatorNoCalculations<T> {
	private int callcounter = 0;
	protected HashMap<T, Long> weights;
	boolean save;
	long lastTime;
	LinkedBlockingQueue<T> queue;
	public AbstractWeightCalculatorNoCalculations(boolean save, LinkedBlockingQueue<T> queue1) {
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
	
		
	public Long getWeight(T t, AbstractComputeMinCostTaskNoCalculations<T> minCostTask) {
		this.callcounter ++;
			
			try {
				queue.put(t);
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			
			if (this.callcounter % 100000 == 0) {
				System.err.println(callcounter + " weights produced : " + (double)(System.currentTimeMillis() - lastTime)/1000 + " seconds");
			}
			//System.out.println(t.toString());
		return 0L;
	}
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}
