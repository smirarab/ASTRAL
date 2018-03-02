package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculator<T> implements Cloneable {
	
	
	private static final boolean TESTRUN = false;
	private int callcounter = 0;
	protected HashMap<T, Long> weights;
	boolean save;
	long lastTime;
	LinkedBlockingQueue<Long> queue;
	boolean done = false;
	int counter;
	public AbstractWeightCalculator(boolean save, LinkedBlockingQueue<Long> queue) {
		this.save = save;
		this.queue = queue;
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

			if(!TESTRUN){
				if(done) {
					
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
					return weight;
				}
				else {
					try {
						weight = queue.take();
					}
					catch(Exception e) {
						e.printStackTrace();
					}
					if(weight == -23) {//random number from CommandLine used as a "poison pill"
						done = true;
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
						return weight;
					}
				}
			}
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

		} else {
			//System.err.println("Found " + t );
		}
		//System.out.println(t.toString());
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