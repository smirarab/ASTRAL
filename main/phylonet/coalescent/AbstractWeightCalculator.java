package phylonet.coalescent;

import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;

import phylonet.tree.model.Tree;

public abstract class AbstractWeightCalculator<T> implements Cloneable {
	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
	
	private static final boolean TESTRUN = false;
	private int callcounter = 0;
	protected HashMap<T, Long> weights;
	boolean save;
	long lastTime;
	ConcurrentLinkedQueue<Long> queue;
	boolean done = false;
	
	public AbstractWeightCalculator(boolean save, ConcurrentLinkedQueue<Long> queue) {
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
				while(queue.isEmpty()) {	
					try {
						if(!done)
							Thread.sleep(100);
						else {
							ICalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
							prepareWeightTask(weigthWork, minCostTask);
							// MP_VERSION: smallWork.fork();
							weight = TESTRUN ? 0 : weigthWork.calculateWeight();
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
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
					
				weight = queue.remove();
					//System.out.println(t.toString());
	
				
				/*
				ICalculateWeightTask<T> weigthWork = getWeightCalculateTask(t);
				prepareWeightTask(weigthWork, minCostTask);
				Long a = weigthWork.calculateWeight();
				if(!weight.equals(a)) {
					try {
						throw new Exception(weight + " " + a);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						System.err.println("");
						e.printStackTrace();
					}
				}
				*/
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
	
	
	protected abstract void prepareWeightTask(ICalculateWeightTask<T> weigthWork,AbstractComputeMinCostTask<T> task);
	
	public abstract ICalculateWeightTask<T> getWeightCalculateTask(T t);
	
	public abstract void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees);
}