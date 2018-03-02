package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

public abstract class AbstractWeightCalculator<T> extends AbstractWeightCalculatorTask<T> implements Cloneable {
	
	boolean save;
	LinkedBlockingQueue<Long> queue;
	boolean done = false;
	
	public AbstractWeightCalculator(boolean save, LinkedBlockingQueue<Long> queue) {
		this.save = save;
		this.queue = queue;
		this.lastTime = System.currentTimeMillis();
	}
	
	public Long getWeight(T t, AbstractComputeMinCostTask<T> minCostTask) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight == null) {


			if(done) {
				
				weight =  calculateWeight(t, minCostTask);
				int count;
				if (save ) {
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
					weight =  calculateWeight(t, minCostTask);
					int count;
					if (save ) {
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

			int count;
			if (save) {
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
	
	
}