package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

public abstract class AbstractWeightCalculatorConsumer<T> extends AbstractWeightCalculator<T> implements Cloneable {
	
	private boolean save;
	private LinkedBlockingQueue<Long> queue;
	private boolean threadingOff = false;
	
	public AbstractWeightCalculatorConsumer(boolean save, LinkedBlockingQueue<Long> queue) {
		this.save = save;
		this.queue = queue;
		this.lastTime = System.currentTimeMillis();
	}
	
	public Long getWeight(T t) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight == null) {


			if(isThreadingOff()) {				
				weight =  calculateWeight(convertToSingletonArray(t))[0];
				if (save ) {
					weights.put(t, weight);
					weights.size();
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
				if(weight == TurnTaskToScores.THEEND) {//random number from CommandLine used as a "poison pill"
					setThreadingOff(true);
					weight =  calculateWeight(convertToSingletonArray(t))[0];
					int count;
					if (save ) {
						weights.put(t, weight);
						count = weights.size();
					} else {
						count = this.callcounter;
					}
					if (count % 100000 == 0) {
						Logging.log("Consumed "+ count +" weights; time (seconds): " + (System.currentTimeMillis() - lastTime)/1000);
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
				Logging.log("Calculated "+ count +" weights; time (seconds): " + (System.currentTimeMillis() - lastTime)/1000);
				lastTime = System.currentTimeMillis();
			}

		} else {
			//Logging.log("Found " + t );
		}
		//System.out.println(t.toString());
		return weight;

	}

	abstract T[] convertToSingletonArray(T t);

	public boolean isThreadingOff() {
		return threadingOff;
	}

	public void setThreadingOff(boolean done) {
		this.threadingOff = done;
		
	}
	
	
}