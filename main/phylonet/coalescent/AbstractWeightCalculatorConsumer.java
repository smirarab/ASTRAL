package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

public abstract class AbstractWeightCalculatorConsumer<T> extends AbstractWeightCalculatorMP<T> implements Cloneable {

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
		if (weight != null) {
			return weight;
		}

		if(isThreadingOff()) { // After the main DP, for computing final score. 		
			weight =  calculateWeight(convertToSingletonArray(t))[0];
			saveWeight(t, weight);
			return weight;
		}

		try {
			weight = queue.take();
		}
		catch(Exception e) {
			throw new RuntimeException(e);
		}
		if(weight == TurnTaskToScores.THEEND) {// After the main DP, for computing final score, we switch to local calculator
			setThreadingOff(true);
			weight =  calculateWeight(convertToSingletonArray(t))[0];
		}
		saveWeight(t, weight);

		return weight;

	}

	private void saveWeight(T t, Long weight) {
		int count;
		if (save ) {
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

	abstract T[] convertToSingletonArray(T t);

	public boolean isThreadingOff() {
		return threadingOff;
	}

	public void setThreadingOff(boolean done) {
		this.threadingOff = done;

	}


}