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
		if (weight != null) {
			return weight;
		}

		if(isThreadingOff()) {				
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
		if(weight == TurnTaskToScores.THEEND) {//random number from CommandLine used as a "poison pill"
			setThreadingOff(true);
			weight =  calculateWeight(convertToSingletonArray(t))[0];
			saveWeight(t, weight);
			return weight;
		}


		saveWeight(t, weight);


		//System.out.println(t.toString());
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