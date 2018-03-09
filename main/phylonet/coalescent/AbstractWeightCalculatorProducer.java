package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

public abstract class AbstractWeightCalculatorProducer<T> extends AbstractWeightCalculatorTask<T>{
	LinkedBlockingQueue<T> queue;
	public AbstractWeightCalculatorProducer(LinkedBlockingQueue<T> queue1) {
		super();
		this.queue = queue1;
	}

	@Override
	public Long getWeight(T t) {
		this.callcounter ++;
			
			try {
				queue.put(t);
			}
			catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
			
			if (this.callcounter % 100000 == 0) {
				System.err.println(callcounter + " weights produced : " + (double)(System.currentTimeMillis() - lastTime)/1000 + " seconds");
				lastTime = System.currentTimeMillis();
			}
			//System.out.println(t.toString());
		return 0L;
	}
	

}
