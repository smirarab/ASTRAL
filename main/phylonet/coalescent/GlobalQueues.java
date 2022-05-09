package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;

public class GlobalQueues {
	
	private LinkedBlockingQueue<Long> queueWeightResults;
	private LinkedBlockingQueue<Iterable<VertexPair>> queueClusterResolutions;

	
	static GlobalQueues instance = new GlobalQueues();
	
	private GlobalQueues() {}
	
	public LinkedBlockingQueue<Long> getQueueWeightResults() {
		return queueWeightResults;
	}


	public void setQueueWeightResults(LinkedBlockingQueue<Long> queueWeightResults) {
		this.queueWeightResults = queueWeightResults;
	}


	public LinkedBlockingQueue<Iterable<VertexPair>> getQueueClusterResolutions() {
		return queueClusterResolutions;
	}


	public void setQueueClusterResolutions(LinkedBlockingQueue<Iterable<VertexPair>> queueClusterResolutions) {
		this.queueClusterResolutions = queueClusterResolutions;
	}


}
