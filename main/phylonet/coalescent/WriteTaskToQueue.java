package phylonet.coalescent;

public class WriteTaskToQueue implements Runnable {
	WQInferenceProducer inf;
	TurnTaskToScores threadgpu;
	public WriteTaskToQueue(WQInferenceProducer inf, TurnTaskToScores threadgpu) {
		this.inf = inf;
		this.threadgpu = threadgpu;
	}
	public void run() {
		inf.inferSpeciesTree();
		try {
			inf.getQueueReadyTripartitions().put(TurnTaskToScores.POISON_PILL);
		}
		catch (Exception e) {
		}
		threadgpu.done = true;
	}
	
}