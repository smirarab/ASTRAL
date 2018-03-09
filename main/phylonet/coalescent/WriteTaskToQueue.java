package phylonet.coalescent;

public class WriteTaskToQueue implements Runnable {
	AbstractInferenceProducer inf;
	TurnTaskToScores threadgpu;
	public WriteTaskToQueue(AbstractInferenceProducer inf, TurnTaskToScores threadgpu) {
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