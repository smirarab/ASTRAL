package phylonet.coalescent;

public class WriteTaskToQueue implements Runnable {
	AbstractInferenceNoCalculations inf;
	TurnTaskToScores threadgpu;
	public WriteTaskToQueue(AbstractInferenceNoCalculations inf, TurnTaskToScores threadgpu) {
		this.inf = inf;
		this.threadgpu = threadgpu;
	}
	public void run() {
		// TODO Auto-generated method stub
		inf.inferSpeciesTree();
		try {
			inf.queue1.put(TurnTaskToScores.POISON_PILL);
		}
		catch (Exception e) {
		}
		threadgpu.done = true;
	}
	
}