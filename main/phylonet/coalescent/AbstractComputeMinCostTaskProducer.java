package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractComputeMinCostTaskProducer<T> extends  AbstractComputeMinCostTask<T>{


	AbstractInferenceProducer<T> inference;
	
	public AbstractComputeMinCostTaskProducer(AbstractInferenceProducer<T> inference, Vertex v) {
		super(inference, v);
		this.inference = inference;
	}

	byte getDoneState() {
		return 3;
	}
	
	byte getOtherDoneState() {
		return 1;
	}
	

	@Override
	public Long getWeight(T t) {
		inference.weightCount++;
			
		try {
			inference.getQueueReadyTripartitions().put(t);
		}
		catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}

		return 0L;
	}

}
