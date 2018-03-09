package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractComputeMinCostTaskProducer<T> extends  AbstractComputeMinCostTask{

	//IClusterCollection clusters;

	//IClusterCollection containedVertecies;
	
	public AbstractComputeMinCostTaskProducer(AbstractInferenceProducer<T> inference, Vertex v) {
		super(inference, v);
		//this.clusters = clusters;
	}

	byte getDoneState() {
		return 3;
	}
	
	byte getOtherDoneState() {
		return 1;
	}

}
