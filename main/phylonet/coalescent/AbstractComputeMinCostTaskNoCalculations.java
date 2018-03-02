package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractComputeMinCostTaskNoCalculations<T> extends  AbstractComputeMinCostTask{

	//IClusterCollection clusters;

	//IClusterCollection containedVertecies;
	
	public AbstractComputeMinCostTaskNoCalculations(AbstractInferenceNoCalculations<T> inference, Vertex v) {
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
