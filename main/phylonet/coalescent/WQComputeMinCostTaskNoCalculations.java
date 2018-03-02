package phylonet.coalescent;

import java.util.List;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQComputeMinCostTaskNoCalculations extends AbstractComputeMinCostTaskNoCalculations<Tripartition>{
	WQDataCollection wqDataCollection;
	
	public WQComputeMinCostTaskNoCalculations(AbstractInferenceNoCalculations<Tripartition> inference, Vertex v) {
		super(inference, v);
		this.wqDataCollection = (WQDataCollection)inference.dataCollection;
	}
	
	protected double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom) {	
		return Wdom;
	}
	
	@Override
	protected long scoreBaseCase(boolean rooted, List trees) {	
		return 0l;
	}

	@Override
	protected AbstractComputeMinCostTaskNoCalculations<Tripartition> newMinCostTask(Vertex v) {
		return new WQComputeMinCostTaskNoCalculations((AbstractInferenceNoCalculations<Tripartition>) inference, v);
	}
	
	@Override
	protected long calculateClusterLevelCost() {
		return 0l;
	}

	@Override
	protected Tripartition STB2T(VertexPair vp) {
		return new Tripartition(vp.cluster1.getCluster(), vp.cluster2.getCluster(), vp.both.getCluster().complementaryCluster());
	}

	@Override
	Long defaultWeightForFullClusters() {
		return 0l;
	}

}
