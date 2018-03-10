package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQComputeMinCostTaskProducer extends AbstractComputeMinCostTaskProducer<Tripartition>{
	WQDataCollection wqDataCollection;
	
	public WQComputeMinCostTaskProducer(AbstractInferenceProducer<Tripartition> inference, Vertex v) {
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
	protected AbstractComputeMinCostTaskProducer<Tripartition> newMinCostTask(Vertex v) {
		return new WQComputeMinCostTaskProducer((AbstractInferenceProducer<Tripartition>) inference, v);
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
