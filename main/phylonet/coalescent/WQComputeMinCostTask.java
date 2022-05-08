package phylonet.coalescent;

import java.util.List;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;

public class WQComputeMinCostTask extends AbstractComputeMinCostTask<Tripartition>{

	WQDataCollection wqDataCollection;
	
	public WQComputeMinCostTask(AbstractInference<Tripartition> inference, Vertex v,
			IClusterCollection clusters) {
		super(inference, (VertexASTRAL3) v, clusters);
		this.wqDataCollection = (WQDataCollection)inference.dataCollection;
	}
	
	protected double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom) {	
		return Wdom;
	}
	
	@Override
	protected long scoreBaseCase(boolean rooted, List<Tree> trees) {	
		return 0l;
	}

	@Override
	protected AbstractComputeMinCostTask<Tripartition> newMinCostTask(
			 Vertex v, IClusterCollection clusters) {
		return new WQComputeMinCostTask(inference, v, clusters);
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
	protected Long defaultWeightForFullClusters() {
		return 0l;
	}

}
