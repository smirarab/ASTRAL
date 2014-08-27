package phylonet.coalescent;

import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQComputeMinCostTask extends AbstractComputeMinCostTask<Tripartition>{

	WQDataCollection wqDataCollection;
	
	public WQComputeMinCostTask(AbstractInference<Tripartition> inference, Vertex v,
			IClusterCollection clusters) {
		super(inference, v, clusters);
		this.wqDataCollection = (WQDataCollection)inference.dataCollection;
	}
	
	protected double adjustWeight(int clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom) {	
		return Wdom;
	}
	
	@Override
	protected int scoreBaseCase(boolean rooted, List<Tree> trees) {	
		return 0;
	}

	@Override
	protected AbstractComputeMinCostTask<Tripartition> newMinCostTask(
			 Vertex v, IClusterCollection clusters) {
		return new WQComputeMinCostTask(inference, v, clusters);
	}
	
	@Override
	protected int calculateClusterLevelCost() {
		return 0;
	}

	@Override
	protected Tripartition STB2T(STBipartition stb) {
		return new Tripartition(stb.cluster1, stb.cluster2, stb.c.complementaryCluster());
	}

	@Override
	Long defaultWeightForFullClusters() {
		return 0l;
	}

}
