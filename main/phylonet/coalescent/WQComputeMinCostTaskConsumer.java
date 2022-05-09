package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeClusterMP;
import phylonet.tree.model.sti.STITreeClusterMP.VertexMP;

public class WQComputeMinCostTaskConsumer extends AbstractComputeMinCostTask<Tripartition>{

	WQDataCollectionMP wqDataCollection;
	//final byte getDoneState = 3;
	//final byte getOtherDoneState = 3;
	VertexMP v;
	
	public WQComputeMinCostTaskConsumer(AbstractInference<Tripartition> inference, VertexMP v) {
		super(inference, inference.dataCollection.clusters);
		this.wqDataCollection = (WQDataCollectionMP)inference.dataCollection;
		this.v = v;
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
	protected AbstractComputeMinCostTask<Tripartition> newMinCostTask( Vertex v,
			IClusterCollection clusters) {
		return new WQComputeMinCostTaskConsumer(inference, (VertexMP) v);
	}



	@Override
	protected long calculateClusterLevelCost() {
		return 0l;
	}



	@Override
	protected Tripartition STB2T(VertexPair vp) {
		return new Tripartition(
				(STITreeClusterMP)vp.cluster1.getCluster(), 
				(STITreeClusterMP)vp.cluster2.getCluster(), 
				(STITreeClusterMP)vp.both.getCluster().complementaryCluster());
	}

	@Override
	protected Long defaultWeightForFullClusters() {
		return 0l;
	}
	
	@Override
	protected Long compute() {
				
		// Already calculated. Don't re-calculate.
		if ( v.isConsDone() ) {
			return v._max_score;
		}
		//
	
		int clusterSize = v.getCluster().getClusterSize();
	
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = scoreBaseCase(inference.isRooted(), inference.trees);

			v.setConsDone(true);
			
			return v._max_score;
		}
	
		LinkedBlockingQueue<VertexPair> clusterResolutions;
		try {
			clusterResolutions = (LinkedBlockingQueue<VertexPair>) GlobalQueues.instance.getQueueClusterResolutions().take();
		} catch (InterruptedException e1) {
			throw new RuntimeException(e1);
		}
	

		Long [] weights  = new Long[clusterResolutions.size()];

		int j=0;
		for (VertexPair bi : clusterResolutions) {
			try {
				Long weight = null;
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
					weight = defaultWeightForFullClusters();
				} else {
					Tripartition t = STB2T(bi);					
					weight =  inference.weightCalculator.getWeight(t);
				}					
				weights[j] = (weight);
				j++;
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
		

                j = 0;
                for (VertexPair bi : clusterResolutions) {
                        try {

                                Long rscore = newMinCostTask((VertexMP) bi.cluster2, this.wqDataCollection.clusters).compute();
                                Long lscore = newMinCostTask((VertexMP) bi.cluster1, this.wqDataCollection.clusters).compute();

                                Long weight = weights[j];
                                j++;

                                long newScore = lscore + rscore + weight;
                                if ( ( newScore < v._max_score) ||
                                                ((newScore == v._max_score) && GlobalMaps.random.nextBoolean()) ) {
                                        continue;
                                }
                                v._max_score = newScore;
                                v._min_lc = (VertexMP) bi.cluster1;
                                v._min_rc =  (VertexMP) bi.cluster2;

                        }
                        catch (Exception e) {
                                throw new RuntimeException(e);
                        }
                }
	
		/**
		 * Should never happen
		 */
		if (v._min_lc == null || v._min_rc == null) {
			Logging.log("WARN: No Resolution found for ( "
					+ v.getCluster().getClusterSize() + " taxa ):\n"
					+ v.getCluster());
			throw new RuntimeException(new CannotResolveException(v.getCluster().toString()));
		}
		
		v.setConsDone((true));
	
		return v._max_score;
	}


}
