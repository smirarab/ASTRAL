package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQComputeMinCostTaskConsumer extends AbstractComputeMinCostTask<Tripartition>{

	WQDataCollection wqDataCollection;
	
	public WQComputeMinCostTaskConsumer(AbstractInference<Tripartition> inference, Vertex v) {
		super(inference, v);
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
	protected AbstractComputeMinCostTask newMinCostTask( Vertex v) {
		return new WQComputeMinCostTaskConsumer(inference, v);
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
	
	double computeMinCost() throws CannotResolveException {
		
		// -2 is used to indicate it cannot be resolved
		if (v._done == 2) {
			throw new CannotResolveException(v.getCluster().toString());
		}
		// Already calculated. Don't re-calculate.
		if ( v._done == this.getDoneState || v._done == 4) {
			return v._max_score;
		}
		//
	
		int clusterSize = v.getCluster().getClusterSize();
	
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
	
			v._max_score = scoreBaseCase(inference.isRooted(), inference.trees);
	
			v._min_lc = (v._min_rc = null);
			if(v._done == getOtherDoneState)
				v._done = 4;
			else
				v._done = getDoneState;
	
			return v._max_score;
		}
	
		LinkedBlockingQueue<VertexPair> clusterResolutions;
		try {
			clusterResolutions = (LinkedBlockingQueue<VertexPair>) this.inference.getQueueClusterResolutions().take();
		} catch (InterruptedException e1) {
			throw new RuntimeException(e1);
		}
	
		long clusterLevelCost = 0;

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
	
				Vertex smallV = bi.cluster1;
				Vertex bigv = bi.cluster2;
	
				AbstractComputeMinCostTask<Tripartition> bigWork = newMinCostTask(bigv);
				AbstractComputeMinCostTask<Tripartition> smallWork = newMinCostTask(smallV);
	
				Double rscore = bigWork.compute();
	
				if (rscore == null) {
					throw new CannotResolveException(bigv.getCluster().toString());
				}
					
				Double lscore = smallWork.compute();
	
				if (lscore == null) {
					throw new CannotResolveException(smallV.getCluster().toString());
				}
	
				Long weight = weights[j];
				j++;
	
				double c = adjustWeight(clusterLevelCost, smallV, bigv, weight);	// Not relevant to ASTRAL				
	
				if ((v._max_score != -1)
						&& (lscore + rscore + c < v._max_score)) {
					continue;
				}
				if (lscore + rscore + c == v._max_score && GlobalMaps.random.nextBoolean()) {
					continue;
				}
				v._max_score = (lscore + rscore + c);
				v._min_lc = smallV;
				v._min_rc = bigv;
				v._c = c;
	
			} catch (CannotResolveException c) {
				c.printStackTrace();
				throw new RuntimeException("cannot resolve");
			}
			catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
	
		/**
		 * Should never happen
		 */
		if (v._min_lc == null || v._min_rc == null) {
			System.err.println("WARN: No Resolution found for ( "
					+ v.getCluster().getClusterSize() + " taxa ):\n"
					+ v.getCluster());
			v._done = 2;
			throw new CannotResolveException(v.getCluster().toString());
		}
		/*
		 * if (clusterSize > 450){
		 * System.out.println(v+" \nis scored "+(v._max_score ) +
		 * " by \n"+v._min_lc + " \n"+v._min_rc); }
		 *//*
		 * if (clusterSize > 5){ counter.addGoodSTB(bestSTB, clusterSize); }
		 */
		//johng23
		//if it's the consumer thread
		
		if(v._done == getOtherDoneState)			
			v._done = 4;
		else
			v._done = getDoneState;
	
		return v._max_score;
	}

}
