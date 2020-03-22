package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractComputeMinCostTask<T> {

	protected AbstractInference<T> inference;
	protected Vertex v;
	protected SpeciesMapper spm;

	final byte getDoneState = 1;
	final byte getOtherDoneState = 3;

	public AbstractComputeMinCostTask(AbstractInference<T> inference, Vertex v) {
		this.inference = inference;
		this.v = v;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
	}

	protected Double compute() {
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}
	
	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
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
					T t = STB2T(bi);					
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
	
				AbstractComputeMinCostTask<T> bigWork = newMinCostTask(bigv);
				AbstractComputeMinCostTask<T> smallWork = newMinCostTask(smallV);
	
				Double rscore = bigWork.compute();
	
				if (rscore == null) {
					throw new CannotResolveException(bigv.getCluster()
							.toString());
				}
					
				Double lscore;
				lscore = smallWork.compute();
	
				if (lscore == null) {
					throw new CannotResolveException(smallV
							.getCluster().toString());
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
				// System.err.println("Warn: cannot resolve: " +
				// c.getMessage());
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

	public Long getWeight(T t) {
		return inference.weightCalculator.getWeight(t);
	}


	abstract Long defaultWeightForFullClusters();

	protected abstract AbstractComputeMinCostTask<T> newMinCostTask(Vertex v);

	protected abstract double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom);

	protected abstract long calculateClusterLevelCost();

	protected abstract long scoreBaseCase(boolean rooted, List<Tree> trees);

	protected abstract T STB2T(VertexPair stb);

	/***
	 * Used in the exact version
	 * @param cluster
	 * @param containedVertecies
	 */
	void addAllPossibleSubClusters(STITreeCluster cluster, IClusterCollection containedVertecies) {
		int size = cluster.getClusterSize();
		for (int i = cluster.getBitSet().nextSetBit(0); i >= 0; i = cluster
				.getBitSet().nextSetBit(i + 1)) {
			STITreeCluster c = new STITreeCluster(cluster);
			c.getBitSet().clear(i);
	
			Vertex nv = c.new Vertex();
			containedVertecies.addCluster(nv, size - 1);
	
			addAllPossibleSubClusters(c, containedVertecies);
		}
	}
	



}
