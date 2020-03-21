package phylonet.coalescent;

import java.util.List;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import java.util.concurrent.*; 
import java.util.*;

public abstract class AbstractComputeMinCostTask<T> {

	protected AbstractInference<T> inference;
	protected Vertex v;
	protected SpeciesMapper spm;

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
	
	byte getDoneState() {
		return 1;
	}
	
	byte getOtherDoneState() {
		return 3;
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
		if ( v._done == this.getDoneState() || v._done == 4) {
			return v._max_score;
		}
		//
	
		int clusterSize = v.getCluster().getClusterSize();
	
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
	
			v._max_score = scoreBaseCase(inference.isRooted(), inference.trees);
	
			v._min_lc = (v._min_rc = null);
			if(v._done == getOtherDoneState())
				v._done = 4;
			else
				v._done = getDoneState();
	
			return v._max_score;
		}
	
		Iterable<VertexPair> clusterResolutions = this.inference.getClusterResolutions(this.v);
	
		long clusterLevelCost = 0;
		/*
		 * System.out.println("xL: "+this.v.getCluster() + " "+xl + " "+
		 * DeepCoalescencesCounter
		 * .getClusterCoalNum(this.inference.trees, this.v.getCluster(),
		 * taxonNameMap, true));
		 */
		//VertexPair bi = null;
		//BlockingQueue<VertexPair> clusterResolutions2  = new LinkedBlockingQueue<VertexPair>();
		BlockingQueue<Long> weights  = new LinkedBlockingQueue<Long>();
		//while ( true ) {
		for (VertexPair bi : clusterResolutions) {
			try {
				//bi = ((BlockingQueue<VertexPair>)clusterResolutions).take(); 
				//clusterResolutions2.put(bi);
				//if (bi == AbstractComputeMinCostTaskProducer.POISON_PILL) break;
				Long weight = null;
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
					weight = defaultWeightForFullClusters();
				} else {
					T t = STB2T(bi);					
					weight =  inference.weightCalculator.getWeight(t);
				}					
				weights.put(weight);
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
		//while ( true ) {
		for (VertexPair bi : clusterResolutions) {
			try {
				//bi = ((BlockingQueue<VertexPair>)clusterResolutions2).take(); 
				//clusterResolutions2.put(bi);
				//if (bi == AbstractComputeMinCostTaskProducer.POISON_PILL) break;
	
				Vertex smallV = bi.cluster1;
				Vertex bigv = bi.cluster2;
	
				AbstractComputeMinCostTask<T> smallWork = newMinCostTask(smallV);
				AbstractComputeMinCostTask<T> bigWork = newMinCostTask(bigv);
	
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
	
				Long weight = weights.take();
	
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
		
		if(v._done == getOtherDoneState())			
			v._done = 4;
		else
			v._done = getDoneState();
	
		return v._max_score;
	}

	public Long getWeight(T t) {
		Long weight;
		weight =  inference.weightCalculator.getWeight(t);
		return weight;
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
