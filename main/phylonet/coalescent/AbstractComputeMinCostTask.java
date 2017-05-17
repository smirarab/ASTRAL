package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

/**
 * This class implements the dynamic programming
 * @author smirarab
 *
 * @param <T>
 */
public abstract class AbstractComputeMinCostTask<T> {

	static STITreeCluster all = null;
	
	AbstractInference<T> inference;
	Vertex v;
	IClusterCollection clusters;
	Double target = 0.0;

	IClusterCollection containedVertecies;
    private SpeciesMapper spm;

	protected Double compute() {
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}
	
	Double computeUpperBound(Vertex v1){
		if (v1._done == 1) return v1._max_score;
		if (v1._done == 3) return v1._upper_bound;
		STITreeCluster c = v1.getCluster();
		if (all == null) all = (new STITreeCluster()).complementaryCluster();
		v1._upper_bound = inference.weightCalculator.getWeight((T) new Tripartition(c, c, all, false), this) / 2.0
				- inference.weightCalculator.getWeight((T) new Tripartition(c, c, c, false), this) / 3.0;
		v1._done = 3;
		return v1._upper_bound;
	}

	public AbstractComputeMinCostTask(AbstractInference<T> inference, Vertex v, 
			IClusterCollection clusters) {
		this.inference = inference;
		this.v = v;
		this.clusters = clusters;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
	}

	private void addComplementaryClusters(int clusterSize) {
		Iterator<Set<Vertex>> it = containedVertecies.getSubClusters().iterator();
		while (it.hasNext()) {
			Collection<Vertex> subClusters = new ArrayList<Vertex>(it.next());
			int i = -1;
			for (Vertex x : subClusters) {
				i = i > 0 ? i : x.getCluster().getClusterSize();
				int complementarySize = clusterSize - i;
				containedVertecies.addCluster(
						getCompleteryVertx(x, v.getCluster()),
						complementarySize);
			}
			if (i < clusterSize * inference.getCD()) {
				return;
			}

		}
	}

	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
	
	private double computeMinCost() throws CannotResolveException {


		// -2 is used to indicate it cannot be resolved
		if (v._done == 2) {
			throw new CannotResolveException(v.getCluster().toString());
		}
		// Already calculated. Don't re-calculate.
		if (v._done == 1) {
			return v._max_score;
		}
		//
		if (v._done == 3 && v._upper_bound < target) {
			return v._upper_bound;
		}
		
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = scoreBaseCase(inference.isRooted(), inference.trees);
			
			v._min_lc = (v._min_rc = null);
			v._done = 1;
			
			return v._max_score;
		}

		boolean tryAnotherTime = false;

		containedVertecies = clusters.getContainedClusters(v);
		
		boolean canSaveWork = true;

		do {
			tryAnotherTime = false;

			Iterable<VertexPair> clusterResolutions;
			
			if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
				clusterResolutions = new ArrayList<VertexPair>();
				Vertex v1 = null;
				int smallestSize = 1;
				while (v1 == null) {
					Set<Vertex> cs = containedVertecies.getSubClusters(smallestSize);
					if (cs.size() != 0)
						v1 = cs.iterator().next();
					else 
						smallestSize++;
				}
				for (Vertex v2: containedVertecies.getSubClusters(GlobalMaps.taxonIdentifier.taxonCount()-smallestSize))
				{
					if (v1.getCluster().isDisjoint(v2.getCluster())) {
						VertexPair vp = new VertexPair(v1, v2, v);
						((ArrayList<VertexPair>) clusterResolutions).add(vp);
						break;
					}
				}
				
			} else {
				if (clusterSize >= GlobalMaps.taxonIdentifier.taxonCount() * inference.getCS()) { //obsolete
					addComplementaryClusters(clusterSize);
				}
				clusterResolutions = containedVertecies.getClusterResolutions();
			}
			
			long clusterLevelCost = 0;
			if (clusterResolutions.iterator().hasNext()) { // Not relevant to ASTRAL
				clusterLevelCost = calculateClusterLevelCost();
			}
			/*
			 * System.out.println("xL: "+this.v.getCluster() + " "+xl + " "+
			 * DeepCoalescencesCounter
			 * .getClusterCoalNum(this.inference.trees, this.v.getCluster(),
			 * taxonNameMap, true));
			 */
			
			for (VertexPair bi : clusterResolutions) {
				try {
					Vertex smallV = bi.cluster1;
					Vertex bigv = bi.cluster2;
					Long weight = null;
					if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
						weight = defaultWeightForFullClusters();
					}
					
					if (weight == null) {
						T t = STB2T(bi);					
						weight =  inference.weightCalculator.getWeight(t, this);
						//System.out.print(weight/Integer.MAX_VALUE);
					}					
					
					double c = adjustWeight(clusterLevelCost, smallV, bigv, weight);	// Not relevant to ASTRAL
					
					Double lscore = computeUpperBound(smallV), rscore = computeUpperBound(bigv);
					AbstractComputeMinCostTask<T> smallWork = newMinCostTask(
							smallV, containedVertecies, v._max_score - c - rscore);
					lscore = smallWork.compute();
					if (lscore == null) throw new CannotResolveException(smallV.getCluster().toString());
					
					AbstractComputeMinCostTask<T> bigWork = newMinCostTask(
							bigv, containedVertecies, v._max_score - c - lscore);
					rscore = bigWork.compute();
					if (rscore == null) throw new CannotResolveException(bigv.getCluster().toString());
					
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
					// System.err.println("Warn: cannot resolve: " +
					// c.getMessage());
				}
			}
			/**
			 * This should be irrelevant if everything in the set X has at least one resolution (which is our goal). 
			 */
			if (v._min_lc == null || v._min_rc == null) {
				if (clusterSize <= 5) {
					addAllPossibleSubClusters(v.getCluster(), containedVertecies);
					tryAnotherTime = true;
				} else if (clusterSize > 1) {
					// System.err.println(maxSubClusters);
					Iterator<Set<Vertex>> it = containedVertecies.getSubClusters().iterator();
					if (it.hasNext()) {
						Collection<Vertex> biggestSubClusters = new ArrayList<STITreeCluster.Vertex>(it.next());
						//if (it.hasNext()) {biggestSubClusters = new ArrayList<STITreeCluster.Vertex>(it.next());}
						int i = -1;
						for (Vertex x : biggestSubClusters) {
							i = i > 0 ? i : x.getCluster().getClusterSize();
							int complementarySize = clusterSize - i;
							if (complementarySize > 1) {
								tryAnotherTime |= containedVertecies
										.addCluster( getCompleteryVertx(x,v.getCluster()), complementarySize);
							}
						}
						/*
						 * if (tryAnotherTime && clusterSize > 10) {
						 * System.err .println("Adding up to " +
						 * biggestSubClusters.size()+" extra |"+i+
						 * "| clusters (complementary of included clusters) for size "
						 * + clusterSize + " : " + v.getCluster()+"\n" +
						 * containedVertecies.getClusterCount()); }
						 */
					}

				}
			}
			/**
			 * We never want this to happen
			 */
			if (tryAnotherTime) {
				System.err.println("... auto expanding: " + this.v.getCluster().getClusterSize()+" "+this.v.getCluster());
			}
		} while (tryAnotherTime);


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

		v._done = 1;
		return v._max_score;
	}

	abstract Long defaultWeightForFullClusters();

	protected abstract AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
	IClusterCollection clusters);

	protected AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
			IClusterCollection clusters, Double target){
		AbstractComputeMinCostTask<T> task = newMinCostTask(v, clusters);
		task.target = target;
		return task;
	}
	
	abstract protected double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom);

    abstract protected long calculateClusterLevelCost();
	
	abstract protected long scoreBaseCase(boolean rooted, List<Tree> trees);
	
	abstract protected T STB2T(VertexPair stb);

	/***
	 * Used in the exact version
	 * @param cluster
	 * @param containedVertecies
	 */
	void addAllPossibleSubClusters(STITreeCluster cluster,
			IClusterCollection containedVertecies) {
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

	public Vertex getCompleteryVertx(Vertex x, STITreeCluster refCluster) {
		STITreeCluster c = x.getCluster();

		STITreeCluster revcluster = new STITreeCluster(refCluster);
		revcluster.getBitSet().xor(c.getBitSet());
		Vertex reverse = revcluster.new Vertex();
		return reverse;
	}

}
