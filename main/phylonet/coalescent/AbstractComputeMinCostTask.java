package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
	static double estimationFactor = 0.8;
	
	AbstractInference<T> inference;
	Vertex v;
	IClusterCollection clusters;
	double target = 0.0;
	
	IClusterCollection containedVertecies;
    private SpeciesMapper spm;

	protected Double compute() {
		return computeMinCost();
	}
	
	Double computeUpperBound(Vertex v1){
		if (v1._done == 1) return v1._max_score;
		if (v1._done == 3 || v1._done == 4 || v1._done == 5) return v1._upper_bound;
		STITreeCluster c = v1.getCluster();
		v1._upper_bound = inference.weightCalculator.getWeight((T) new Tripartition(c, c, c, false), this);
		v1._done = 3;
		return v1._upper_bound;
	}

	Double estimateUpperBound(Vertex v1){
		if (v1._done == 5) return v1._max_score;
		if (v1._done == 3) return v1._upper_bound;
		STITreeCluster c = v1.getCluster();
		v1._upper_bound = inference.weightCalculator.getWeight((T) new Tripartition(c, c, c, false), this);
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
	
	private double computeMinCost() {
		// Already calculated. Don't re-calculate.
		if (v._done == 1) {
			return v._max_score;
		}
		//
		if (( v._done == 3 || v._done == 5 ) && v._upper_bound <= target) {
			return v._upper_bound;
		}
		
		if (v._done == 0){
			inference.bestSoFar = estimateMinCost();
			System.err.println("Sub-optimal score: " + (long) inference.bestSoFar / 4);
		}
		
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			
			v._min_lc = (v._min_rc = null);
			v._done = 1;
			
			return v._max_score;
		}

		containedVertecies = clusters.getContainedClusters(v);

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
		
		ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
		
		for (VertexPair bi : clusterResolutions){
			if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
			else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
			bi.upperbound = computeUpperBound(bi.cluster1) + computeUpperBound(bi.cluster2) + bi.weight;
			clusterResolutionArrayList.add(bi);
		}
		
		Collections.sort(clusterResolutionArrayList);
		
		for (VertexPair bi : clusterResolutionArrayList) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;
			
			double lscore = computeUpperBound(smallV), rscore = computeUpperBound(bigv);
			AbstractComputeMinCostTask<T> smallWork = newMinCostTask(
					smallV, containedVertecies, v._max_score - bi.weight - rscore);
			lscore = smallWork.computeMinCost();
			
			AbstractComputeMinCostTask<T> bigWork = newMinCostTask(
					bigv, containedVertecies, v._max_score - bi.weight - lscore);
			rscore = bigWork.computeMinCost();
			
			if (lscore + rscore + bi.weight < v._max_score) {
				continue;
			}
			v._max_score = (lscore + rscore + bi.weight);
			v._min_lc = smallV;
			v._min_rc = bigv;
			v._c = bi.weight;
		}

		v._done = 1;
		return v._max_score;
	}

	private double estimateMinCost(){
		// Already calculated. Don't re-calculate.
		if (v._done == 1 || v._done == 5) {
			return v._max_score;
		}
		//
		if (v._done == 3 && v._upper_bound * estimationFactor <= target) {
			return v._upper_bound * estimationFactor;
		}
		
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			v._min_lc = (v._min_rc = null);
			v._done = 1;
			
			return v._max_score;
		}

		containedVertecies = clusters.getContainedClusters(v);

		boolean canSaveWork = true;
		
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
		
		ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
		
		for (VertexPair bi : clusterResolutions){
			if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
			else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
			bi.upperbound = computeUpperBound(bi.cluster1) + computeUpperBound(bi.cluster2) + bi.weight;
			clusterResolutionArrayList.add(bi);
		}
		
		Collections.sort(clusterResolutionArrayList);
		
		for (VertexPair bi : clusterResolutions) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;

			double lscore = estimateUpperBound(smallV), rscore = estimateUpperBound(bigv);
			AbstractComputeMinCostTask<T> smallWork = newMinCostTask(
					smallV, containedVertecies, v._max_score - bi.weight - rscore);
			lscore = smallWork.estimateMinCost();
			
			AbstractComputeMinCostTask<T> bigWork = newMinCostTask(
					bigv, containedVertecies, v._max_score - bi.weight - lscore);
			rscore = bigWork.estimateMinCost();
			
			canSaveWork = (canSaveWork && smallV._done == 1 && bigv._done == 1);
			if (lscore + rscore + bi.weight < v._max_score) {
				continue;
			}
			v._max_score = (lscore + rscore + bi.weight);
			v._min_lc = smallV;
			v._min_rc = bigv;
			v._c = bi.weight;
		}
		
		if (v._max_score / estimationFactor < v._upper_bound) v._upper_bound = v._max_score / estimationFactor;
		if (canSaveWork) v._done = 1;
		else v._done = 5;
		return v._max_score;
	}
	abstract Long defaultWeightForFullClusters();

	protected abstract AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
	IClusterCollection clusters);

	protected AbstractComputeMinCostTask<T> newMinCostTask(Vertex v, 
			IClusterCollection clusters, double target){
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
