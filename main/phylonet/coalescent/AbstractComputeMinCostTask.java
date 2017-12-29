package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.HashMap;
import java.util.Random;

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
	
	AbstractInference<T> inference;
	Vertex v;
	IClusterCollection clusters;
	double target = 0.0;
	boolean noPrecomputation = true;
	
	IClusterCollection containedVertecies;
    private SpeciesMapper spm;

	protected Double compute() {
		Random rnd = new Random(System.currentTimeMillis());
		int n = GlobalMaps.taxonIdentifier.taxonCount();
		long[] hash1 = new long[n], hash2 = new long[n];
		Vertex[][] arr = new Vertex[n + 1][];
		HashMap<Long, Vertex> map = new HashMap<Long, Vertex>();
		boolean succeed = false;
		while (!succeed) {
			map = new HashMap<Long, Vertex>();
			succeed = true;
			for (int i = 0; i < n; i++) {
				hash1[i] = rnd.nextLong();
				hash2[i] = rnd.nextLong();
			}
			for (int i = 1; i <= n; i++) {
				if (clusters.getSubClusters(i) == null) {
					arr[i] = new Vertex[0];
					continue;
				}
				arr[i] = new Vertex[clusters.getSubClusters(i).size()];
				int j = 0;
				for (Vertex v: clusters.getSubClusters(i)) {
					arr[i][j++] = v;
				}
			}
			for (int i = 1; i <= n; i++) {
				for (int j = 0; j < arr[i].length; j++) {
					long h1 = 0, h2 = 0;
					Vertex v = arr[i][j];
					BitSet b = v.getCluster().getBitSet();
					for (int k = 0; k < n; k++) {
						if (b.get(k)) {
							h1 += hash1[k];
							h2 ^= hash2[k];
						}
					}
					if (map.containsKey(h1)) {
						succeed = false;
						break;
					}
					map.put(h1, v);
					v.hash1 = h1;
					v.hash2 = h2;
					v.clusterSize = i;
				}
				if (!succeed) break;
			}
		}
		
		ArrayList<VertexPair> todolist = new ArrayList<VertexPair>();
		computeWeights(todolist, arr, map);
		Polytree.PTNative.compute(todolist);
		return computeMinCost();
	}
	
	Double computeUpperBound(Vertex v1){
		if (noPrecomputation) return 1e10;
		if (v1._done == 3) return v1._max_score;
		if (v1._done == 2) {
			if (v1._upper_bound < v1._estimated * inference.estimationFactor) return v1._upper_bound;
			return v1._estimated * inference.estimationFactor;
		}
		if (v1._done == 1) return v1._upper_bound;
		STITreeCluster c = v1.getCluster();
		v1._upper_bound = inference.weightCalculator.getWeight((T) new Tripartition(c, c, c, false), this);
		v1._done = 1;
		return v1._upper_bound;
	}

	Double estimateUpperBound(Vertex v1){
		if (noPrecomputation) return 1e10;
		if (v1._done == 3) return v1._max_score;
		if (v1._done == 2) return v1._estimated;
		if (v1._done == 1) return v1._upper_bound;
		STITreeCluster c = v1.getCluster();
		v1._upper_bound = inference.weightCalculator.getWeight((T) new Tripartition(c, c, c, false), this);
		v1._done = 1;
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
	
	private void computeWeights(ArrayList<VertexPair> todolist, Vertex[][] arr, HashMap<Long, Vertex> map) {
		if (v._done == 4) {
			return;
		}
		
		int clusterSize = v.getCluster().getClusterSize();
		v.clusterResolutions = new ArrayList<VertexPair>();
		
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			v._max_score = 0;
			v._min_lc = (v._min_rc = null);
			v._done = 4;
			
			return;
		}
		if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
			Vertex v1 = null, v2 = null;
			for (int i = 0; i < arr[1].length; i++) {
				v1 = arr[1][i];
				v2 = map.get(v.hash1 - v1.hash1);
				if (v2 != null) break;
			}
			v.clusterResolutions.add(new VertexPair(v1, v2, v));
		}
		else {
			for (int i = 1; i <= clusterSize / 2; i++) {
				Vertex[] ar = (arr[i].length < arr[clusterSize - i].length) ? arr[i] : arr[clusterSize - i];
				for (Vertex v1: ar) {
					Vertex v2 = map.get(v.hash1 - v1.hash1);
					if (v2 != null) {
						if (v.hash2 == (v1.hash2 ^ v2.hash2)
							&& v.clusterSize == v1.clusterSize + v2.clusterSize) {
							if (v1.clusterSize != v2.clusterSize || v1.hash1 < v2.hash1) {
								v.clusterResolutions.add(new VertexPair(v1, v2, v));
							}
						}
					}
				}
			}
		}
		for (VertexPair bi : v.clusterResolutions){
			if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
			else todolist.add(bi);
			
			newMinCostTask(bi.cluster1, containedVertecies).computeWeights(todolist, arr, map);
			newMinCostTask(bi.cluster2, containedVertecies).computeWeights(todolist, arr, map);
		}
		v._done = 4;
		return;
	}
	
	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
	
	private double computeMinCost() {
		// Already calculated. Don't re-calculate.
		if (v._done == 3) {
			return v._max_score;
		}
		/*
		if (!noPrecomputation && v._done == 0){
			double greedyScore = greedy();
			System.err.println("Greedy score: " + (long) greedyScore / 4);
			estimateUpperBound(v);
			inference.estimationFactor = v._upper_bound / greedyScore;
			System.err.println("estimationFactor: " + inference.estimationFactor);
			double estimateScore = estimateMinCost();
			System.err.println("Sub-optimal score: " + (long) estimateScore / 4);
		}
		
		//
		if (computeUpperBound(v) <= target) {
			return computeUpperBound(v);
		}
		*/
		int clusterSize = v.getCluster().getClusterSize();
		
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			
			v._min_lc = (v._min_rc = null);
			v._done = 3;
			
			return v._max_score;
		}
		/*
		Iterable<VertexPair> clusterResolutions;
		containedVertecies = clusters.getContainedClusters(v);
		
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
		
		if (v.clusterResolutions != null) clusterResolutions = v.clusterResolutions;
		else {
			ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
			
			for (VertexPair bi : clusterResolutions){
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
				else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
				computeUpperBound(bi.cluster1);
				computeUpperBound(bi.cluster2);
				bi.upperbound = bi.cluster1._upper_bound + bi.cluster2._upper_bound + bi.weight;
				clusterResolutionArrayList.add(bi);
			}
			
			Collections.sort(clusterResolutionArrayList);
			
			clusterResolutions = clusterResolutionArrayList;
		}*/
		for (VertexPair bi : v.clusterResolutions) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;
			
			double lscore = computeUpperBound(smallV), rscore = computeUpperBound(bigv);
			AbstractComputeMinCostTask<T> smallWork = newMinCostTask(
					smallV, containedVertecies); //, v._max_score - bi.weight - rscore);
			lscore = smallWork.computeMinCost();
			
			AbstractComputeMinCostTask<T> bigWork = newMinCostTask(
					bigv, containedVertecies); //, v._max_score - bi.weight - lscore);
			rscore = bigWork.computeMinCost();
			
			if (lscore + rscore + bi.weight <= v._max_score) {
				continue;
			}
			v._max_score = (lscore + rscore + bi.weight);
			v._min_lc = smallV;
			v._min_rc = bigv;
			v._c = bi.weight;
		}
		//v.clusterResolutions = null;
		v._done = 3;
		return v._max_score;
	}

	private double estimateMinCost(){
		estimateUpperBound(v);
		// Already calculated. Don't re-calculate.
		if (v._done == 3) {
			return v._max_score;
		}
		if (v._done == 2) {
			return v._estimated;
		}
		//
		if (v._done == 1 && v._upper_bound <= target * inference.estimationFactor) {
			return v._upper_bound / inference.estimationFactor;
		}
		
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			v._estimated = 0;
			v._min_lc = (v._min_rc = null);
			v._done = 3;
			
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
		
		if (v.clusterResolutions != null) clusterResolutions = v.clusterResolutions;
		else {
			ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
			
			for (VertexPair bi : clusterResolutions){
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
				else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
				estimateUpperBound(bi.cluster1);
				estimateUpperBound(bi.cluster2);
				bi.upperbound = bi.cluster1._upper_bound + bi.cluster2._upper_bound + bi.weight;
				clusterResolutionArrayList.add(bi);
			}
			
			Collections.sort(clusterResolutionArrayList);
			clusterResolutions = clusterResolutionArrayList;
			v.clusterResolutions = clusterResolutionArrayList;
		}
		
		for (VertexPair bi : clusterResolutions) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;

			double lscore = estimateUpperBound(smallV), rscore = estimateUpperBound(bigv);
			AbstractComputeMinCostTask<T> smallWork = newMinCostTask(
					smallV, containedVertecies, v._estimated - bi.weight - rscore);
			lscore = smallWork.estimateMinCost();
			
			AbstractComputeMinCostTask<T> bigWork = newMinCostTask(
					bigv, containedVertecies, v._estimated - bi.weight - lscore);
			rscore = bigWork.estimateMinCost();
			
			canSaveWork = (canSaveWork && smallV._done == 3 && bigv._done == 3);
			if (lscore + rscore + bi.weight <= v._estimated) {
				continue;
			}
			v._estimated = (lscore + rscore + bi.weight);
			v._min_lc = smallV;
			v._min_rc = bigv;
			v._c = bi.weight;
		}
		v._done = 2;
		if (canSaveWork) {
			v._done = 3;
			v._max_score = v._estimated;
		}
		return v._estimated;
	}
	
	private double greedy(){
		double result = -1e18;
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			v._estimated = 0;
			v._min_lc = (v._min_rc = null);
			v._done = 3;
			
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
		
		if (v.clusterResolutions != null) clusterResolutions = v.clusterResolutions;
		else {
			ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
			
			for (VertexPair bi : clusterResolutions){
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
				else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
				estimateUpperBound(bi.cluster1);
				estimateUpperBound(bi.cluster2);
				bi.upperbound = bi.cluster1._upper_bound + bi.cluster2._upper_bound + bi.weight;
				clusterResolutionArrayList.add(bi);
			}
			
			Collections.sort(clusterResolutionArrayList);
			clusterResolutions = clusterResolutionArrayList;
			v.clusterResolutions = clusterResolutionArrayList;
		}
		
		for (VertexPair bi : clusterResolutions) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;

			double lscore = estimateUpperBound(smallV), rscore = estimateUpperBound(bigv);
			AbstractComputeMinCostTask<T> smallWork = newMinCostTask(
					smallV, containedVertecies, result - bi.weight - rscore);
			lscore = smallWork.greedy();
			
			AbstractComputeMinCostTask<T> bigWork = newMinCostTask(
					bigv, containedVertecies, result - bi.weight - lscore);
			rscore = bigWork.greedy();
			
			if (lscore + rscore + bi.weight <= result) {
				continue;
			}
			result = (lscore + rscore + bi.weight);
			break;
		}
		return result;
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
