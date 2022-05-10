package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;

public class WQComputeMinCostTask extends AbstractComputeMinCostTask<Tripartition>{

	WQDataCollection wqDataCollection;
	
	protected VertexASTRAL3 v;
	
	public WQComputeMinCostTask(AbstractInference<Tripartition> inference, Vertex v,
			IClusterCollection clusters) {
		super(inference, clusters);
		this.wqDataCollection = (WQDataCollection)inference.dataCollection;
		this.v = (VertexASTRAL3) v;
	}
	
	public WQComputeMinCostTask(AbstractInference<Tripartition> inference, 
			IClusterCollection clusters) {
		super(inference,  clusters);
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
		return new Tripartition(vp.cluster1.getCluster(), 
				vp.cluster2.getCluster(), 
				vp.both.getCluster().complementaryCluster());
	}

	@Override
	protected Long defaultWeightForFullClusters() {
		return 0l;
	}
	
	/**
	 * This is the dynamic programming
	 * @return
	 * @throws CannotResolveException
	 */
	protected long computeMinCost() throws CannotResolveException {
		// Already calculated. Don't re-calculate.
		if (v.isDone() == 3) {
			return v._max_score;
		}
		
		if (v.isDone() == 0){
			long greedyScore = greedy();
			Logging.log("Greedy score: " + (long) greedyScore / 4);
			estimateUpperBound(v);
			inference.estimationFactor = v.get_upper_bound() / (greedyScore+1.0);
			Logging.log("estimationFactor: " + inference.estimationFactor);
			long estimateScore = estimateMinCost();
			Logging.log("Sub-optimal score: " + (long) estimateScore / 4);
		}
		
		//
		if (computeUpperBound(v) <= target) {
			return computeUpperBound(v);
		}
		
		int clusterSize = v.getCluster().getClusterSize();
	
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			
			v._min_lc = (v._min_rc = null);
			v.setDone(3);
			
			return v._max_score;
		}
		
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
				addComplementaryClusters(clusterSize, this.v);
			}
			clusterResolutions = containedVertecies.getClusterResolutions();
		}
		
		if (v.getClusterResolutions() != null) clusterResolutions = v.getClusterResolutions();
		else {
			ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
			
			for (VertexPair bi : clusterResolutions){
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
				else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
				computeUpperBound(bi.cluster1);
				computeUpperBound(bi.cluster2);
				bi.upperbound = ((VertexASTRAL3) bi.cluster1).get_upper_bound() + ((VertexASTRAL3) bi.cluster2).get_upper_bound() + bi.weight;
				clusterResolutionArrayList.add(bi);
			}
			
			Collections.sort(clusterResolutionArrayList);
			
			clusterResolutions = clusterResolutionArrayList;
		}
		for (VertexPair bi : clusterResolutions) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;
			
			long lscore = computeUpperBound(smallV), rscore = computeUpperBound(bigv);
			AbstractComputeMinCostTask smallWork = newMinCostTask(
					smallV, containedVertecies, v._max_score - bi.weight - rscore);
			lscore = smallWork.computeMinCost();
			
			AbstractComputeMinCostTask bigWork = newMinCostTask(
					bigv, containedVertecies, v._max_score - bi.weight - lscore);
			rscore = bigWork.computeMinCost();
			
			if (lscore + rscore + bi.weight <= v._max_score) {
				continue;
			}
			v._max_score = (lscore + rscore + bi.weight);
			v._min_lc = smallV;
			v._min_rc = bigv;
			v.set_c(bi.weight);
		}
		v.setClusterResolutions(null);
		v.setDone(3);
		return v._max_score;
	}
	
	long estimateMinCost(){
		estimateUpperBound(v);
		// Already calculated. Don't re-calculate.
		if (v.isDone() == 3) {
			return v._max_score;
		}
		if (v.isDone() == 2) {
			return v.get_estimated();
		}
		//
		if (v.isDone() == 1 && v.get_upper_bound() <= target * inference.estimationFactor) {
			return (long) (v.get_upper_bound() / inference.estimationFactor);
		}
		
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			v.set_estimated(0);
			v._min_lc = (v._min_rc = null);
			v.setDone(3);
			
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
				addComplementaryClusters(clusterSize, this.v);
			}
			clusterResolutions = containedVertecies.getClusterResolutions();
		}
		
		if (v.getClusterResolutions() != null) clusterResolutions = v.getClusterResolutions();
		else {
			ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
			
			for (VertexPair bi : clusterResolutions){
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
				else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
				estimateUpperBound(bi.cluster1);
				estimateUpperBound(bi.cluster2);
				bi.upperbound = ((VertexASTRAL3) bi.cluster1).get_upper_bound() + ((VertexASTRAL3) bi.cluster2).get_upper_bound() + bi.weight;
				clusterResolutionArrayList.add(bi);
			}
			
			Collections.sort(clusterResolutionArrayList);
			clusterResolutions = clusterResolutionArrayList;
			v.setClusterResolutions(clusterResolutionArrayList);
		}
		
		for (VertexPair bi : clusterResolutions) {
			VertexASTRAL3 smallV = (VertexASTRAL3) bi.cluster1;
			VertexASTRAL3 bigv = (VertexASTRAL3) bi.cluster2;

			long lscore = estimateUpperBound(smallV), rscore = estimateUpperBound(bigv);
			WQComputeMinCostTask smallWork = (WQComputeMinCostTask) newMinCostTask(
					smallV, containedVertecies, v.get_estimated() - bi.weight - rscore);
			lscore = smallWork.estimateMinCost();
			
			WQComputeMinCostTask bigWork = (WQComputeMinCostTask) newMinCostTask(
					bigv, containedVertecies, v.get_estimated() - bi.weight - lscore);
			rscore = bigWork.estimateMinCost();
			
			canSaveWork = (canSaveWork && smallV.isDone() == 3 && bigv.isDone() == 3);
			if (lscore + rscore + bi.weight <= v.get_estimated()) {
				continue;
			}
			v.set_estimated((lscore + rscore + bi.weight));
			v._min_lc = smallV;
			v._min_rc = bigv;
			v.set_c(bi.weight);
		}
		v.setDone(2);
		if (canSaveWork) {
			v.setDone(3);
			v._max_score = v.get_estimated();
		}
		return v.get_estimated();
	}
	
	protected long greedy(){
		long result = Long.MIN_VALUE;
		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
			
			v._max_score = 0;
			v.set_estimated(0);
			v._min_lc = (v._min_rc = null);
			v.setDone(3);
			
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
				addComplementaryClusters(clusterSize, this.v);
			}
			clusterResolutions = containedVertecies.getClusterResolutions();
		}
		
		if (v.getClusterResolutions() != null) clusterResolutions = v.getClusterResolutions();
		else {
			ArrayList<VertexPair> clusterResolutionArrayList = new ArrayList<VertexPair>();
			
			for (VertexPair bi : clusterResolutions){
				if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) bi.weight = defaultWeightForFullClusters();
				else bi.weight = inference.weightCalculator.getWeight(STB2T(bi), this);
				estimateUpperBound(bi.cluster1);
				estimateUpperBound(bi.cluster2);
				bi.upperbound = ((VertexASTRAL3) bi.cluster1).get_upper_bound() + ((VertexASTRAL3) bi.cluster2).get_upper_bound() + bi.weight;
				clusterResolutionArrayList.add(bi);
			}
			
			Collections.sort(clusterResolutionArrayList);
			clusterResolutions = clusterResolutionArrayList;
			v.setClusterResolutions(clusterResolutionArrayList);
		}
		
		for (VertexPair bi : clusterResolutions) {
			Vertex smallV = bi.cluster1;
			Vertex bigv = bi.cluster2;

			long lscore = estimateUpperBound(smallV), rscore = estimateUpperBound(bigv);
			WQComputeMinCostTask smallWork = (WQComputeMinCostTask) newMinCostTask(
					smallV, containedVertecies, result - bi.weight - rscore);
			lscore = smallWork.greedy();
			
			WQComputeMinCostTask bigWork = (WQComputeMinCostTask) newMinCostTask(
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
	
	Long computeUpperBound(Vertex v1v){
		VertexASTRAL3 v1 = (VertexASTRAL3) v1v;
		if (v1.isDone() == 3) return v1._max_score;
		if (v1.isDone() == 2) {
			if (v1.get_upper_bound() < v1.get_estimated() * inference.estimationFactor) return v1.get_upper_bound();
			return (long) (v1.get_estimated() * inference.estimationFactor);
		}
		if (v1.isDone() == 1) return v1.get_upper_bound();
		STITreeCluster c = v1.getCluster();
		v1.set_upper_bound(inference.weightCalculator.getWeight(new Tripartition(c, c, c, false), this));
		v1.setDone((byte) 1);
		return v1.get_upper_bound();
	}

	Long estimateUpperBound(Vertex v1v){
		VertexASTRAL3 v1 = (VertexASTRAL3) v1v;
		if (v1.isDone() == 3) return v1._max_score;
		if (v1.isDone() == 2) return v1.get_estimated();
		if (v1.isDone() == 1) return v1.get_upper_bound();
		STITreeCluster c = v1.getCluster();
		v1.set_upper_bound(inference.weightCalculator.getWeight(new Tripartition(c, c, c, false), this));
		v1.setDone((byte) 1);
		return v1.get_upper_bound();
	}


}
