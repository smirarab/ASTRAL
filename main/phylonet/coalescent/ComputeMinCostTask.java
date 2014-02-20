package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import phylonet.coalescent.DuplicationWeightCounter.CalculateWeightTask;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class ComputeMinCostTask {

	/**
	 * 
	 */
	private static final long serialVersionUID = 244989909835073096L;
	private MGDInference_DP inference;
	private Vertex v;
	private ClusterCollection clusters;

	protected Double compute() {
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}

	public ComputeMinCostTask(MGDInference_DP inference, Vertex v, ClusterCollection clusters) {
		this.inference = inference;
		this.v = v;
		this.clusters = clusters;
	}

	// final int maxEL = 10000000;
	ClusterCollection containedVertecies;

	private void add_complementary_clusters(int clusterSize) {
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
			if (i < clusterSize * inference.CD) {
				return;
			}

		}
	}

	private double computeMinCost() throws CannotResolveException {

		boolean rooted = inference.rooted;
		List<Tree> trees = inference.trees;
		DuplicationWeightCounter counter = inference.counter;

		// -2 is used to indicate it cannot be resolved
		if (v._done == 2) {
			throw new CannotResolveException(v.getCluster().toString());
		}
		// Already calculated. Don't re-calculate.
		if (v._done == 1) {
			return v._max_score;
		}
		//

		int clusterSize = v.getCluster().getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1) {
			int _el_num = -1;
			if (inference.optimizeDuploss == 3) {
				if (GlobalMaps.taxonNameMap == null) {
					_el_num = DeepCoalescencesCounter.getClusterCoalNum(trees,
							v.getCluster(), rooted);
					// System.out.println(v + " XL is " + _el_num);
				} else {
					_el_num = DeepCoalescencesCounter.getClusterCoalNumMap(trees,
							v.getCluster(), rooted);
				}
			} else {
				_el_num = 0;
			}

			// v._min_cost = 0;
			v._max_score = -_el_num;
			v._min_lc = (v._min_rc = null);
			v._done = 1;
			return v._max_score;
		}

		List<Integer> El = new ArrayList<Integer>();
		for (int k = 0; k < trees.size(); k++)
			El.add(null);

		boolean tryAnotherTime = false;

		containedVertecies = clusters.getContainedClusters(v.getCluster());

		do {
			tryAnotherTime = false;

			if (clusterSize >= inference.counter.stTaxa.length * inference.CS) {
				add_complementary_clusters(clusterSize);
			}
			Collection<STBipartition> clusterResolutions = containedVertecies
					.getClusterResolutions();
			int clusterDLCost = 0;
			if (inference.optimizeDuploss == 3 && !clusterResolutions.isEmpty()) {
				clusterDLCost = inference.counter.calculateDLstdClusterCost(
						this.v.getCluster(), inference.trees);
			}
			/*
			 * System.out.println("xL: "+this.v.getCluster() + " "+xl + " "+
			 * DeepCoalescencesCounter
			 * .getClusterCoalNum(this.inference.trees, this.v.getCluster(),
			 * taxonNameMap, true));
			 */
			for (STBipartition bi : clusterResolutions) {
				try {
					Vertex smallV = containedVertecies.getVertexForCluster(bi.cluster1);
					Vertex bigv = containedVertecies.getVertexForCluster(bi.cluster2);
					ComputeMinCostTask smallWork = new ComputeMinCostTask(
							inference, smallV, containedVertecies);
					ComputeMinCostTask bigWork = new ComputeMinCostTask(
							inference, bigv, containedVertecies);
					CalculateWeightTask weigthWork = null;

					// MP_VERSION: smallWork.fork();
					Double rscore = bigWork.compute();

					if (rscore == null) {
						// MP_VERSION: weigthWork.cancel(false);
						// MP_VERSION: smallWork.cancel(false);
						throw new CannotResolveException(bigv.getCluster()
								.toString());
					}

					Double lscore;
					// MP_VERSION: lscore = smallWork.join();
					lscore = smallWork.compute();

					if (lscore == null) {
						// MP_VERSION: weigthWork.cancel(false);
						throw new CannotResolveException(smallV
								.getCluster().toString());
					}
					// MP_VERSION: w = weigthWork.join();

					Integer Wdom = counter.getCalculatedBiPartitionDPWeight(bi);

					if (Wdom == null) {
						weigthWork = counter.new CalculateWeightTask(bi,
								containedVertecies);
						// MP_VERSION: smallWork.fork();
						Wdom = weigthWork.compute();
					}

					Integer e = 0;
					double c;
					if (inference.optimizeDuploss == 3) {
						
						int OverlappingGeneCount = 0;
						int someSideMissingXLCount = 0;
						int bothSidesPresentGeneCount = 0;
						for (int k = 0; k < inference.trees.size(); k++) {
							STITreeCluster treeAll = inference.counter.treeAlls.get(k);
							Tree tree = inference.trees.get(k);
							boolean pDisJoint = smallV.getCluster().isDisjoint(treeAll);
							boolean qDisJoint = bigv.getCluster().isDisjoint(treeAll);
							if (pDisJoint || qDisJoint) {
								someSideMissingXLCount +=  GlobalMaps.taxonNameMap == null ?
									DeepCoalescencesCounter.getClusterCoalNum_rooted(tree, this.v.getCluster()):
									DeepCoalescencesCounter.getClusterCoalNum_rootedMap(tree, this.v.getCluster());
							}
							if (!pDisJoint && !qDisJoint) {
								bothSidesPresentGeneCount += 1;
							}
							if (!pDisJoint || !qDisJoint) {
								OverlappingGeneCount += 1;
							}							
						}
						
						c = - ( clusterDLCost - 3 * Wdom - OverlappingGeneCount - someSideMissingXLCount 
								+ inference.DLbdWeigth * (someSideMissingXLCount + OverlappingGeneCount + bothSidesPresentGeneCount) );
					} else {
						c = Wdom;
					}					

					if ((v._max_score != -1)
							&& (lscore + rscore + c < v._max_score)) {
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
			if (v._min_lc == null || v._min_rc == null) {
				if (clusterSize <= 5) {
					addAllPossibleSubClusters(v.getCluster(),
							containedVertecies);
					tryAnotherTime = true;
				} else if (clusterSize > 1) {
					// System.err.println(maxSubClusters);
					Iterator<Set<Vertex>> it = containedVertecies
							.getSubClusters().iterator();
					if (it.hasNext()) {
						Collection<Vertex> biggestSubClusters = new ArrayList<Vertex>(
								it.next());
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
		} while (tryAnotherTime);


		if (v._min_lc == null || v._min_rc == null) {
			if (MGDInference_DP._print) {
				System.err.println("WARN: No Resolution found for ( "
						+ v.getCluster().getClusterSize() + " taxa ):\n"
						+ v.getCluster());
			}
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
		v._done = 1;
		return v._max_score;
	}

	int calculateDLbdAdjustment(Vertex smallV, Vertex bigv) {
		int e = 0;
		for (int k = 0; k < inference.trees.size(); k++) {
			STITreeCluster treeAll = inference.counter.treeAlls.get(k);
			boolean pDisJoint = smallV.getCluster().isDisjoint(treeAll);
			boolean qDisJoint = bigv.getCluster().isDisjoint(treeAll);
			if (!pDisJoint && !qDisJoint) {
				e += 1;
			}
			// System.err.println("E for " + v.getCluster() + " is "+e +
			// " and k is  " + k);
			// System.out.println(bigv + "|" + smallV+ " Total is " + e +
			// " extra is "+extraTerms + " for tree" +k);
		}
		return e;
	}

	void addAllPossibleSubClusters(STITreeCluster cluster,
			ClusterCollection containedVertecies) {
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
		// int size = reverse._cluster.getClusterSize();
		return reverse;
	}

	/*
	 * public boolean addAllPossibleSubClusters(STITreeCluster cluster) { int
	 * size = cluster.getClusterSize(); boolean ret = false; for (int i =
	 * cluster.getCluster().nextSetBit(0); i>=0 ;i =
	 * cluster.getCluster().nextSetBit(i+1)){ STITreeCluster c = new
	 * STITreeCluster(cluster); c.getCluster().clear(i); ret |= addToClusters(c,
	 * size-1); ret |= addAllPossibleSubClusters(c); } return ret; }
	 */

	/*
	 * public boolean addAllPossibleSubClusters(STITreeCluster cluster) { int
	 * size = cluster.getClusterSize(); boolean ret = false; for (int i =
	 * cluster.getCluster().nextSetBit(0); i>=0 ;i =
	 * cluster.getCluster().nextSetBit(i+1)){ STITreeCluster c = new
	 * STITreeCluster(cluster); c.getCluster().clear(i); ret |= addToClusters(c,
	 * size-1); ret |= addAllPossibleSubClusters(c); } return ret; }
	 */

}
