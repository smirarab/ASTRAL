package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.RecursiveTask;

import phylonet.coalescent.DuplicationWeightCounter.CalculateWeightTask;
import phylonet.coalescent.DuplicationWeightCounter.STBipartition;
import phylonet.coalescent.MGDInference_DP.TaxonNameMap;
import phylonet.coalescent.MGDInference_DP.Vertex;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

public class ComputeMinCostTask extends RecursiveTask<Integer> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 244989909835073096L;
	private MGDInference_DP inference;
	private Vertex v;

	@Override
	protected Integer compute() {
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}
	}

	public ComputeMinCostTask(MGDInference_DP inference, Vertex v) {
		this.inference = inference;
		this.v = v;
	}
	
	//final int maxEL = 10000000;
	
	private int computeMinCost() throws CannotResolveException {

		boolean rooted = inference.rooted;
		List<Tree> trees = inference.trees;
		DuplicationWeightCounter counter = inference.counter;
		TaxonNameMap taxonNameMap = inference.taxonNameMap;
		
		// -2 is used to indicate it cannot be resolved
		if (v._done == 2) {
			throw new CannotResolveException(v._cluster.toString());
		}
		// Already calculated. Don't re-calculate.
		if (v._done == 1) {
			return v._max_score;
		}
		//

		int clusterSize = v._cluster.getClusterSize();

		// SIA: base case for singelton clusters.
		if (clusterSize <= 1) {
			int _el_num = -1;
			// SIA: TODO: this is 0, right?
			if (inference.optimizeDuploss == 3) {
				if (taxonNameMap == null) {
					_el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, rooted);
				} else {
					_el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, taxonNameMap, rooted);
				}
			} else {
				_el_num = 0;
			}

			//v._min_cost = 0;
			v._max_score = - _el_num;
			v._min_lc = (v._min_rc = null);
			v._done = 1;
			return v._max_score;
		}
		Set<STBipartition> clusterBiPartitions = counter
				.getClusterBiPartitions(v._cluster);

		// STBipartition bestSTB = null;
		if (inference.fast) {
			
			fast_STB_based_inference(trees, counter,
					clusterBiPartitions);
			
		} else {
			List<Integer> El = new ArrayList<Integer>();
			for (int k = 0; k < trees.size(); k++) El.add(null);
			
			boolean tryAnotherTime = false;
			
			// First find what clusters are contained in this cluster
			ClusterCollection containedVertecies = inference.clusters.getContainedClusters(v._cluster);
			/*for (int i = 1; i <= (clusterSize / 2); i++) {
				List<Vertex> leftList = new ArrayList<Vertex>(
						inference.clusters.get(i));
				HashSet<Vertex> leftSet = new HashSet<Vertex>();
				containedVertecies.put(i, leftSet);
				for (Vertex smallV : leftList) {
					if (!v._cluster.containsCluster(smallV._cluster)) {
						continue;
					}
					leftSet.add(smallV);
					List<Vertex> rightList = new ArrayList<Vertex>(
							inference.clusters.get(clusterSize - i));
					HashSet<Vertex> rightSet = new HashSet<Vertex>();
					containedVertecies.put(clusterSize - i, rightSet);
					for (Vertex bigv : rightList) {
						if (!v._cluster.containsCluster(bigv._cluster)) {
							continue;
						}
						rightSet.add(bigv);
					}
				}
			}*/
				
			do {
				tryAnotherTime = false;
				
	
				for (STBipartition bi : containedVertecies.getClusterResolutions(v._cluster)) {
					try {
							Vertex smallV = containedVertecies.getVertexForCluster(bi.cluster1);
							Vertex bigv = containedVertecies.getVertexForCluster(bi.cluster2);
							ComputeMinCostTask smallWork = new ComputeMinCostTask(
									inference, smallV);
							ComputeMinCostTask bigWork = new ComputeMinCostTask(
									inference, bigv);
							CalculateWeightTask weigthWork = null;

							/*STBipartition bi = new STBipartition(
									smallV._cluster, bigv._cluster,
									v._cluster);*/

							Integer w = counter
									.getCalculatedBiPartitionDPWeight(bi);
							if (w == null) {
								weigthWork = counter.new CalculateWeightTask(
										bi);
								// MP_VERSION: smallWork.fork();
								w = weigthWork.compute();									
							}

							// MP_VERSION: smallWork.fork();
							Integer rscore = bigWork.compute();
							
							if (rscore == null) {
								// MP_VERSION: weigthWork.cancel(false);
								// MP_VERSION: smallWork.cancel(false);
								throw new CannotResolveException(
										bigv._cluster.toString());
							}

							Integer lscore;
							// MP_VERSION: lscore = smallWork.join();
							lscore = smallWork.compute();

							if (lscore == null) {
								// MP_VERSION: 	weigthWork.cancel(false);
								throw new CannotResolveException(
										smallV._cluster.toString());
							}
							// MP_VERSION: w = weigthWork.join();

							Integer e = 0;
							// If in duploss mode, need to get MDC cost as
							// well
							if (inference.optimizeDuploss == 3) {
								for (int k = 0; k < trees.size(); k++) {
									Tree tr = trees.get(k);
									STITreeCluster treeAll = inference.counter.treeAlls
											.get(k);
									if (smallV._cluster.isDisjoint(treeAll)
											|| bigv._cluster
													.isDisjoint(treeAll)) {
										// System.err
										// .println("skipping "+bi+" for "
										// +treeAll);
										continue;
									}
									if (El.get(k) == null) {
										if (taxonNameMap == null) {
											El.set(k,
													DeepCoalescencesCounter
															.getClusterCoalNum_rooted(
																	tr,
																	v._cluster));
										} else {
											El.set(k,
													DeepCoalescencesCounter
															.getClusterCoalNum_rooted(
																	tr,
																	v._cluster,
																	taxonNameMap));
										}
									} else {
										// System.err
										// .println("Used cached");
									}
									e += El.get(k);
									// System.err.println("E for " +
									// v._cluster + " is "+e + " and k is  "
									// + k);
								}
							} else {
								e = 0;
							}

							int c = inference.optimizeDuploss * w - e;

							if ((v._max_score != -1)
									&& (lscore + rscore + c < v._max_score)) {
								continue;
							}
							v._max_score = (lscore + rscore + c);
							v._min_lc = smallV;
							v._min_rc = bigv;
							v._c = c;

							//break; // Already found the only pair of
									// clusters whose union is v's cluster.
						} catch (CannotResolveException c) {
							// System.err.println("Warn: cannot resolve: " +
							// c.getMessage());
						}
					}
				if (v._min_lc == null || v._min_rc == null) {
					if (clusterSize <= 5) {
						counter.addAllPossibleSubClusters(v._cluster,
							containedVertecies);
						tryAnotherTime = true;
					} else if (clusterSize > 1) {
						/*if (clusterSize > 20) {
							System.err
								.println("Adding extra clusters (complementary of included clusters) for size "
										+ clusterSize + " : " + v._cluster);
						}*/
	
						//System.err.println(maxSubClusters);
						Iterator<Set<Vertex>> it = containedVertecies.getSubClusters();
						if (it.hasNext()) {
							Collection<Vertex> biggestSubClusters = new ArrayList<Vertex>(it.next());
							for (Vertex x : biggestSubClusters) {
								int i = x._cluster.getClusterSize();
								int complementarySize  = clusterSize - i;						
				
								tryAnotherTime |= containedVertecies.addCluster(counter.getCompleteryVertx(x, v._cluster),complementarySize);
							}
						}
						
					}
				}
			} while (tryAnotherTime); 
		}

		if (v._min_lc == null || v._min_rc == null) {
			if (MGDInference_DP._print) {
				System.err.println("WARN: No Resolution found for ( "
						+ v._cluster.getClusterSize() + " taxa ):\n"
						+ v._cluster);
			}
			v._done = 2;
			throw new CannotResolveException(v._cluster.toString());
		}
/*		if (clusterSize > 450){
			System.out.println(v+" \nis scored "+(v._max_score ) + " by \n"+v._min_lc + " \n"+v._min_rc);
		}
*/		/*
		 * if (clusterSize > 5){ counter.addGoodSTB(bestSTB, clusterSize); }
		 */
		v._done = 1;
		return v._max_score ;
	}

	private void fast_STB_based_inference( 
			List<Tree> trees, DuplicationWeightCounter counter,
			Set<STBipartition> clusterBiPartitions)
			throws CannotResolveException {
		TaxonNameMap taxonNameMap = inference.taxonNameMap;


		List<Integer> El = new ArrayList<Integer>();
		for (int k = 0; k < trees.size(); k++)
			El.add(null);

		if (clusterBiPartitions == null) {
			System.err.println("Warn: the following cluster ( "
					+ v._cluster.getClusterSize()
					+ " taxa ) has no STBs:\n" + v._cluster);
			v._max_score = -2;
			throw new CannotResolveException(v._cluster.toString());
		}
		for (STBipartition stb : clusterBiPartitions) {

			Vertex lv = inference.clusters.getVertexForCluster(stb.cluster1);
			Vertex rv = inference.clusters.getVertexForCluster(stb.cluster2);

			/*
			 * if (lv == null || rv == null) {
			 * //System.out.println("There is no STB for one half of : " +
			 * stb); continue; }
			 */

			try {

				// vertexStack.push(lv);
				ComputeMinCostTask worker1 = new ComputeMinCostTask(
						inference, lv);
				ComputeMinCostTask worker2 = new ComputeMinCostTask(
						inference, rv);
				CalculateWeightTask worker3 = null;

				worker1.fork();

				STBipartition bi = new STBipartition(lv._cluster,
						rv._cluster, v._cluster);

				Integer w = counter.getCalculatedBiPartitionDPWeight(bi);
				if (w == null) {
					worker3 = counter.new CalculateWeightTask(bi);
					worker3.fork();
				}
				Integer e = 0;
				// If in duploss mode, need to get MDC cost as well
				if (inference.optimizeDuploss == 3) {
					for (int k = 0; k < trees.size(); k++) {
						Tree tr = trees.get(k);
						STITreeCluster treeAll = inference.counter.treeAlls
								.get(k);
						if (rv._cluster.isDisjoint(treeAll)
								|| lv._cluster.isDisjoint(treeAll)) {
							// System.err
							// .println("skipping "+bi+" for " +treeAll);
							continue;
						}
						if (El.get(k) == null) {
							if (taxonNameMap == null) {
								El.set(k, DeepCoalescencesCounter
										.getClusterCoalNum_rooted(tr,
												v._cluster));
							} else {
								El.set(k, DeepCoalescencesCounter
										.getClusterCoalNum_rooted(tr,
												v._cluster, taxonNameMap));
							}
						} else {
							// System.err
							// .println("Used cached");
						}
						e += El.get(k);
						// System.err.println("E for " + v._cluster +
						// " is "+e + " and k is  " + k);
					}
				} else {
					e = 0;
				}

				Integer rscore = worker2.compute();
				if (rscore == null) {
					throw new CannotResolveException(rv._cluster.toString());
				}

				Integer lscore = worker1.join();
				if (lscore == null) {
					throw new CannotResolveException(lv._cluster.toString());
				}

				if (w == null) {
					w = worker3.join();
				}
				// vertexStack.pop();
				// vertexStack.push(rv);
				// vertexStack.pop();

				int c = inference.optimizeDuploss * w - e;

				if ((v._max_score != -1)
						&& (lscore + rscore + c <= v._max_score)) {
					continue;
				}
				v._max_score = (lscore + rscore + c);
/*				v._min_cost = inference.sigmaNs
						- (c + lv._max_score + rv._max_score - 2 * maxEL);*/
				// stem.out.println(maxEL - (z*w + lv._max_score +
				// rv._max_score));
				v._min_lc = lv;
				v._min_rc = rv;
				v._c = c;
			} catch (CannotResolveException c) {
				System.err.println("Warn: cannot resolve: "
						+ c.getMessage());
			}

			// bestSTB = stb;
		}
	}

}
