package phylonet.coalescent;

import java.util.HashSet;
import java.util.List;
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
	protected Integer compute(){		
		try {
			return computeMinCost();
		} catch (CannotResolveException e) {
			return null;
		}		
	}
	
	public ComputeMinCostTask (MGDInference_DP inference,Vertex v) {
		this.inference = inference;
		this.v = v;
	}
	
	private int computeMinCost() throws CannotResolveException {
		int maxEL = inference.maxEL;
		TaxonNameMap taxonNameMap = inference.taxonNameMap;
		boolean rooted = inference.rooted;
		List<Tree> trees = inference.trees;
		DuplicationWeightCounter counter = inference.counter;
		
		if (v._max_score == -2) {
			throw new CannotResolveException(v._cluster.toString());
		}
		// Already calculated. Don't re-calculate.
		if (v._max_score != -1) {
			return v._max_score - maxEL;
		}		
		// If in duploss mode, need to get MDC cost as well
		if (inference.optimizeDuploss == 3) {
			if (v._el_num == -1) {
				if (taxonNameMap == null) {
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, rooted);
				} else {
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, taxonNameMap, rooted);
				}
			}
		} else {
			v._el_num = 0;
		}
		
		int clusterSize = v._cluster.getClusterSize();
		
		// SIA: base case for singelton clusters.
		{
			if (clusterSize <= 1) {
				// SIA: TODO: this is 0, right?
				v._min_cost = 0;
				v._max_score = maxEL - v._el_num;
				v._min_lc = (v._min_rc = null);
				return v._max_score - maxEL;
			}
		}
		Set<STBipartition> clusterBiPartitions = counter
				.getClusterBiPartitions(v._cluster);

		// STBipartition bestSTB = null;
	 	if (inference.fast) {	
			if (clusterBiPartitions == null) {
				if (v._cluster.getClusterSize() <= 3) {
					for (int j=0; j< v._cluster.getClusterSize(); j++){
						STITreeCluster c1 = new STITreeCluster(v._cluster.getTaxa());
						c1.addLeaf(v._cluster.getClusterLeaves()[j]);
						STITreeCluster c2 = new STITreeCluster(v._cluster.getTaxa());
						for (int i = 0; i < v._cluster.getClusterSize(); i++) {
							if (i != j){
								c2.addLeaf(v._cluster.getClusterLeaves()[i]);
							}
						}
						if (inference.clusterToVertex.containsKey(c2)) {
							STBipartition stb = new STBipartition(c1, c2, v._cluster);
							if (clusterBiPartitions == null)
								clusterBiPartitions = new HashSet<STBipartition>(3);
							clusterBiPartitions.add(stb);
							System.err.println("Adding: " + stb);
						}
					}
				}
			}
	
			if (clusterBiPartitions == null) {
				System.err.println("Warn: the following cluster ( " + v._cluster.getClusterSize()+" taxa ) has no STBs:\n"
						+ v._cluster);
				v._max_score = -2;
				throw new CannotResolveException(v._cluster.toString());
			}
			for (STBipartition stb : clusterBiPartitions) {
	
				Vertex lv = inference.clusterToVertex.get(stb.cluster1);
				Vertex rv = inference.clusterToVertex.get(stb.cluster2);
								
				/*
				 * if (lv == null || rv == null) {
				 * //System.out.println("There is no STB for one half of : " + stb);
				 * continue; }
				 */
	
				try {
	
					//vertexStack.push(lv);
					ComputeMinCostTask worker1 = new ComputeMinCostTask(inference, lv);
					ComputeMinCostTask worker2 = new ComputeMinCostTask(inference, rv);
					CalculateWeightTask worker3 = null;
					
					worker1.fork();
					
					STBipartition bi = new STBipartition( lv._cluster,
							rv._cluster, v._cluster);
					
					Integer w = counter.getCalculatedBiPartitionDPWeight(bi);					
					if (w == null){						
						worker3 = counter.new CalculateWeightTask(bi);
						worker3.fork();
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
					//vertexStack.pop();
					//vertexStack.push(rv);					
					//vertexStack.pop();						
	
					int c = inference.optimizeDuploss * w - v._el_num;
						
					if ((v._max_score != -1)
							&& (lscore + rscore + c + maxEL <= v._max_score)) {
						continue;
					}
					v._max_score = (lscore + rscore + c) + maxEL;
					v._min_cost = inference.sigmaNs
							- (c + lv._max_score + rv._max_score - 2 * maxEL);
					// stem.out.println(maxEL - (z*w + lv._max_score +
					// rv._max_score));
					v._min_lc = lv;
					v._min_rc = rv;
					v._c = c;				
				} catch (CannotResolveException c) {
					System.err.println("Warn: cannot resolve: " + c.getMessage());
				}
	
				// bestSTB = stb;
			}
		} else {

			for (int i = 1; i <= (clusterSize / 2); i++) {
				Set<Vertex> leftList = inference.clusters.get(i);
				if (leftList != null) {
					for (Vertex smallV : leftList) {
						if (v._cluster.containsCluster(smallV._cluster)) {						
	
							Set<Vertex> rightList = inference.clusters.get(clusterSize - i);
							if (rightList != null) {							
								for (Vertex bigv : rightList) {
									if (! (smallV._cluster.isDisjoint(bigv._cluster) && v._cluster.containsCluster(bigv._cluster)) ) {
										continue;
									}
									try {
										ComputeMinCostTask smallWork = new ComputeMinCostTask(inference, smallV);
										ComputeMinCostTask bigWork = new ComputeMinCostTask(inference, bigv);
										CalculateWeightTask weigthWork = null;
																				
										STBipartition bi = new STBipartition( smallV._cluster,
												bigv._cluster, v._cluster);
										
										Integer w = counter.getCalculatedBiPartitionDPWeight(bi);					
										if (w == null){						
											weigthWork = counter.new CalculateWeightTask(bi);
											if (false) {
												smallWork.fork();
											} else {
												w = weigthWork.compute();
											}
										}
										
										if (false) {
											smallWork.fork();
										} 
										
										Integer rscore = bigWork.compute();
										if (rscore == null) {
											if (false) weigthWork.cancel(false);
											if (false) smallWork.cancel(false);
											throw new CannotResolveException(bigv._cluster.toString());
										}									
										
										Integer lscore;
										if (false) {										
											lscore = smallWork.join();
										} else {
											lscore = smallWork.compute();
										}
										
										if (lscore == null) {
											if (false) {
												weigthWork.cancel(false);
											}											
											throw new CannotResolveException(smallV._cluster.toString());
										}
										if (w == null) {
											w = weigthWork.join();
										}

	
										int c = inference.optimizeDuploss * w - v._el_num;
	
										if ((v._max_score != -1)
												&& (lscore + rscore + c + maxEL < v._max_score)) {
											continue;
										}
										v._max_score = (lscore + rscore + c) + maxEL;
										v._min_cost = inference.sigmaNs
												- (c + smallV._max_score + bigv._max_score - 2 * maxEL);
										// stem.out.println(maxEL - (z*w + lv._max_score +
										// rv._max_score));
										v._min_lc = smallV;
										v._min_rc = bigv;
										v._c = c;
	
										break;	// Already found the only pair of clusters whose union is v's cluster.
									} catch (CannotResolveException c) {										
										//System.err.println("Warn: cannot resolve: " + c.getMessage());											
									}
								}								
							}
						}
					}
				}
			}
	
		}
		
		if (v._min_lc == null || v._min_rc == null) {
			System.err.println("WARN: No Resolution found for ( " + v._cluster.getClusterSize()+" taxa ):\n"
					+ v._cluster);
			v._max_score = -2;
			throw new CannotResolveException(v._cluster.toString());
		}

		//System.out.println(v+" is scored "+(v._max_score - maxEL));
		/*
		 * if (clusterSize > 5){ counter.addGoodSTB(bestSTB, clusterSize); }
		 */
		return v._max_score - maxEL;
	}


	
}
