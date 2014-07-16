package phylonet.coalescent;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class WQDataCollection extends DataCollection<Tripartition> {


	List<STITreeCluster> treeAllClusters = new ArrayList<STITreeCluster>();

	Tripartition [] finalTripartitions = null;
	int [] finalCounts = null;	
	
	Integer [] geneTreesAsInts;

	float [][] distMatrix;
	private int n;
	private int algorithm;
	
	public WQDataCollection( WQClusterCollection clusters, int alg) {
		this.clusters = clusters;
		this.algorithm = alg;
	}
	
	
	void traverseTrees(List<Tree> trees, boolean fromGeneTrees, int n,
			Map<Tripartition, Integer> geneTreeTripartitonCount) {
		
		for (Tree tr : trees) {
	
			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
			STITreeCluster gtAll = new STITreeCluster();
			String[] gtLeaves = tr.getLeaves();
			for (int i = 0; i < gtLeaves.length; i++) {
				gtAll.addLeaf(GlobalMaps.taxonIdentifier.taxonId(gtLeaves[i]));
			}
			treeAllClusters.add(gtAll);
			for (TNode node : tr.postTraverse()) {				
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());
					
					STITreeCluster cluster = new STITreeCluster();
					Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
					cluster.addLeaf(taxonID);

					addBipartition(gtAll.getBitSet(), cluster,
						cluster.complementaryCluster());

					nodeToSTCluster.put(node, cluster);

				} else {
					int childCount = node.getChildCount();
					
					if (childCount >3 || (childCount == 3 && node != tr.getRoot()) ) {
					    if (fromGeneTrees) {
					        throw new RuntimeException(
								"not a bifurcating tree: " + tr + "\n"
										+ node);
					    } 
					}
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster((BitSet) bs.clone());

					//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
					nodeToSTCluster.put(node, cluster);
					
					int size = cluster.getClusterSize();

					STITreeCluster remaining = cluster.complementaryCluster();
					remaining.getBitSet().and(gtAll.getBitSet());
					
					if (addBipartition(gtAll.getBitSet(), cluster, remaining) && !fromGeneTrees) {
					    System.err.println("Extra bipartition added: " + cluster +" | "+remaining);
					}

					if (fromGeneTrees) {
					    if (childCount == 2 ) {
					        if (size != n) {
					            tryAddingTripartition( childbslist[0],  childbslist[1], 
					                    remaining, node, geneTreeTripartitonCount);
					        }
					    } else if (childCount == 3) {
					        tryAddingTripartition(childbslist[0], childbslist[1], childbslist[2] , 
					                node, geneTreeTripartitonCount);

					    } else {
					        throw new RuntimeException("hmmm?");
					        /*
					         * if (childCount == 2) { STITreeCluster l_cluster =
					         * childbslist[0];
					         * 
					         * STITreeCluster r_cluster = childbslist[1];
					         * 
					         * STITreeCluster allMinuslAndr_cluster =
					         * treeComplementary(null this should be
					         * gtCluster?,leaves);
					         * 
					         * STITreeCluster lAndr_cluster = cluster;
					         * 
					         * if (allMinuslAndr_cluster.getClusterSize() != 0) { //
					         * add Vertex STBs tryAddingSTB(l_cluster, r_cluster,
					         * cluster, node, true); tryAddingSTB( r_cluster,
					         * allMinuslAndr_cluster, null, node, true);
					         * tryAddingSTB(l_cluster, allMinuslAndr_cluster, null,
					         * node, true);
					         * 
					         * // Add the Edge STB tryAddingSTB(lAndr_cluster,
					         * allMinuslAndr_cluster, null, node, true); }
					         * 
					         * } else if (childCount == 3 && node.isRoot()) {
					         * STITreeCluster l_cluster = childbslist[0];
					         * 
					         * STITreeCluster m_cluster = childbslist[1];
					         * 
					         * STITreeCluster r_cluster = childbslist[2];
					         * 
					         * tryAddingSTB(l_cluster, r_cluster, null, node, true);
					         * tryAddingSTB(r_cluster, m_cluster, null, node, true);
					         * tryAddingSTB(l_cluster, m_cluster, null, node, true);
					         * } else { throw new
					         * RuntimeException("None bifurcating tree: "+ tr+ "\n"
					         * + node); }
					         */
					    }
					}
				}
			}

		}

	}

	private boolean addBipartition(BitSet gtAllBS,
			STITreeCluster c1, STITreeCluster c2) {
	    
	    boolean added = false;
	    
		STITreeCluster c1c = new STITreeCluster (c1);
		STITreeCluster c2c = new STITreeCluster (c2);
		if (distMatrix != null) {
			BitSet b1c = c1c.getBitSet();
			BitSet b2c = c2c.getBitSet();
			for (int i = gtAllBS.nextClearBit(0); i < n ; i = gtAllBS.nextClearBit(i+1)) {
				float dist1 = 0, dist2 = 0;
				for (int j = b1c.nextSetBit(0); j >= 0 ; j = b1c.nextSetBit(j+1)) {
					dist1 = Math.max(dist1 ,distMatrix[i][j]);
				}
				for (int j = b2c.nextSetBit(0); j >= 0 ; j = b2c.nextSetBit(j+1)) {
					dist2 = Math.max(dist2 ,distMatrix[i][j]);
				}
				if (dist1 > dist2) {
					b1c.set(i);
				} else {
					b2c.set(i);
				}
			}
		}
		
        // TODO: should this be treated differently?
        if (c1.getClusterSize() == 1) {
            //spm.addMissingIndividuals(c1.getBitSet());
            added |= addToClusters(c1, c1.getClusterSize());
        }
        
		SpeciesMapper spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
		spm.addMissingIndividuals(c1c.getBitSet());
		spm.addMissingIndividuals(c2c.getBitSet());
		
		int mainSize = c1c.getClusterSize();
		added |= addToClusters(c1c, mainSize);	
		int remSize  = c2c.getClusterSize();
		added |= addToClusters(c2c, remSize);
		
		c1c = c1c.complementaryCluster();
		c2c = c2c.complementaryCluster();

		mainSize = c1c.getClusterSize();
		added |= addToClusters(c1c, mainSize);  
		remSize  = c2c.getClusterSize();
		added |= addToClusters(c2c, remSize);

        return added;
        
	}

    
	public void computeTreePartitions(Inference<Tripartition> inference) {

		int k = inference.trees.size();
		n = GlobalMaps.taxonIdentifier.taxonCount();
		
		int haveMissing = 0;
		for (Tree tree : inference.trees) {
			if (tree.getLeafCount() != n) {
				haveMissing++;
			}
		}
		int [][] denom = new int [n][n];
		Deque<List<Integer>> stack = new ArrayDeque<List<Integer>>();
		if (haveMissing > k/3) {
			distMatrix = new float[n][n];
			for (Tree tree : inference.trees) {
				Integer treeall = tree.getLeafCount();
				for (TNode node : tree.postTraverse()) {
					if (node.isLeaf()) {
						ArrayList<Integer> tmp = new ArrayList<Integer>();
						tmp.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
						stack.push(tmp);
					} else if (node.isRoot() && node.getChildCount() == 3){
						List<Integer> left = stack.pop();
						List<Integer> middle = stack.pop();
						List<Integer> right = stack.pop();
						updateDistanceForTwoNodes(treeall, left, right);
						updateDistanceForTwoNodes(treeall, left, middle);
						updateDistanceForTwoNodes(treeall, middle, right);
						left.addAll(right);
						left.addAll(middle);
						stack.push(left);
					} else {
						List<Integer> left = stack.pop();
						List<Integer> right = stack.pop();
						updateDistanceForTwoNodes(treeall, left, right);
						left.addAll(right);
						right.clear();
						stack.push(left);
					}
					if (node.isRoot()) {
						List<Integer> all = stack.pop();
						int c = all.size() - 2;
						for (Integer l : all) {
							for (Integer r : all) {
								denom[l][r] += c*(c-1)/2;
								denom[r][l] = denom[l][r];
							}
						}
					}
				}
			}
			
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					if (denom[i][j] == 0)
						distMatrix[i][j] = 0;
					 else
					        distMatrix[i][j] = distMatrix[i][j] / denom[i][j];
					distMatrix[j][i] = distMatrix[i][j];
				}
			}
		}

		

		Map<Tripartition, Integer> geneTreeTripartitonCount = new HashMap<Tripartition, Integer>(k * n);
		//geneTreeInvalidSTBCont = new HashMap<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer>();
		// geneTreeRootSTBs = new HashMap<Tripartition, Integer>(k*n);
		// needed for fast version
		// clusterToSTBs = new HashMap<STITreeCluster, Set<Tripartition>>(k*n);

		STITreeCluster all = new STITreeCluster();
		all.getBitSet().set(0, n);
		addToClusters(all, GlobalMaps.taxonIdentifier.taxonCount());

		
		traverseTrees(inference.trees, true, n, geneTreeTripartitonCount);
		
		int s = 0;
		for (Integer c : geneTreeTripartitonCount.values()) {
			s += c;
		}
		System.err.println("Number of gene trees: " + k);
		System.err.println("Tripartitons in gene trees (count): "
				+ geneTreeTripartitonCount.size());
		System.err.println("Tripartitons in gene trees (sum): " + s);
		
		if (this.algorithm == -1) {
			this.algorithm = (n <= 32 || (geneTreeTripartitonCount.size() < k*6)) ? 2 : 1;
		}
		
		if (this.algorithm == 2) {
			System.err.println("Using tripartition-based weight calculation.");
		
			finalTripartitions = new Tripartition[geneTreeTripartitonCount.size()];
			finalCounts = new int[geneTreeTripartitonCount.size()];
			int i = 0;
			for (Entry<Tripartition, Integer> entry : geneTreeTripartitonCount.entrySet()){
				finalTripartitions[i] = entry.getKey();
				finalCounts[i] = entry.getValue();
				i++;
			}
		} else {
			System.err.println("Using tree-based weight calculation.");
			List<Integer> temp = new ArrayList<Integer>(); 
			
			for (Tree tr : inference.trees) {
				int internalNodes = 0;
				for (TNode node : tr.postTraverse()) {
					if (node.isLeaf()) {
						if (internalNodes != 0) {
							temp.add(-internalNodes);
							internalNodes = 0;
						}
						temp.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
					} else {
						internalNodes ++;
					}
					if (node.isRoot()) {
						temp.add(-internalNodes);
						temp.add(Integer.MIN_VALUE);
					}
				}
			}
			geneTreesAsInts = temp.toArray(new Integer[]{});
		}
		System.err.println("Number of Clusters: " + clusters.getClusterCount());

		inference.weightCalculator.initializeWeightContainer(
				geneTreeTripartitonCount.size() * 2);
		// System.err.println("sigma n is "+sigmaN);

	}

	private void updateDistanceForTwoNodes(Integer treeall, List<Integer> left,
			List<Integer> right) {
		int c = treeall - left.size() - right.size();
		c = c*(c-1)/2;
		for (Integer l : left) {
			for (Integer r : right) {
				distMatrix[l][r] += c;
				distMatrix[r][l] = distMatrix[l][r];
			}
		}
	}

	public void addExtraBipartitionsByInput(ClusterCollection extraClusters,
			List<Tree> trees, boolean extraTreeRooted) {

		traverseTrees(trees, false, GlobalMaps.taxonIdentifier.taxonCount(), null);
		int s = extraClusters.getClusterCount();
		/*
		 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
		 */
		System.err
				.println("Number of Clusters after additions from extra trees: "
						+ s);
	}

	private void tryAddingTripartition(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster remaining, TNode node,
			 Map<Tripartition, Integer> geneTreeTripartitonCount) {

			Tripartition trip = new Tripartition(l_cluster, r_cluster, remaining);
			geneTreeTripartitonCount.put(trip,
					geneTreeTripartitonCount.containsKey(trip) ? 
							geneTreeTripartitonCount.get(trip) + 1 : 1);
	}



}
