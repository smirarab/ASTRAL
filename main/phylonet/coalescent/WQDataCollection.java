package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class WQDataCollection extends DataCollection<Tripartition> {

	String[] gtTaxa;
	String[] stTaxa;

	List<STITreeCluster> treeAllClusters = new ArrayList<STITreeCluster>();

	Tripartition [] finalTripartitions = null;
	int [] finalCounts = null;	
	
	Integer [] geneTreesAsInts;

	
	public WQDataCollection(String[] gtTaxa, String[] stTaxa, WQClusterCollection clusters) {
		this.gtTaxa = gtTaxa;
		this.stTaxa = stTaxa;
		this.clusters = clusters;
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
					String nodeName = getSpeciesName(node.getName());
					
					STITreeCluster cluster = new STITreeCluster();
					Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
					cluster.addLeaf(taxonID);

					addToClusters(cluster, 1);
					
					addToClusters(cluster.complementaryCluster(), n - 1);

					nodeToSTCluster.put(node, cluster);
					

					// nodeToGTCluster.put(node, gtCluster);

				} else {
					int childCount = node.getChildCount();
					
					if (childCount >3 || (childCount == 3 && node != tr.getRoot()) ) {
						throw new RuntimeException(
								"not a bifurcating tree: " + tr + "\n"
										+ node);
					}
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(stTaxa.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster((BitSet) bs.clone());

					//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));

					int size = cluster.getClusterSize();

					addToClusters(cluster, size);
					nodeToSTCluster.put(node, cluster);
					
					STITreeCluster remaining = cluster.complementaryCluster();
					remaining.getBitSet().and(gtAll.getBitSet());
					int remSize  = remaining.getClusterSize();
					
					if (remSize != 0) {
						addToClusters(remaining, remSize);
					}


					if (childCount == 2 ) {
						if (size != n) {
							tryAddingTripartition( childbslist[0],  childbslist[1], 
									remaining, node, fromGeneTrees, geneTreeTripartitonCount);
						}
					} else if (childCount == 3) {
						tryAddingTripartition(childbslist[0], childbslist[1], childbslist[2] , 
								node, fromGeneTrees, geneTreeTripartitonCount);

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
						 */}
				}
			}

		}

	}

	public void computeTreePartitions(Inference<Tripartition> inference) {

		int k = inference.trees.size();
		int n = stTaxa.length;

		Map<Tripartition, Integer> geneTreeTripartitonCount = new HashMap<Tripartition, Integer>(k * n);
		//geneTreeInvalidSTBCont = new HashMap<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer>();
		// geneTreeRootSTBs = new HashMap<Tripartition, Integer>(k*n);
		// needed for fast version
		// clusterToSTBs = new HashMap<STITreeCluster, Set<Tripartition>>(k*n);

		STITreeCluster all = new STITreeCluster();
		all.getBitSet().set(0, n);
		addToClusters(all, stTaxa.length);

		
		traverseTrees(inference.trees, true, n, geneTreeTripartitonCount);
		
		int s = 0;
		for (Integer c : geneTreeTripartitonCount.values()) {
			s += c;
		}
		System.err.println("Number of gene trees: " + k);
		System.err.println("Tripartitons in gene trees (count): "
				+ geneTreeTripartitonCount.size());
		System.err.println("Tripartitons in gene trees (sum): " + s);
		
		if (n <= 32 || (geneTreeTripartitonCount.size() < k*8)) {
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

	public void addExtraBipartitionsByInput(ClusterCollection extraClusters,
			List<Tree> trees, boolean extraTreeRooted) {

		traverseTrees(trees, false, stTaxa.length, null);
		int s = extraClusters.getClusterCount();
		/*
		 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
		 */
		System.err
				.println("Number of Clusters After additions from extra Trees: "
						+ s);
	}

	private void tryAddingTripartition(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster remaining, TNode node,
			boolean fromGeneTrees, Map<Tripartition, Integer> geneTreeTripartitonCount) {
		
		if (fromGeneTrees) {
			Tripartition trip = new Tripartition(l_cluster, r_cluster, remaining);
			geneTreeTripartitonCount.put(trip,
					geneTreeTripartitonCount.containsKey(trip) ? 
							geneTreeTripartitonCount.get(trip) + 1 : 1);
		}
	}

	// static public int cnt = 0;


}
