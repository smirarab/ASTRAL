package phylonet.coalescent;

import java.util.AbstractMap;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;

import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class DuplicationWeightCounter {

	HashMap<STBipartition, Integer> weights;
	String[] gtTaxa;
	String[] stTaxa;

	// private List<Set<STBipartition>> X;

	private Map<STBipartition, Integer> geneTreeSTBCount;
	private Map<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer> geneTreeInvalidSTBCont;

	private boolean rooted;

	private ClusterCollection clusters;

	// Used for dup loss calculations
	List<STITreeCluster> treeAlls = new ArrayList<STITreeCluster>();

	HashMap<STBipartition, Set<STBipartition>> alreadyWeigthProcessed = new HashMap<STBipartition, Set<STBipartition>>();

	public DuplicationWeightCounter(String[] gtTaxa, String[] stTaxa,
			boolean rooted,
			ClusterCollection clusters2) {
		this.gtTaxa = gtTaxa;
		this.stTaxa = stTaxa;
		this.rooted = rooted;
		this.clusters = clusters2;
	}

	private String getSpeciesName(String geneName) {
		String stName = geneName;
		if (GlobalMaps.taxonNameMap != null) {
			stName = GlobalMaps.taxonNameMap.getTaxonName(geneName);
		}
		return stName;
	}

	private boolean addToClusters(STITreeCluster c, int size,
			boolean geneTreeCluster) {
		Vertex nv = c.new Vertex();
		return clusters.addCluster(nv, size);
	}

	Integer calculateHomomorphicCost(List<Integer> El, STITreeCluster cluster,
			Vertex smallV, Vertex bigv, List<Tree> trees) {
		Integer e = 0;
		for (int k = 0; k < trees.size(); k++) {
			Tree tr = trees.get(k);
			STITreeCluster treeAll = treeAlls.get(k);
			if (smallV.getCluster().isDisjoint(treeAll)
					|| bigv.getCluster().isDisjoint(treeAll)) {
				continue;
			}
			if (El.get(k) == null) {
				if (GlobalMaps.taxonNameMap == null) {
					El.set(k, DeepCoalescencesCounter.getClusterCoalNum_rooted(
							tr, cluster));
				} else {
					El.set(k, DeepCoalescencesCounter.getClusterCoalNum_rootedMap(
							tr, cluster));
				}
			}
			e += El.get(k);
		}
		return e;
	}

	/*
	 * Calculates the cost of a cluster based on the standard definition
	 * */
	int calculateDLstdClusterCost(STITreeCluster cluster, List<Tree> trees) {
		/*
		 * if (XLweights.containsKey(cluster)) { return XLweights.get(cluster);
		 * }
		 */
		int weight = 0;
		for (Entry<STBipartition, Integer> entry : geneTreeSTBCount.entrySet()) {
			STBipartition otherSTB = entry.getKey();
			/*if (cluster.containsCluster(otherSTB.c))
				continue;*/
			boolean c1 = cluster.containsCluster(otherSTB.cluster1);
			boolean c2 = cluster.containsCluster(otherSTB.cluster2);
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (Entry<SimpleEntry<STITreeCluster, STITreeCluster>, Integer> entry : geneTreeInvalidSTBCont.entrySet()) {
			SimpleEntry<STITreeCluster, STITreeCluster> otherSTB = entry.getKey();
			boolean c1 = cluster.containsCluster(otherSTB.getKey());
			boolean c2 = cluster.containsCluster(otherSTB.getValue());
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (STITreeCluster treeAll : treeAlls) {
			if (cluster.containsCluster(treeAll)) {
				weight++;
			}
		}
		int ret = weight;
		// XLweights.put(cluster, ret);
		return ret;
	}

	double computeTreeSTBipartitions(MGDInference_DP inference) {

		double unweigthedConstant = 0;
		double weightedConstant = 0;
		int k = inference.trees.size();
		String[] leaves = stTaxa;
		int n = leaves.length;
		boolean duploss = (inference.optimizeDuploss == 3);		

		geneTreeSTBCount = new HashMap<STBipartition, Integer>(k * n);
		geneTreeInvalidSTBCont = new HashMap<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer>();
		// geneTreeRootSTBs = new HashMap<STBipartition, Integer>(k*n);
		// needed for fast version
		// clusterToSTBs = new HashMap<STITreeCluster, Set<STBipartition>>(k*n);

		STITreeCluster all = new STITreeCluster();
		for (int i = 0; i < stTaxa.length; i++) {
			all.addLeaf(i);
		}
		addToClusters(all, leaves.length, false);

		for (int t = 0; t < inference.trees.size(); t++) {
			Tree tr = inference.trees.get(t);

			STITreeCluster allInducedByGT = new STITreeCluster();

			String[] gtLeaves = tr.getLeaves();
			for (int i = 0; i < gtLeaves.length; i++) {
				String l = gtLeaves[i];
				allInducedByGT.addLeaf(getSpeciesName(l));
			}
			treeAlls.add(allInducedByGT);
			int allInducedByGTSize = allInducedByGT.getClusterSize();

			weightedConstant += duploss ? 
					2 * (allInducedByGTSize - 1) : 0;
					
			unweigthedConstant += (tr.getLeafCount() - 1);

			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
			// Map<TNode,STITreeCluster> nodeToGTCluster = new HashMap<TNode,
			// STITreeCluster>(n);

			for (TNode node : tr.postTraverse()) {				
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					String nodeName = getSpeciesName(node.getName());
					
					STITreeCluster cluster = new STITreeCluster();
					cluster.addLeaf(nodeName);

					addToClusters(cluster, 1, true);

					nodeToSTCluster.put(node, cluster);

					// nodeToGTCluster.put(node, gtCluster);

					if (!rooted) {
						throw new RuntimeException("Unrooted not implemented.");
						// tryAddingSTB(cluster, treeComplementary(null
						// /*gtCluster*/,leaves), null, node, true);
					}
				} else {
					int childCount = node.getChildCount();
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(leaves.length);
					// BitSet gbs = new BitSet(leaves.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
						// gbs.or(nodeToGTCluster.get(child).getBitSet());
					}

					// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
					// gtCluster.setCluster(gbs);

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster((BitSet) bs.clone());

					int size = cluster.getClusterSize();

					addToClusters(cluster, size, true);
					nodeToSTCluster.put(node, cluster);
					// nodeToGTCluster.put(node, gtCluster);

					if (rooted) {

						if (index > 2) {
							throw new RuntimeException(
									"None bifurcating tree: " + tr + "\n"
											+ node);
						}

						STITreeCluster l_cluster = childbslist[0];
						STITreeCluster r_cluster = childbslist[1];

						tryAddingSTB(l_cluster, r_cluster, cluster, node, true);

					} else {
						throw new RuntimeException("Unrooted not implemented.");
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

		int s = 0;
		for (Integer c : geneTreeSTBCount.values()) {
			s += c;
		}
		System.err.println("STBs in gene trees (count): "
				+ geneTreeSTBCount.size());
		System.err.println("STBs in gene trees (sum): " + s);

		s = clusters.getClusterCount();

		System.err.println("Number of Clusters: " + s);

		weights = new HashMap<STBipartition, Integer>(
				geneTreeSTBCount.size() * 2);
		// System.err.println("sigma n is "+sigmaN);

		if (inference.DLbdWeigth == -1) {			
			inference.DLbdWeigth = (weightedConstant + 2*k + 0.0D) / 2*(k*n);
			System.out.println("Estimated bd weight = " + inference.DLbdWeigth);
		}
			
		return (unweigthedConstant + (1 - inference.DLbdWeigth) * weightedConstant);
	}

	void addAllPossibleSubClusters(STITreeCluster cluster, int size) {
		BitSet bs = (BitSet) cluster.getBitSet().clone();
		bs.clear(0, size);
		while (true) {
			int tsb = bs.nextClearBit(0);
			if (tsb >= size) {
				break;
			}
			bs.set(tsb);
			bs.clear(0, tsb);
			STITreeCluster c = new STITreeCluster();
			c.setCluster((BitSet) bs.clone());
			addToClusters(c, c.getClusterSize(), false);
		}
		System.err
				.println("Number of Clusters After Adding All possible clusters: "
						+ clusters.getClusterCount());
	}

	void addExtraBipartitionsByInput(ClusterCollection extraClusters,
			List<Tree> trees, boolean extraTreeRooted) {

		String[] leaves = stTaxa;
		int n = leaves.length;

		// STITreeCluster all = extraClusters.getTopVertex().getCluster();

		for (Tree tr : trees) {

			/*
			 * String[] treeLeaves = tr.getLeaves(); STITreeCluster treeAll =
			 * new STITreeCluster(treeLeaves); for (int i = 0; i <
			 * treeLeaves.length; i++) { String l = treeLeaves[i];
			 * treeAll.addLeaf(l); }
			 */
			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(
					n);

			for (Iterator<TNode> nodeIt = tr.postTraverse().iterator(); nodeIt
					.hasNext();) {
				TNode node = nodeIt.next();
				if (node.isLeaf()) {
					String treeName = node.getName();
					String nodeName = getSpeciesName(treeName);

					STITreeCluster tb = new STITreeCluster();
					tb.addLeaf(nodeName);

					nodeToSTCluster.put(node, tb);

					if (!extraTreeRooted) {
						// TODO: fix the following (first null)
						throw new RuntimeException("Unrooted not implemented.");
						// tryAddingSTB(tb, tb.complementaryCluster(),
						// all,node,false);
					}
				} else {
					int childCount = node.getChildCount();
					STITreeCluster childbslist[] = new STITreeCluster[childCount];
					BitSet bs = new BitSet(leaves.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster((BitSet) bs.clone());

					addToClusters(cluster, cluster.getClusterSize(), false);
					nodeToSTCluster.put(node, cluster);

					if (extraTreeRooted) {

						if (index > 2) {
							throw new RuntimeException(
									"None bifurcating tree: " + tr + "\n"
											+ node);
						}

						STITreeCluster l_cluster = childbslist[0];

						STITreeCluster r_cluster = childbslist[1];

						tryAddingSTB(l_cluster, r_cluster, cluster, node, false);
					} else {
						throw new RuntimeException("Unrooted not implemented.");
						/*
						 * if (childCount == 2) { STITreeCluster l_cluster =
						 * childbslist[0];
						 * 
						 * STITreeCluster r_cluster = childbslist[1]; // Fix the
						 * following (first null) STITreeCluster
						 * allMinuslAndr_cluster = treeComplementary(null,
						 * leaves);
						 * 
						 * STITreeCluster lAndr_cluster = cluster;
						 * 
						 * // add Vertex STBs tryAddingSTB( l_cluster,
						 * r_cluster, cluster,node,false); if
						 * (allMinuslAndr_cluster.getClusterSize() != 0) {
						 * tryAddingSTB( r_cluster, allMinuslAndr_cluster,
						 * null,node,false); tryAddingSTB( l_cluster,
						 * allMinuslAndr_cluster, null,node,false);
						 * tryAddingSTB( lAndr_cluster, allMinuslAndr_cluster,
						 * all,node,false); }
						 * 
						 * } else if (childCount == 3 && node.isRoot()) {
						 * STITreeCluster l_cluster = childbslist[0];
						 * 
						 * STITreeCluster m_cluster = childbslist[1];
						 * 
						 * STITreeCluster r_cluster = childbslist[2];
						 * 
						 * tryAddingSTB( l_cluster, r_cluster, null,node,false);
						 * tryAddingSTB( r_cluster, m_cluster, null,node,false);
						 * tryAddingSTB( l_cluster, m_cluster, null,node,false);
						 * } else { throw new
						 * RuntimeException("None bifurcating tree: "+ tr+ "\n"
						 * + node); }
						 */
					}
				}
			}

		}
		int s = extraClusters.getClusterCount();
		/*
		 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
		 */
		System.err
				.println("Number of Clusters After additions from extra Trees: "
						+ s);
	}

	private void tryAddingSTB(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster cluster, TNode node,
			boolean fromGeneTrees) {
		// System.err.println("before adding: " + STBCountInGeneTrees);
		// System.err.println("Trying: " + l_cluster + "|" + r_cluster);
		int size = cluster.getClusterSize();
		if (l_cluster.isDisjoint(r_cluster)) {

			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);
			((STINode) node).setData(stb);
			if (fromGeneTrees) {
				clusters.addGeneTreeSTB(stb, size);
				// gtNodeToSTBs.put(node,stb);
				// addSTBToX(stb,size);
				// System.out.println(stb + " hashes to " + stb.hashCode());
				// if (! hash.containsKey(stb.hashCode()))
				// hash.put(stb.hashCode(), new HashSet<STBipartition>());
				// hash.get(stb.hashCode()).add(stb);
				geneTreeSTBCount.put(
						stb,
						geneTreeSTBCount.containsKey(stb) ? geneTreeSTBCount
								.get(stb) + 1 : 1);
			}

			/*
			 * if (size == allInducedByGTSize){ if (!
			 * geneTreeRootSTBs.containsKey(stb)) { geneTreeRootSTBs.put(stb,
			 * 1); } else { geneTreeRootSTBs.put(stb,
			 * geneTreeRootSTBs.get(stb)+1); } }
			 */
		} else {
			AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster> stb = 
					l_cluster.getBitSet().cardinality() > r_cluster.getBitSet().cardinality() ? 
					new AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>(l_cluster, r_cluster) : 
					new AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>(r_cluster, l_cluster);
			geneTreeInvalidSTBCont.put(stb,	geneTreeInvalidSTBCont.containsKey(stb) ? 
					geneTreeInvalidSTBCont.get(stb) + 1 : 
					1);
			// System.err.println("Adding only to extra");
			// This case could happen for multiple-copy
			BitSet and = (BitSet) l_cluster.getBitSet().clone();
			and.and(r_cluster.getBitSet());

			BitSet l_Minus_r = (BitSet) and.clone();
			l_Minus_r.xor(l_cluster.getBitSet());
			STITreeCluster lmr = new STITreeCluster();
			lmr.setCluster(l_Minus_r);

			BitSet r_Minus_l = (BitSet) and.clone();
			r_Minus_l.xor(r_cluster.getBitSet());
			STITreeCluster rml = new STITreeCluster();
			rml.setCluster(r_Minus_l);

			if (!rml.getBitSet().isEmpty()) {
				addToClusters(rml, rml.getClusterSize(), false);
				// addSTBToX( new STBipartition(l_cluster, rml, cluster),size);
			}
			if (!lmr.getBitSet().isEmpty()) {
				addToClusters(lmr, lmr.getClusterSize(), false);
				// addSTBToX(new STBipartition(lmr, r_cluster, cluster), size);
			}
		}
	}

	public Integer getCalculatedBiPartitionDPWeight(STBipartition bi) {
		if (!weights.containsKey(bi)) {
			// weights.put(bi,calculateMissingWeight(bi));
			return null;
		}
		return weights.get(bi);
	}

	// static public int cnt = 0;

	void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {

		if (rooted && GlobalMaps.taxonNameMap == null && stTaxa.length > trees.size()) {
			calculateWeightsByLCA(trees, trees);
			if (extraTrees != null) {
				calculateWeightsByLCA(extraTrees, trees);
			}
		}

	}

	void calculateWeightsByLCA(List<Tree> stTrees, List<Tree> gtTrees) {

		for (Tree stTree : stTrees) {
			SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(stTree);
			for (Tree gtTree : gtTrees) {
				Stack<TNode> stack = new Stack<TNode>();
				for (TNode gtNode : gtTree.postTraverse()) {
					if (gtNode.isLeaf()) {
						stack.push(stTree.getNode(gtNode.getName()));
					} else {
						TNode rightLCA = stack.pop();
						TNode leftLCA = stack.pop();
						// If gene trees are incomplete, we can have this case
						if (rightLCA == null || leftLCA == null) {
							stack.push(null);
							continue;
						}
						TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
						stack.push(lca);
						if (lca != leftLCA && lca != rightLCA) {
							// LCA in stTree dominates gtNode in gene tree
							// gtTree
							STBipartition stSTB = (STBipartition) ((STINode) lca)
									.getData();
							STBipartition gtSTB = (STBipartition) ((STINode) gtNode)
									.getData();
							Set<STBipartition> alreadyProcessedSTBs = alreadyWeigthProcessed
									.get(gtSTB);

							if (alreadyProcessedSTBs == null) {
								alreadyProcessedSTBs = new HashSet<STBipartition>(
										gtTrees.size() / 4);
								alreadyWeigthProcessed.put(gtSTB,
										alreadyProcessedSTBs);
							}

							if (alreadyProcessedSTBs.contains(stSTB)) {
								continue;
							}

							weights.put(
									stSTB,
									(weights.containsKey(stSTB) ? weights
											.get(stSTB) : 0)
											+ geneTreeSTBCount.get(gtSTB));
							alreadyProcessedSTBs.add(stSTB);
						}
					}
				}
			}
		}
	}

	class CalculateWeightTask {

		/**
		 * 
		 */
		private static final long serialVersionUID = -2614161117603289345L;
		private STBipartition stb;
		private ClusterCollection containedClusterCollection;

		public CalculateWeightTask(STBipartition stb,
				ClusterCollection collection) {
			this.stb = stb;
			this.containedClusterCollection = collection;
		}

		int calculateMissingWeight() {
			// System.err.print("Calculating weight for: " + biggerSTB);
			int weight = 0;
			for (STBipartition smallerSTB : containedClusterCollection
					.getContainedGeneTreeSTBs()) {
				if (smallerSTB == stb || smallerSTB.isDominatedBy(stb)) {
					weight += geneTreeSTBCount.get(smallerSTB);
				}
			}
			// System.err.print(" ... " + weight);
			if (!rooted) {
				throw new RuntimeException("Unrooted not implemented.");
				/*
				 * for (STBipartition rootSTB : geneTreeRootSTBs.keySet()) { int
				 * c = geneTreeRootSTBs.get(rootSTB); STBipartition inducedSTB =
				 * biggerSTB.getInducedSTB(rootSTB.c); if
				 * (inducedSTB.equals(rootSTB)){ weight -= 2 * c;
				 * //System.err.print(" .. (" + rootSTB +" )" +c+" "+ weight);
				 * if (inducedSTB.cluster1.getClusterSize() != 1 &&
				 * inducedSTB.cluster2.getClusterSize() != 1) { weight -= 2 * c;
				 * //System.err.print(" . " + weight); } } }
				 */
			}
			weights.put(stb, weight);
			// System.err.println("Weight of " + biggerSTB + " is " + weight);
			return weight;
		}

		protected Integer compute() {
			return calculateMissingWeight();
		}

	}

	/*
	 * public void addGoodSTB (STBipartition good, int size) {
	 * goodSTBs.get(size).add(good); }
	 */
	/*
	 * public Set<STBipartition> getClusterBiPartitions(STITreeCluster cluster)
	 * {
	 * 
	 * return clusterToSTBs.get(cluster); }
	 */

	/*
	 * public boolean addCompleteryVertx(Vertex x, STITreeCluster refCluster) {
	 * STITreeCluster c = x._cluster; Vertex reverse = new Vertex();
	 * reverse._cluster = new STITreeCluster(refCluster);
	 * reverse._cluster.getCluster().xor(c.getCluster()); int size =
	 * reverse._cluster.getClusterSize(); if
	 * (!clusters.get(size).contains(reverse)){ clusters.get(size).add(reverse);
	 * return true; //System.err.println("Clusters: "+clusters); } return false;
	 * }
	 */
	/*
	 * private void addSTBToX(STBipartition stb, int size) {
	 * //System.err.println("Adding to X: "+stb+" "+stb.c
	 * +" "+clusterToSTBs.containsKey(stb.c)); if
	 * (clusterToSTBs.containsKey(stb.c) &&
	 * clusterToSTBs.get(stb.c).contains(stb)){ return; } //int size =
	 * stb.c.getClusterSize(); // TODO: following line is algorithmically
	 * harmless, // but inefficient. is it necessary?
	 * //geneTreeSTBCount.put(stb, 0); addToClusters(stb.c, size, false); //
	 * Following needed for Fast //Set<STBipartition> stbs =
	 * clusterToSTBs.get(c); //stbs = (stbs== null)? new
	 * HashSet<STBipartition>() : stbs; //stbs.add(stb); //clusterToSTBs.put(c,
	 * stbs); //System.err.println("X updated: "+STBCountInGeneTrees);
	 * //System.err.println("X updated: "+clusterToSTBs); }
	 */

	/*
	 * private STITreeCluster treeComplementary(STITreeCluster gtCluster,
	 * String[] leaves){ //System.err.print("Tree complementary of "+gtCluster);
	 * STITreeCluster newGTCluster = gtCluster.complementaryCluster();
	 * //System.err.println(" is: "+newGTCluster.getCluster()); STITreeCluster
	 * newSTCluster = new STITreeCluster(); for (String s :
	 * newGTCluster.getClusterLeaves()) {
	 * newSTCluster.addLeaf(getSpeciesName(s)); }
	 * //System.err.println("Tree complementary of "
	 * +gtCluster+" is: "+newSTCluster); return newSTCluster; }
	 */

	/*
	 * private STITreeCluster treeComplementary(List<String> treeNames, Cluster
	 * c , TaxonNameMap taxonMap){ HashSet<String> set = new HashSet<String> ();
	 * set.add(cluster); return treeComplementary(treeNames, set, taxonMap); }
	 */

	/*
	 * void addExtraBipartitionsByHeuristics(ClusterCollection clusters2) {
	 * //goodSTBs = X; //if (true) return; int added = 0; for (int i=1;
	 * i<goodSTBs.size(); i++) { Set<STBipartition> curr_set = goodSTBs.get(i);
	 * for (STBipartition stb1:curr_set) { //if (Math.random() < 0.70) continue;
	 * for (int j=i; j<goodSTBs.size(); j++) { Set<STBipartition> other_set =
	 * goodSTBs.get(j); //if (Math.random() < 0.70) continue; for (STBipartition
	 * stb2:other_set) { //System.out.println(stb1 +" **AND** " + stb2); if
	 * (stb1.cluster1.getClusterSize() < 3 || stb1.cluster2.getClusterSize() < 3
	 * || stb2.cluster1.getClusterSize() < 3 || stb2.cluster2.getClusterSize() <
	 * 3) { if (tryToAdd(stb1,stb2,bipToAddToX) != null) added++; } }
	 * System.err.println(bipToAddToX.size() + " " + i); } } }
	 * 
	 * for (STBipartition stb: bipToAddToX) { //System.err.println( "Adding: " +
	 * stb); addSTBToX(clusters, stb); } System.out.println("\n\nAdded " +
	 * added+ " bipartitions:\n");
	 * 
	 * int s = 0; for (Integer c: clusters.keySet()){ s +=
	 * clusters.get(c).size(); }
	 * System.out.println("Number of Clusters After Addition: " +s);
	 * 
	 * }
	 * 
	 * 
	 * private STBipartition tryAddingExtraSTB_AndreRule(STBipartition stb1,
	 * STBipartition stb2, Set<STBipartition> bipToAddToX) { if
	 * (stb1.equals(stb2)) return null; if ( stb1.isDominatedBy(stb2) ||
	 * stb2.isDominatedBy(stb1) ) return null;
	 * 
	 * if ( stb1.c.isDisjoint(stb2.c) ) return null;
	 * 
	 * if ( stb1.cluster1.isDisjoint(stb2.cluster2) &&
	 * stb1.cluster2.isDisjoint(stb2.cluster1)) { STITreeCluster cl1 = new
	 * STITreeCluster(stb1.cluster1); cl1 = cl1.merge(stb2.cluster1);
	 * STITreeCluster cl2 = new STITreeCluster(stb1.cluster2); cl2 =
	 * cl2.merge(stb2.cluster2); STITreeCluster cl = new STITreeCluster(stb1.c);
	 * cl = cl.merge(stb2.c); STBipartition r = new STBipartition(cl1,cl2,cl);
	 * bipToAddToX.add(r); return r; } else if (
	 * stb1.cluster1.isDisjoint(stb2.cluster1) &&
	 * stb1.cluster2.isDisjoint(stb2.cluster2) ) { STITreeCluster cl1 = new
	 * STITreeCluster(stb1.cluster1); cl1 = cl1.merge(stb2.cluster2);
	 * STITreeCluster cl2 = new STITreeCluster(stb1.cluster2); cl2 =
	 * cl2.merge(stb2.cluster1); STITreeCluster cl = new STITreeCluster(stb1.c);
	 * cl = cl.merge(stb2.c); STBipartition r = new STBipartition(cl1,cl2,cl);
	 * bipToAddToX.add(r); return r; } return null; }
	 */
}
