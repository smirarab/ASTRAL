package phylonet.coalescent;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class DLDataCollection extends AbstractDataCollection<STBipartition>{

	double sigmaNs;

	// private List<Set<STBipartition>> X;

	Map<STBipartition, Integer> geneTreeSTBCount;
	Map<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer> geneTreeInvalidSTBCont;

	boolean rooted;

	// Used for dup loss calculations
	List<STITreeCluster> treeAlls = new ArrayList<STITreeCluster>();

	HashMap<STBipartition, Set<STBipartition>> alreadyWeigthProcessed = new HashMap<STBipartition, Set<STBipartition>>();

	public DLDataCollection(boolean rooted, DLClusterCollection clusters) {
		this.rooted = rooted;
		this.clusters = clusters;
	}

	public void formSetX(AbstractInference<STBipartition> inference) {

		double unweigthedConstant = 0;
		double weightedConstant = 0;
		int k = inference.trees.size();
		int n = GlobalMaps.taxonIdentifier.taxonCount();
		boolean duploss = (((DLInference)inference).getOptimizeDuploss() == 3);		

		geneTreeSTBCount = new HashMap<STBipartition, Integer>(k * n);
		geneTreeInvalidSTBCont = new HashMap<AbstractMap.SimpleEntry<STITreeCluster, STITreeCluster>, Integer>();
		// geneTreeRootSTBs = new HashMap<STBipartition, Integer>(k*n);
		// needed for fast version
		// clusterToSTBs = new HashMap<STITreeCluster, Set<STBipartition>>(k*n);

		STITreeCluster all = new STITreeCluster(GlobalMaps.taxonIdentifier);
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++) {
			all.addLeaf(i);
		}
		addToClusters(all, GlobalMaps.taxonIdentifier.taxonCount());

		for (int t = 0; t < inference.trees.size(); t++) {
			Tree tr = inference.trees.get(t);

			STITreeCluster allInducedByGT = new STITreeCluster(GlobalMaps.taxonIdentifier);

			String[] gtLeaves = tr.getLeaves();
			for (int i = 0; i < gtLeaves.length; i++) {
				allInducedByGT.addLeaf(
						GlobalMaps.taxonIdentifier.taxonId(GlobalMaps.taxonNameMap.getTaxonName(gtLeaves[i])));
			}
			treeAlls.add(allInducedByGT);
			int allInducedByGTSize = allInducedByGT.getClusterSize();

			weightedConstant += duploss ? 2 * (allInducedByGTSize - 1) : 0;
					
			unweigthedConstant += (tr.getLeafCount() - 1);

			Map<TNode, STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
			// Map<TNode,STITreeCluster> nodeToGTCluster = new HashMap<TNode,
			// STITreeCluster>(n);

			for (TNode node : tr.postTraverse()) {				
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					String nodeName = GlobalMaps.taxonNameMap.getTaxonName(node.getName());
					
					STITreeCluster cluster = new STITreeCluster(GlobalMaps.taxonIdentifier);
					cluster.addLeaf(GlobalMaps.taxonIdentifier.taxonId(nodeName));

					addToClusters(cluster, 1);

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
					BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					// BitSet gbs = new BitSet(leaves.length);
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
						// gbs.or(nodeToGTCluster.get(child).getBitSet());
					}

					// STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
					// gtCluster.setCluster(gbs);

					STITreeCluster cluster = new STITreeCluster(GlobalMaps.taxonIdentifier);
					cluster.setCluster((BitSet) bs.clone());

					int size = cluster.getClusterSize();

					addToClusters(cluster, size);
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

		inference.weightCalculator.initializeWeightContainer(
				geneTreeSTBCount.size() * 2);
		//inference.weightCalculatorNoCalculations.initializeWeightContainer(geneTreeSTBCount.size() * 2);
		// System.err.println("sigma n is "+sigmaN);

		if (inference.getDLbdWeigth() == -1) {			
			inference.setDLbdWeigth((weightedConstant + 2*k + 0.0D) / 2*(k*n));
			System.out.println("Estimated bd weight = " + inference.getDLbdWeigth());
		}
			
		System.err.println("Sigma N: " + sigmaNs);
		
		sigmaNs = (unweigthedConstant + (1 - inference.getDLbdWeigth()) * weightedConstant);
	}


	public void addExtraBipartitionsByInput(
			List<Tree> trees, boolean extraTreeRooted) {

		int n = GlobalMaps.taxonIdentifier.taxonCount();

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
					String nodeName = GlobalMaps.taxonNameMap.getTaxonName(treeName);

					STITreeCluster tb = new STITreeCluster(GlobalMaps.taxonIdentifier);
					tb.addLeaf(GlobalMaps.taxonIdentifier.taxonId(nodeName));

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
					BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					int index = 0;
					for (TNode child: node.getChildren()) {
						childbslist[index++] = nodeToSTCluster.get(child);
						bs.or(nodeToSTCluster.get(child).getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster(GlobalMaps.taxonIdentifier);
					cluster.setCluster((BitSet) bs.clone());

					addToClusters(cluster, cluster.getClusterSize());
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
				((DLClusterCollection)clusters).addGeneTreeSTB(stb, size);
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
			STITreeCluster lmr = new STITreeCluster(GlobalMaps.taxonIdentifier);
			lmr.setCluster(l_Minus_r);

			BitSet r_Minus_l = (BitSet) and.clone();
			r_Minus_l.xor(r_cluster.getBitSet());
			STITreeCluster rml = new STITreeCluster(GlobalMaps.taxonIdentifier);

			rml.setCluster(r_Minus_l);

			if (!rml.getBitSet().isEmpty()) {
				addToClusters(rml, rml.getClusterSize());
				// addSTBToX( new STBipartition(l_cluster, rml, cluster),size);
			}
			if (!lmr.getBitSet().isEmpty()) {
				addToClusters(lmr, lmr.getClusterSize());
				// addSTBToX(new STBipartition(lmr, r_cluster, cluster), size);
			}
		}
	}


	long maxPossibleScore(Tripartition trip) {
		return 0;
	}

	@Override
	public void removeExtraBipartitionsByInput(List<Tree> extraTrees, boolean extraTreeRooted) {
		// TODO Auto-generated method stub
		
	}
}
