package phylonet.dl;

import java.util.AbstractMap.SimpleEntry;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;

import phylonet.coalescent.AbstractInference;
import phylonet.coalescent.AbstractWeightCalculator;
import phylonet.coalescent.GlobalMaps;
import phylonet.coalescent.WQInference;
import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;


class DLWeightCalculator extends AbstractWeightCalculator<STBipartition>{
	
	private DLDataCollection dataCollection;
	private DLInference inference;
	
	public DLWeightCalculator(AbstractInference<STBipartition> inference) {
		super(true);
		dataCollection = (DLDataCollection) inference.dataCollection;
		this.inference = (DLInference) inference;
	}
	
	Integer calculateHomomorphicCost(List<Integer> El, STITreeCluster cluster,
			Vertex smallV, Vertex bigv, List<Tree> trees) {
		Integer e = 0;
		for (int k = 0; k < trees.size(); k++) {
			Tree tr = trees.get(k);
			STITreeCluster treeAll = dataCollection.treeAlls.get(k);
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
		for (Entry<STBipartition, Integer> entry : dataCollection.geneTreeSTBCount.entrySet()) {
			STBipartition otherSTB = entry.getKey();
			/*if (cluster.containsCluster(otherSTB.c))
				continue;*/
			boolean c1 = cluster.containsCluster(otherSTB.cluster1);
			boolean c2 = cluster.containsCluster(otherSTB.cluster2);
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (Entry<SimpleEntry<STITreeCluster, STITreeCluster>, Integer> entry : dataCollection.geneTreeInvalidSTBCont.entrySet()) {
			SimpleEntry<STITreeCluster, STITreeCluster> otherSTB = entry.getKey();
			boolean c1 = cluster.containsCluster(otherSTB.getKey());
			boolean c2 = cluster.containsCluster(otherSTB.getValue());
			if ((c1 && !c2) || (c2 && !c1)) {
				weight += entry.getValue();
			}
		}
		for (STITreeCluster treeAll : dataCollection.treeAlls) {
			if (cluster.containsCluster(treeAll)) {
				weight++;
			}
		}
		int ret = weight;
		// XLweights.put(cluster, ret);
		return ret;
	}

	public void preCalculateWeights(List<Tree> trees, List<Tree> extraTrees) {

		if (inference.isRooted() && GlobalMaps.taxonNameMap == null && GlobalMaps.taxonIdentifier.taxonCount() > trees.size()) {
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
							Set<STBipartition> alreadyProcessedSTBs = dataCollection.alreadyWeigthProcessed
									.get(gtSTB);

							if (alreadyProcessedSTBs == null) {
								alreadyProcessedSTBs = new HashSet<STBipartition>(
										gtTrees.size() / 4);
								dataCollection.alreadyWeigthProcessed.put(gtSTB,
										alreadyProcessedSTBs);
							}

							if (alreadyProcessedSTBs.contains(stSTB)) {
								continue;
							}
							//TODO: this should happen in abstract class?
							weights.put(
									stSTB,
									(weights.containsKey(stSTB) ? weights
											.get(stSTB) : 0)
											+ dataCollection.geneTreeSTBCount.get(gtSTB));
							alreadyProcessedSTBs.add(stSTB);
						}
					}
				}
			}
		}
	}

	

	@Override
	protected
	Long calculateWeight(STBipartition stb) {
		// TODO: this is supposed to be
		// DLClusterCollection containedClusterCollection = (DLClusterCollection) task.containedVertecies;
		// We removed that because task is not passed around in new interface. Makes it much slower. 
		DLClusterCollection containedClusterCollection = null;
		// System.err.print("Calculating weight for: " + biggerSTB);
		Long weight = 0l;
		for (STBipartition smallerSTB : containedClusterCollection.getContainedGeneTreeSTBs()) {
			if (smallerSTB == stb || smallerSTB.isDominatedBy(stb)) {
				weight += dataCollection.geneTreeSTBCount.get(smallerSTB);
			}
		}
		// System.err.print(" ... " + weight);
		if (!dataCollection.rooted) {
			throw new RuntimeException("Unrooted not implemented.");
			/*
			 * for (STBipartition rootSTB : geneTreeRootSTBs.keySet()) { int c =
			 * geneTreeRootSTBs.get(rootSTB); STBipartition inducedSTB =
			 * biggerSTB.getInducedSTB(rootSTB.c); if
			 * (inducedSTB.equals(rootSTB)){ weight -= 2 * c;
			 * //System.err.print(" .. (" + rootSTB +" )" +c+" "+ weight); if
			 * (inducedSTB.cluster1.getClusterSize() != 1 &&
			 * inducedSTB.cluster2.getClusterSize() != 1) { weight -= 2 * c;
			 * //System.err.print(" . " + weight); } } }
			 */
		}
		// System.err.println("Weight of " + biggerSTB + " is " + weight);
		return weight;
	}

	@Override
	public void setupGeneTrees(WQInference wqInference) {
		// TODO Auto-generated method stub
		
	}
	
	public Long[] calculateWeight(STBipartition[] ts) {
		Long [] ret = new Long[ts.length];
		int i = 0;
		for (STBipartition t: ts)
			ret[i++] = this.calculateWeight(t);
		return ret;
	}

}
