package phylonet.coalescent;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

public class WQInference extends AbstractInference<Tripartition> {
	
	int forceAlg = -1;
	long maxpossible;
	
	public WQInference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution, boolean duploss, int alg, int addExtra) {
		super(rooted, extrarooted, trees, extraTrees, exactSolution, addExtra);
		this.forceAlg = alg;
	}

	public void scoreGeneTree(Tree st) {
		long startTime = System.currentTimeMillis();

		mapNames();

		IClusterCollection clusters = newClusterCollection();


		dataCollection = newCounter(clusters);
		weightCalculator = newWeightCalculator();

		dataCollection.computeTreePartitions(this);
		Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
		long sum = 0l;
		for (TNode node: st.postTraverse()) {
			if (node.isLeaf()) {
				String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

				STITreeCluster cluster = new STITreeCluster();
				Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
				cluster.addLeaf(taxonID);
				STITreeCluster remaining = cluster.complementaryCluster();

				stack.add(cluster);

			} else {
				ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				for (TNode child: node.getChildren()) {
					STITreeCluster pop = stack.pop();
					childbslist.add(pop);
					bs.or(pop.getBitSet());
				}
				
				STITreeCluster cluster = new STITreeCluster();
				cluster.setCluster((BitSet) bs.clone());

				//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
				stack.add(cluster);


				STITreeCluster remaining = cluster.complementaryCluster();
				if (remaining.getClusterSize() != 0) {
					childbslist.add(remaining);
				}
				for (int i = 0; i < childbslist.size(); i++) {
					for (int j = i+1; j < childbslist.size(); j++) {
						for (int k = j+1; k < childbslist.size(); k++) {
							sum += weightCalculator.getWeight(
									new Tripartition(childbslist.get(i),  childbslist.get(j), childbslist.get(k))
									, null);
						}
					}					       
				}
				if (childbslist.size() > 3) {
					for (STITreeCluster chid :childbslist) {
						System.err.print(chid.getClusterSize()+" ");
					}
					System.err.println(" (polytomy)");
				}
			}
		}
		System.out.println(sum/4l);
	}


	@Override
	Long getTotalCost(Vertex all) {
		System.err.println("Normalized score: " + all._max_score/4./this.maxpossible);
		return (long) (all._max_score/4l);
	}


	@Override
	AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(AbstractInference<Tripartition> dlInference,
			Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTask( (WQInference) dlInference, all,  clusters);
	}

	IClusterCollection newClusterCollection() {
		return new WQClusterCollection(GlobalMaps.taxonIdentifier.taxonCount(), new HashMap<STITreeCluster, STITreeCluster.Vertex>());
	}
	
	AbstractDataCollection<Tripartition> newCounter(IClusterCollection clusters) {
		return new WQDataCollection((WQClusterCollection)clusters, this.forceAlg, this);
	}



	@Override
	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		return new WQWeightCalculator(this);
	}

}
