package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;
import java.util.Map.Entry;

import phylonet.coalescent.Utils.ClusterComparator;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

public class GreedyConsensus {

	protected GreedyConsensus() {};

	public static GreedyConsensus instance = new GreedyConsensus();
	/**
	 * Greedy consensus
	 * @param trees
	 * @param randomize
	 * @param taxonIdentifier
	 * @param keepclusters should we keep clusters as node objects
	 * @return
	 */
	public Tree greedyConsensus(Iterable<Tree> trees, boolean randomize,
			TaxonIdentifier taxonIdentifier, boolean keepclusters) {
		return greedyConsensus(trees,new double[]{0d}, randomize, 1, taxonIdentifier, keepclusters).iterator().next();
	}

	/**
	 * Greedy consensus
	 * @param trees
	 * @param randomize
	 * @param taxonIdentifier
	 * @param threshold
	 * @param geneTreeSkipProb 
	 * @return
	 */
	/*
    public static final Tree greedyConsensus(Iterable<Tree> trees, boolean randomize,
    		TaxonIdentifier taxonIdentifier, double threshold, double geneTreeKeepProb) {
    	return greedyConsensus(trees,new double[]{threshold}, randomize, 1, taxonIdentifier, geneTreeKeepProb).iterator().next();
    }
	 */
	/***
	 * Greedy consensus with a set of thresholds
	 * @param trees
	 * @param thresholds
	 * @param randomzie
	 * @param repeat
	 * @param taxonIdentifier
	 * @param keepclusters should we keep clusters as node objects
	 * @return
	 */
	public Collection<Tree> greedyConsensus(Iterable<Tree> trees, 
			double[] thresholds, boolean randomzie, int repeat, 
			TaxonIdentifier taxonIdentifier, boolean keepclusters) {
		Logging.logTimeMessage("Utils 219-222: ");

		List<Tree> outTrees = new ArrayList<Tree>();

		greedyConsensus(trees, thresholds, randomzie, repeat, taxonIdentifier, keepclusters, outTrees);

		return outTrees;
	}

	protected void greedyConsensus(Iterable<Tree> trees, double[] thresholds, boolean randomzie, int repeat,
			TaxonIdentifier taxonIdentifier, boolean keepclusters, List results) {
		HashMap<STITreeCluster, Integer> count = new HashMap<STITreeCluster, Integer>();
		int treecount = countClusteres(trees, taxonIdentifier, count);

		Logging.logTimeMessage("Utils 240-243: " );

		for (int gi = 0; gi < repeat; gi++) {
			TreeSet<Entry<STITreeCluster,Integer>> countSorted = new 
					TreeSet<Entry<STITreeCluster,Integer>>(new ClusterComparator(randomzie, taxonIdentifier.taxonCount()));

			countSorted.addAll(count.entrySet());

			int ti = thresholds.length - 1;
			double threshold = thresholds[ti];
			List<STITreeCluster> clusters = new ArrayList<STITreeCluster>(); 
			for (Entry<STITreeCluster, Integer> entry : countSorted) {
				if (threshold > (entry.getValue()+.0d)/treecount) {	
					results.add(0,processOneSet(taxonIdentifier, keepclusters, clusters));
					ti--;
					if (ti < 0) {
						break;
					}
					threshold = thresholds[ti];
				}
				clusters.add(entry.getKey());
			}
			while (ti >= 0) {
				results.add(0, processOneSet(taxonIdentifier, keepclusters, clusters));
				ti--;
			}
		}
	}

	protected Object processOneSet(TaxonIdentifier taxonIdentifier, boolean keepclusters,
			List<STITreeCluster> clusters) {
		Object buildTreeFromClusters = Utils.buildTreeFromClusters(clusters, taxonIdentifier, keepclusters);
		return buildTreeFromClusters;
	}


	protected int countClusteres(Iterable<Tree> trees, TaxonIdentifier taxonIdentifier,
			HashMap<STITreeCluster, Integer> count) {
		int treecount = 0;
		for (Tree tree : trees) {
			treecount++;
			List<STITreeCluster> geneClusters = Utils.getGeneClusters(tree, taxonIdentifier); //taxoncount changes
			for (STITreeCluster cluster: geneClusters) {

				if (count.containsKey(cluster)) {
					count.put(cluster, count.get(cluster) + 1);
					continue;
				}
				STITreeCluster comp = cluster.complementaryCluster();
				if (count.containsKey(comp)) {
					count.put(comp, count.get(comp) + 1);
					continue;
				}
				count.put(cluster, 1);
			}
		}
		return treecount;
	}

}
