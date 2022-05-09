package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import phylonet.coalescent.Utils.ClusterComparator;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

public class GreedyConsensusMP extends GreedyConsensus{
	
	protected GreedyConsensusMP() {};
	
	public static GreedyConsensus instance = new GreedyConsensusMP();
	
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
	@Override
    public Collection<Tree> greedyConsensus(Iterable<Tree> trees, 
    		double[] thresholds, boolean randomzie, int repeat, 
    		TaxonIdentifier taxonIdentifier, boolean keepclusters) {
    	Logging.logTimeMessage("Utils 219-222: ");
			
    	List<Tree> outTrees = new ArrayList<Tree>();
        HashMap<STITreeCluster, Integer> count = new HashMap<STITreeCluster, Integer>();
        int treecount = countClusteres(trees, taxonIdentifier, count);
        
        Logging.logTimeMessage("Utils 240-243: " );
			
        ArrayList<Future<Tree>> futures = new ArrayList<Future<Tree>>();
        for (int gi = 0; gi < repeat; gi++) {
        	TreeSet<Entry<STITreeCluster,Integer>> countSorted = new 
        			TreeSet<Entry<STITreeCluster,Integer>>(new ClusterComparator(randomzie, taxonIdentifier.taxonCount()));
        
	        countSorted.addAll(count.entrySet());
	        
	        int ti = thresholds.length - 1;
	        double threshold = thresholds[ti];
	        List<STITreeCluster> clusters = new ArrayList<STITreeCluster>();   
	        for (Entry<STITreeCluster, Integer> entry : countSorted) {
	        	if (threshold > (entry.getValue()+.0d)/treecount) {	
	        		List<STITreeCluster> clusterCopy = new ArrayList<STITreeCluster>(clusters);
	        		futures.add(Threading.submit(new greedyConsensusLoop(taxonIdentifier, keepclusters, clusterCopy)));
	        		ti--;
	        		if (ti < 0) {
	        			break;
	        		}
	        		threshold = thresholds[ti];
	        	}
	    		clusters.add(entry.getKey());
	        }
	        while (ti >= 0) {
        		List<STITreeCluster> clusterCopy = new ArrayList<STITreeCluster>(clusters);
        		futures.add(Threading.submit(new greedyConsensusLoop(taxonIdentifier, keepclusters, clusterCopy)));
	    		ti--;
	        }
        }
        for(int i = 0; i < futures.size(); i++) {
        	try {
				outTrees.add(0, futures.get(i).get());
			} catch (InterruptedException e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			} catch (ExecutionException e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
        }
        Logging.logTimeMessage("Utils 269-272: " );
			
        return outTrees;
    }
	
    public static class greedyConsensusLoop implements Callable<Tree>{
		TaxonIdentifier taxonIdentifier;
		boolean keepclusters;
		List<STITreeCluster> clusters;
		public greedyConsensusLoop(TaxonIdentifier taxonIdentifier, boolean keepclusters, List<STITreeCluster> clusters) {
			this.keepclusters = keepclusters;
			this.taxonIdentifier = taxonIdentifier;
			this.clusters = clusters;
		}
		public Tree call() {
			return Utils.buildTreeFromClusters(clusters, taxonIdentifier, keepclusters);
		}
	}

}
