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
						
        ArrayList futures = new ArrayList();
        List<Tree> outTrees =  new ArrayList<Tree>();
        
        super.greedyConsensus(trees, thresholds, randomzie, repeat, taxonIdentifier, keepclusters, futures);
        
        for(int i = 0; i < futures.size(); i++) {
        	try {
				outTrees.add(((Future<Tree>) futures.get(i)).get());
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
	
	protected Object processOneSet(TaxonIdentifier taxonIdentifier, boolean keepclusters,
			List<STITreeCluster> clusters) {
		List<STITreeCluster> clusterCopy = new ArrayList<STITreeCluster>(clusters);
		return Threading.submit(new greedyConsensusLoop(taxonIdentifier, keepclusters, clusterCopy));
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
