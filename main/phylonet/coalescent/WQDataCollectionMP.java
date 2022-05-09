package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;

/**
 * Sets up the set X
 * 
 * @author smirarab
 * 
 */
public class WQDataCollectionMP extends WQDataCollection
	implements Cloneable {

	static TaxonIdentifierMP taxid = (TaxonIdentifierMP) GlobalMaps.taxonIdentifier;

	public WQDataCollectionMP(HashClusterCollection clusters,
			AbstractInference inference) {
		super(clusters, inference);
	}

	/**
	 * Adds extra bipartitions added by user using the option -e and -f
	 */
	public void addExtraBipartitionsByInput(List<Tree> extraTrees,
			boolean extraTreeRooted) {

		// List<Tree> completedExtraGeeneTrees = new ArrayList<Tree>();
		ArrayList<Future> res = new ArrayList();
		for (Tree tr : extraTrees) {
			res.add(Threading.submit(new addExtraBipartitionByInputLoop(tr)));
		}
		
		for (Future f : res)
			try {
				f.get();
			} catch (InterruptedException | ExecutionException e) {
				throw new RuntimeException(e);
			}
	}

	@Override
	protected void secondRoundSampling(int secondRoundSampling, List<Tree>[] allGreedies, ArrayList<Tree> baseTrees) {
		CountDownLatch latch = new CountDownLatch(secondRoundSampling*allGreedies.length);
		for (int ii = 0; ii < secondRoundSampling; ii++) {
			for (int j = 0; j < allGreedies.length; j++) {

				Threading.execute(new FormSetXLoop(allGreedies[j].get(ii), baseTrees, latch));

			}			
			Logging.log("------------------------------");
		}
		try {
			latch.await();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Logging.logTimeMessage(" WQDataCollection 728-731: ");
	}

	public class FormSetXLoop extends WQDataCollection.FormSetXLoop {

		CountDownLatch latch;

		public FormSetXLoop(Tree tree,
				ArrayList<Tree> baseTrees, CountDownLatch latch) {
			super(tree,baseTrees);
			this.latch = latch;
		}

		public void run() {
			super.run();
			latch.countDown();
		}

	}
	
	@Override
	ArrayList addExtraBipartitionByHeuristics(Collection<Tree> contractedTrees,
			TaxonIdentifier tid, int polylimit) {

		ArrayList stringOutput = super.addExtraBipartitionByHeuristics(contractedTrees, tid, polylimit);
		
		for(int i = 0; i < stringOutput.size(); i++) {
			if(stringOutput.get(i) instanceof String) {
				System.err.print(stringOutput.get(i));
			} else {
				try {
					System.err.print(((Future<String>)stringOutput.get(i)).get());
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}
			}
		}
		return stringOutput;
	}
	
	@Override
	protected void myLog(Object log, ArrayList stringOutput) {
		stringOutput.add(log);
	}
	
	@Override
	protected Object invokeRunner(addExtraBipartitionByHeuristicsLoop callable) {
		return Threading.submit(callable);
	}


	public long[] getAllArray() {
		int counter = 0;
		int wordLength =  (taxid.taxonCount() / 64 + 1);
		long[] allArray = new long[this.treeAllClusters.size() * wordLength];
		for (int i = 0; i < this.treeAllClusters.size(); i++) {
			for (int j = wordLength - 1; j >= 0; j--)
				allArray[counter++] = this.treeAllClusters
						.get(i).getBitSet().words[j];
		}
		return allArray;
	}
}
