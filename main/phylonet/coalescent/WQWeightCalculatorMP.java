package phylonet.coalescent;


import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;

/**
 * Knows how to compute the score of a given tripartition
 * @author smirarab
 *
 */
class WQWeightCalculatorMP extends WQWeightCalculator{

	private TraversalWeightCalculator tmpalgorithm;
	
	private LinkedBlockingQueue<Long> queue;
	private boolean threadingOff = false;

	public WQWeightCalculatorMP(AbstractInference<Tripartition> inference, LinkedBlockingQueue<Long> queue2) {
		super(inference);
		this.queue = queue2;
		this.lastTime = System.currentTimeMillis();
		this.dataCollection = (WQDataCollectionMP) inference.dataCollection;
		this.inference = (WQInferenceConsumerMP) inference;

		//this.algorithm = new TraversalWeightCalculator();
		this.algorithm = new CondensedTraversalWeightCalculator();
		tmpalgorithm = new TraversalWeightCalculator();
		//tmpalgorithm.setupGeneTrees((WQInference) inference);


	}
	
	
	public int getCalculatedWeightCount() {
		return this.callcounter;
}
	
	public Long calculateWeight(Tripartition t) {
		return this.calculateWeight(new Tripartition[] {(Tripartition) t})[0];
	}

	public Long getWeight(Tripartition t) {
		this.callcounter ++;
		Long weight = getCalculatedWeight(t);
		if (weight != null) {
			return weight;
		}

		if(isThreadingOff()) { // After the main DP, for computing final score. 		
			weight =  calculateWeight(t);
			saveWeight(t, weight);
			return weight;
		}

		try {
			weight = queue.take();
		}
		catch(Exception e) {
			throw new RuntimeException(e);
		}
		if(weight == TurnTaskToScores.THEEND) {// After the main DP, for computing final score, we switch to local calculator
			setThreadingOff(true);
			weight =  calculateWeight(t);
		}
		saveWeight(t, weight);

		return weight;

	}


	public boolean isThreadingOff() {
		return threadingOff;
	}

	public void setThreadingOff(boolean done) {
		this.threadingOff = done;

	}
	

	/**
	 * one of ASTRAL-III way of calculating weights
	 * Should be memory efficient
	 * @author chaoszhang
	 *
	 */
	class CondensedTraversalWeightCalculator extends WeightCalculatorAlgorithm {
		Polytree polytree;

		@Override
		Long[] calculateWeight(Tripartition[] trip) {
			return polytree.WQWeightByTraversal(trip);
		}
		
		@Override
		Long calculateWeight(Tripartition t) {
			return polytree.WQWeightByTraversal(t);
		}

		/***
		 * Each gene tree is represented as a list of integers, using positive numbers
		 * for leaves, where the number gives the index of the leaf. 
		 * We use negative numbers for internal nodes, where the value gives the number of children. 
		 * Minus infinity is used for separating different genes. 
		 */
		@Override
		void setupGeneTrees(WQInference inference) {
			polytree = new Polytree(inference.trees, (WQDataCollectionMP) dataCollection);
		}


	}


	/**
	 * ASTRAL-II way of calculating weights 
	 * @author smirarab
	 *
	 */
	class TraversalWeightCalculator extends WQWeightCalculator.TraversalWeightCalculator {

		public int maxHeight;

		/***
		 * Each gene tree is represented as a list of integers, using positive numbers
		 * for leaves, where the number gives the index of the leaf. 
		 * We use negative numbers for internal nodes, where the value gives the number of children. 
		 * Minus infinity is used for separating different genes. 
		 */
		@Override
		void setupGeneTrees(WQInference inference) {
			Logging.log("Using tree-based weight calculation.");
			List<Integer> temp = new ArrayList<Integer>(); 

			Stack<Integer> stackHeight = new Stack<Integer>();
			maxHeight = 0;
			for (Tree tr :  inference.trees) {
				/**
				 * Traverse tree and 1) build geneTreesAsInts, 2) compute maxHeight
				 */

				for (TNode node : tr.postTraverse()) {
					if (node.isLeaf()) {                        
						temp.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
						stackHeight.push(0);
					} else {
						temp.add(-node.getChildCount());
						int h = 0;
						for (int i = 0; i < node.getChildCount(); i++) {
							int childheight = stackHeight.pop();
							if(childheight > h)
								h = childheight;
						}
						h++;
						stackHeight.push(h);
					}
					if (node.isRoot()) {
						temp.add(Integer.MIN_VALUE);
						stackHeight.clear();
					}
					if(stackHeight.size()>maxHeight) {
						maxHeight = stackHeight.size();
					}
				}

				//Logging.log(tr);
			}
			geneTreesAsInts = new int[temp.size()];
			int i = 0;
			for (int v : temp) {
				geneTreesAsInts[i++] = v;
			}

		}

		public int[] geneTreesAsInts() {

			return this.geneTreesAsInts;
		}

	}

	
	/**
	 * Each algorithm will have its own data structure for gene trees
	 * @param wqInference
	 */
	@Override
	public void setupGeneTrees(WQInference wqInference) {
		tmpalgorithm.setupGeneTrees((WQInferenceConsumerMP) wqInference);
		this.algorithm.setupGeneTrees(wqInference);
	}

	//TODO: this is algorithm-specific should not be exposed. Fix.
	public int[] geneTreesAsInts() {
		return (tmpalgorithm).geneTreesAsInts;
	}
	@Override
	public Long[] calculateWeight(Tripartition[] t) {
		return this.algorithm.calculateWeight(t);
	}

	//TODO: this is algorithm-specific should not be exposed. Fix. 
	public int maxHeight() {
		return ((TraversalWeightCalculator)tmpalgorithm).maxHeight;
	}




}
