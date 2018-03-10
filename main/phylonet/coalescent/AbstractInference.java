package phylonet.coalescent;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.util.Collapse;


/***
 * Type T corresponds to a tripartition in ASTRAL
 * @author smirarab
 *
 * @param <T>
 */
public abstract class AbstractInference<T> implements Cloneable{
	
	protected List<Tree> trees;
	protected List<Tree> extraTrees = null;


	Collapse.CollapseDescriptor cd = null;
	
	AbstractDataCollection<T> dataCollection;
	AbstractWeightCalculatorTask<T> weightCalculator;
	AbstractWeightCalculatorProducer<T> weightCalculatorNoCalculations;
	Boolean done = false;
	protected Options options;
	DecimalFormat df;
	
	private LinkedBlockingQueue<Long> queueWeightResults;
	private LinkedBlockingQueue<Iterable<VertexPair>> queueClusterResolutions;

	double estimationFactor = 0;
	

	public AbstractInference(Options options, List<Tree> trees,
			List<Tree> extraTrees) {
		super();
		this.options = options;
		this.trees = trees;
		this.extraTrees = extraTrees;
		
		this.initDF();
		
	}

	private void initDF() {
		df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		DecimalFormatSymbols dfs = DecimalFormatSymbols.getInstance();
		dfs.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(dfs);
	}
	
	public boolean isRooted() {
		return options.isRooted();
	}
	
	public Iterable<VertexPair> getClusterResolutions(Vertex v) {
			Iterable<VertexPair> ret = null;
			try{
				ret = getQueueClusterResolutions().take();
			}
			catch (Exception e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
			return ret;

	}
	
	protected Collapse.CollapseDescriptor doCollapse(List<Tree> trees) {
		Collapse.CollapseDescriptor cd = Collapse.collapse(trees);
		return cd;
	}

	protected void restoreCollapse(List<Solution> sols, Collapse.CollapseDescriptor cd) {
		for (Solution sol : sols) {
			Tree tr = sol._st;
			Collapse.expand(cd, (MutableTree) tr);
			for (TNode node : tr.postTraverse())
				if (((STINode) node).getData() == null)
					((STINode) node).setData(Integer.valueOf(0));
		}
	}

	private int getResolutionsNumber(int nodeNumber) {
		int total = 1;
		for (int i = 3; i <= nodeNumber; i++) {
			total *= (2 * i - 3);
		}
		return total;
	}

	//TODO: Check whether this is in the right class
	public void mapNames() {
		HashMap<String, Integer> taxonOccupancy = new HashMap<String, Integer>();
		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}
        for (Tree tr : trees) {
            String[] leaves = tr.getLeaves();
            for (int i = 0; i < leaves.length; i++) {
                GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
                taxonOccupancy.put(leaves[i], Utils.increment(taxonOccupancy.get(leaves[i])));
            }
        }
        
        GlobalMaps.taxonNameMap.checkMapping(trees);

		System.err.println("Number of taxa: " + GlobalMaps.taxonIdentifier.taxonCount()+
		        " (" + GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesCount() +" species)"
		);
		System.err.println("Taxa: " + GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesNames());
		System.err.println("Taxon occupancy: " + taxonOccupancy.toString());
	}

	/***
	 * Scores a given tree. 
	 * @param scorest
	 * @param initialize
	 * @return
	 */
	public abstract double scoreSpeciesTreeWithGTLabels(Tree scorest, boolean initialize) ;

	/***
	 * This implements the dynamic programming algorithm
	 * @param clusters
	 * @return
	 */
	List<Solution> findTreesByDP(IClusterCollection clusters) {
		List<Solution> solutions = new ArrayList<Solution>();
		/*
		 * clusterToVertex = new HashMap<STITreeCluster, Vertex>(); for
		 * (Set<Vertex> vs: clusters.values()) { for (Vertex vertex : vs) {
		 * clusterToVertex.put(vertex._cluster,vertex); } } Vertex all =
		 * (Vertex) clusters.get(Integer .valueOf(stTaxa.length)).toArray()[0];
		 * computeMinCost(clusters, all, sigmaN, counter,trees, taxonMap);
		 * 
		 * System.out.println("first round finished, adding new STBs");
		 * counter.addExtraBipartitions(clusters, stTaxa);
		 */
/*		clusterToVertex = new HashMap<STITreeCluster, Vertex>(sigmaNs);
		for (Set<Vertex> vs : clusters.values()) {
			for (Vertex vertex : vs) {
				vertex._max_score = -1;
				clusterToVertex.put(vertex._cluster, vertex);
			}
		}
*/
		Vertex all = (Vertex) clusters.getTopVertex();

		System.err.println("Size of largest cluster: " +all.getCluster().getClusterSize());
		
		try {
			//vertexStack.push(all);
			AbstractComputeMinCostTask<T> allTask = newComputeMinCostTask(this,all);
			//ForkJoinPool pool = new ForkJoinPool(1);
			allTask.compute();
			double v = all._max_score;
			if (v == Integer.MIN_VALUE) {
				throw new CannotResolveException(all.getCluster().toString());
			}
		} catch (CannotResolveException e) {
			System.err.println("Was not able to build a fully resolved tree. Not" +
					"enough clusters present in input gene trees ");
			e.printStackTrace();
			System.exit(1);
		}
		
		Logging.logTimeMessage("AbstractInference 193: " );
		
		System.err.println("Total Number of elements: " + countWeights());

		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		List<Double> coals = new LinkedList<Double>();
		Stack<Vertex> minVertices = new Stack<Vertex>();
		if (all._min_rc != null) {
			minVertices.push(all._min_rc);
		}
		if (all._min_lc != null) {
			minVertices.push(all._min_lc);
		}
		if (all._subcl != null) {
			for (Vertex v : all._subcl) {
				minVertices.push(v);
			}
		}		
		SpeciesMapper spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
		while (!minVertices.isEmpty()) {
			Vertex pe = (Vertex) minVertices.pop();
			STITreeCluster stCluster = spm.
					getSTClusterForGeneCluster(pe.getCluster());
			//System.out.println(pe._min_rc);
			//System.out.println(pe._min_lc);
			minClusters.add(stCluster);
			//System.out.println(pe.getCluster().getClusterSize()+"\t"+pe._max_score);
			// int k = sigmaNs/(stTaxa.length-1);

			if ( !GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleSP(pe.getCluster().getBitSet()) && (pe._min_lc == null || pe._min_rc == null))
				System.err.println("hmm; this shouldn't have happened: "+ pe);
			
			if (pe._min_rc != null) {
				minVertices.push(pe._min_rc);
			}
			if (pe._min_lc != null) {
				minVertices.push(pe._min_lc);
			}
			if (pe._min_lc != null && pe._min_rc != null) {
				coals.add(pe._c);
			} else {
				coals.add(0D);
			}
			if (pe._subcl != null) {
				for (Vertex v : pe._subcl) {
					minVertices.push(v);
				}
			}
		}
		Solution sol = new Solution();
		if ((minClusters == null) || (minClusters.isEmpty())) {
			System.err.println("WARN: empty minClusters set.");
			STITree<Double> tr = new STITree<Double>();
			for (String s : GlobalMaps.taxonIdentifier.getAllTaxonNames()) {
				((MutableTree) tr).getRoot().createChild(s);
			}
			sol._st = tr;
		} else {
			sol._st = Utils.buildTreeFromClusters(minClusters, spm.getSTTaxonIdentifier(), false);
		}

		/* HashMap<TNode,BitSet> map = new HashMap<TNode,BitSet>();
		for (TNode node : sol._st.postTraverse()) {
			BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
			if (node.isLeaf()) {
				bs.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
				map.put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = map.get(child);
					bs.or(childCluster);
				}
				map.put(node, bs);
			}
//            System.err.println("Node: "+node);
			STITreeCluster c = new STITreeCluster();
			c.setCluster(bs);
//            System.err.println("m[0]: "+((STITreeCluster)minClusters.get(0)).toString2());
//            System.err.println("C: "+c.toString2());
//            System.err.println("Equals: "+((STITreeCluster)minClusters.get(0)).equals(c));
			if (c.getClusterSize() == GlobalMaps.taxonIdentifier.taxonCount()) {
				((STINode<Double>) node).setData(Double.valueOf(0));
			} else {
				int pos = minClusters.indexOf(c);                                
				((STINode<Double>) node).setData((Double) coals.get(pos));
			}
		}*/

		Long cost = getTotalCost(all);
		sol._totalCoals = cost;
		solutions.add(sol);
		Logging.logTimeMessage("AbstractInference 283: ");
			
        System.err.println("Final optimization score: " + cost);
        
		return (List<Solution>) (List<Solution>) solutions;
	}

	public int countWeights() {
		return weightCalculator.getCalculatedWeightCount();
	}
	
	/**
	 * Sets up data structures before starting DP
	 */
	void setup() {
		this.setupSearchSpace();
		this.initializeWeightCalculator();
		this.setupMisc();
	}
	
	/***
	 * Creates the set X 
	 */
	private void setupSearchSpace() {
		long startTime = System.currentTimeMillis();

		mapNames();

		dataCollection = newCounter(newClusterCollection());
		weightCalculator = newWeightCalculator();
		
		/**
		 * Forms the set X by adding from gene trees and
		 * by adding using ASTRAL-II hueristics
		 */
		dataCollection.formSetX(this); //TODO: is this necessary. 

		
		if (options.isExactSolution()) {
	          System.err.println("calculating all possible bipartitions ...");
		    dataCollection.addAllPossibleSubClusters(this.dataCollection.clusters.getTopVertex().getCluster());
		}

	      
		if (extraTrees != null && extraTrees.size() > 0) {		
	        System.err.println("calculating extra bipartitions from extra input trees ...");
			dataCollection.addExtraBipartitionsByInput(extraTrees,options.isExtrarooted());
			int s = this.dataCollection.clusters.getClusterCount();
			/*
			 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
			 */
			System.err.println("Number of Clusters after additions from extra trees: "
					+ s);
		}
		
		if (this.options.isOutputSearchSpace()) {
			for (Set<Vertex> s: dataCollection.clusters.getSubClusters()) {
				for (Vertex v : s) {
					System.out.println(v.getCluster());
				}
			}
		}

		//counter.addExtraBipartitionsByHeuristics(clusters);

		Logging.logTimeMessage("" );
			
		System.err.println("partitions formed in "
			+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
		
		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		
		// Obsolete 
		weightCalculator.preCalculateWeights(trees, extraTrees);

		System.err.println("Dynamic Programming starting after "
				+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
		
	}
	
	public List<Solution> inferSpeciesTree() {
		
		List<Solution> solutions;		
		
		solutions = findTreesByDP(this.dataCollection.clusters);

/*		if (GlobalMaps.taxonNameMap == null && rooted && extraTrees == null && false) {
			restoreCollapse(solutions, cd);
			}*/

		return (List<Solution>) solutions;
	}

	protected Object semiDeepCopy() {
		try {
			AbstractInference<T> clone =  (AbstractInference<T>) super.clone();
			clone.dataCollection = (AbstractDataCollection<T>) this.dataCollection.clone();
			clone.weightCalculator = (AbstractWeightCalculatorConsumer<T>) this.weightCalculator.clone();
			return clone;
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
			throw new RuntimeException("unexpected error");
		}
	}

	abstract void initializeWeightCalculator();

	abstract void setupMisc();

	abstract IClusterCollection newClusterCollection();
	
	abstract AbstractDataCollection<T> newCounter(IClusterCollection clusters);
	
	abstract AbstractWeightCalculatorTask<T> newWeightCalculator();

	abstract AbstractComputeMinCostTask<T> newComputeMinCostTask(AbstractInference<T> dlInference,
			Vertex all);
	
	abstract Long getTotalCost(Vertex all);
	
	public double getDLbdWeigth() {
		return options.getDLbdWeigth();
	}

	
	public double getCS() {
		return options.getCS();
	}

	public double getCD() {
		return options.getCD();
	}

	
    public int getAddExtra() {
        return options.getAddExtra();
    }

	public int getBranchAnnotation() {
		return this.options.getBranchannotation();
	}

	public boolean shouldOutputCompleted() {
		
		return options.isOutputCompletedGenes();
	}

	public void setDLbdWeigth(double d) {
		options.setDLbdWeigth(d);
	}

	public LinkedBlockingQueue<Iterable<VertexPair>> getQueueClusterResolutions() {
		return queueClusterResolutions;
	}

	public void setQueueClusterResolutions(LinkedBlockingQueue<Iterable<VertexPair>> queueClusterResolutions) {
		this.queueClusterResolutions = queueClusterResolutions;
	}

	public LinkedBlockingQueue<Long> getQueueWeightResults() {
		return queueWeightResults;
	}

	public void setQueueWeightResults(LinkedBlockingQueue<Long> queueWeightResults) {
		this.queueWeightResults = queueWeightResults;
	}
}
