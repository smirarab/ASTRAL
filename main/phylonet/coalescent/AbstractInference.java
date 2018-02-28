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
import java.util.concurrent.PriorityBlockingQueue;

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
	
	//protected boolean rooted = true;
	//protected boolean extrarooted = true;
	protected List<Tree> trees;
	protected List<Tree> extraTrees = null;
	//protected boolean exactSolution;
	
	//protected String[] gtTaxa;
	//protected String[] stTaxa;

	Collapse.CollapseDescriptor cd = null;
	
	AbstractDataCollection<T> dataCollection;
	AbstractWeightCalculator<T> weightCalculator;
	AbstractWeightCalculatorNoCalculations<T> weightCalculatorNoCalculations;
	Boolean done = false;
//	private int addExtra;
//	public boolean outputCompleted;
//	boolean searchSpace;
//	private boolean run;
	protected Options options;
	DecimalFormat df;
	

	LinkedBlockingQueue<Long> queue2;
	public LinkedBlockingQueue<Iterable<VertexPair>> queue4;

	double estimationFactor = 0;
	

	public AbstractInference(Options options, List<Tree> trees,
			List<Tree> extraTrees) {
		super();
		this.options = options;
		this.trees = trees;
		this.extraTrees = extraTrees;
		
		df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		DecimalFormatSymbols dfs = DecimalFormatSymbols.getInstance();
		dfs.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(dfs);
		
	}

	public boolean isRooted() {
		return options.isRooted();
	}
	
	public Iterable<VertexPair> getClusterResolutions2() {
			Iterable<VertexPair> ret = null;
			try{
				ret = queue4.take();
			}
			catch (Exception e) {
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
		
		
		//if (CommandLine._print) {
			//System.err.println("Weights are: "
				//	+ counter.weights);
		//}
		//System.out.println("domination calcs:" + counter.cnt);
		if(CommandLine.timerOn) {
	       	System.err.println("TIME TOOK FROM LAST NOTICE AbstractInference 193: " + (double)(System.nanoTime()-CommandLine.timer)/1000000000);
			CommandLine.timer = System.nanoTime();
		}
		System.err.println("Total Number of elements weighted: "+ weightCalculator.getCalculatedWeightCount());

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
		if(CommandLine.timerOn) {
	       	System.err.println("TIME TOOK FROM LAST NOTICE AbstractInference 283: " + (double)(System.nanoTime()-CommandLine.timer)/1000000000);
			CommandLine.timer = System.nanoTime();
		}
        System.err.println("Final optimization score: " + cost);
        
		return (List<Solution>) (List<Solution>) solutions;
	}
	
	/**
	 * Sets up data structures before starting DP
	 */
	void setup() {
		this.setupSearchSpace();
		this.initializeWeightCalculator();
		this.setupMisc();
	}
	
	abstract void initializeWeightCalculator();

	

	/***
	 * Creates the set X 
	 */
	private void setupSearchSpace() {
		long startTime = System.currentTimeMillis();

		mapNames();

		dataCollection = newCounter(newClusterCollection());
		weightCalculator = newWeightCalculator();
		
		/**
		 * Fors the set X by adding from gene trees and
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

		    //johng23
	   if(CommandLine.timerOn) {
			System.err.println("TIME TOOK FROM LAST NOTICE: " + (double)(System.nanoTime()-CommandLine.timer)/1000000000);
			CommandLine.timer = System.nanoTime();
		    }
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
	
	abstract void setupMisc();

	public List<Solution> inferSpeciesTree() {
		
		List<Solution> solutions;		
		
		solutions = findTreesByDP(this.dataCollection.clusters);

/*		if (GlobalMaps.taxonNameMap == null && rooted && extraTrees == null && false) {
			restoreCollapse(solutions, cd);
			}*/

		return (List<Solution>) solutions;
	}

	abstract IClusterCollection newClusterCollection();
	
	abstract AbstractDataCollection<T> newCounter(IClusterCollection clusters);
	
	abstract AbstractWeightCalculator<T> newWeightCalculator();

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
	
	protected Object semiDeepCopy() {
		try {
			AbstractInference<T> clone =  (AbstractInference<T>) super.clone();
			clone.dataCollection = (AbstractDataCollection<T>) this.dataCollection.clone();
			clone.weightCalculator = (AbstractWeightCalculator<T>) this.weightCalculator.clone();
			return clone;
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
			throw new RuntimeException("unexpected error");
		}
	}
}
