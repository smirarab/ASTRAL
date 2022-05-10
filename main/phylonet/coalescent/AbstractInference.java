package phylonet.coalescent;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;
import phylonet.tree.util.Collapse;

/***
 * Type T corresponds to a tripartition in ASTRAL
 * @author smirarab
 *
 * @param <T>
 */
public abstract class AbstractInference<T> {


	public List<Tree> trees;
	public List<Tree> extraTrees = null;
	public List<Tree> toRemoveExtraTrees = null;
	protected boolean removeExtraTree;


	Collapse.CollapseDescriptor cd = null;
	
	public AbstractDataCollection<T> dataCollection;
	public AbstractWeightCalculator<T> weightCalculator;

	public Options options;
	protected DecimalFormat df;
	
	double estimationFactor = 0;
	
	public AbstractInference(Options options, List<Tree> trees,
			List<Tree> extraTrees, List<Tree> toRemoveExtraTrees) {
		super();
		this.options = options;
		this.trees = trees;
		this.extraTrees = extraTrees;
		this.removeExtraTree = options.isRemoveExtraTree();
		this.toRemoveExtraTrees = toRemoveExtraTrees;
		
		this.initDF();

	}


	/**
	 * Sets up data structures before starting DP
	 */
	public void setup() {
		this.setupSearchSpace();
		this.initializeWeightCalculator();
		this.setupMisc();
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
	
	protected Collapse.CollapseDescriptor doCollapse(List<Tree> trees) {
		Collapse.CollapseDescriptor cd = Collapse.collapse(trees);
		return cd;
	}

	protected void restoreCollapse(List<Solution> sols, Collapse.CollapseDescriptor cd) {
		for (Solution sol : sols) {
			Tree tr = sol.getTree();
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

		Logging.log("Number of taxa: " + GlobalMaps.taxonIdentifier.taxonCount()+
		        " (" + GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesCount() +" species)");
		Logging.log("Taxa: " + GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesNames());
		Logging.log("Taxon occupancy: " + taxonOccupancy.toString());
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
	protected List<Solution> findTreesByDP(IClusterCollection clusters) {

		Vertex all = (Vertex) clusters.getTopVertex();

		Logging.log("Size of largest cluster: " +all.getCluster().getClusterSize());

		AbstractComputeMinCostTask<T> allTask = newComputeMinCostTask(this,all,clusters);
		allTask.compute();
		
		List<Solution> solutions = processSolutions(all);

		return (List<Solution>) (List<Solution>) solutions;
	}


	List<Solution> processSolutions( Vertex all) {
		List<Solution> solutions = new ArrayList<Solution>();
		
		try {
			if ( all._max_score == Integer.MIN_VALUE) {
				throw new CannotResolveException(all.getCluster().toString());
			}
		} catch (CannotResolveException e) {
			Logging.log("Was not able to build a fully resolved tree. Not" +
					"enough clusters present in input gene trees ");
			e.printStackTrace();
			System.exit(1);
		}
		
		Logging.logTimeMessage("AbstractInference 193: " );

		Logging.log("Total Number of elements weighted: "+ countWeights());

		List<STITreeCluster> minClusters = new LinkedList<STITreeCluster>();
		//List<Double> coals = new LinkedList<Double>();
		Stack<Vertex> minVertices = new Stack<Vertex>();
		if (all._min_rc != null) {
			minVertices.push(all._min_rc);
		}
		if (all._min_lc != null) {
			minVertices.push(all._min_lc);
		}
		SpeciesMapper spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
		while (!minVertices.isEmpty()) {
			Vertex pe =  minVertices.pop();
			STITreeCluster stCluster = spm.
					getSTClusterForGeneCluster(pe.getCluster());
			//System.out.println(pe._min_rc);
			//System.out.println(pe._min_lc);
			minClusters.add(stCluster);
			//System.out.println(pe.getCluster().getClusterSize()+"\t"+pe._max_score);
			// int k = sigmaNs/(stTaxa.length-1);

			if ( !GlobalMaps.taxonNameMap.getSpeciesIdMapper().isSingleSP(pe.getCluster().getBitSet()) && (pe._min_lc == null || pe._min_rc == null))
				Logging.log("hmm; this shouldn't have happened: "+ pe);
			
			if (pe._min_rc != null) {
				minVertices.push(pe._min_rc);
			}
			if (pe._min_lc != null) {
				minVertices.push(pe._min_lc);
			}
//			if (pe._min_lc != null && pe._min_rc != null) {
//				coals.add(pe.get_c());
//			} else {
//				coals.add(0D);
//			}
//			if (pe.getSubcl() != null) {
//				for (Vertex v : pe.getSubcl()) {
//					minVertices.push(v);
//				}
//			}
		}
		Solution sol = new Solution();
		if ((minClusters == null) || (minClusters.isEmpty())) {
			Logging.log("WARN: empty minClusters set.");
			STITree<Double> tr = new STITree<Double>();
			for (String s : GlobalMaps.taxonIdentifier.getAllTaxonNames()) {
				((MutableTree) tr).getRoot().createChild(s);
			}
			sol.setTree(tr);
		} else {
			sol.setTree( Utils.buildTreeFromClusters(minClusters, spm.getSTTaxonIdentifier(), false) );
		}

		Long cost = getTotalCost(all);
		sol.setCoalNum(cost);
		solutions.add(sol);
        Logging.log("Optimization score: " + cost);
		return solutions;
	}
	
	protected abstract void initializeWeightCalculator();
	
	public int countWeights() {
		return weightCalculator.getCalculatedWeightCount();
	}

	/***
	 * Creates the set X 
	 */
	private void setupSearchSpace() {
		long startTime = System.currentTimeMillis();

		mapNames();

		dataCollection = Factory.instance.newCounter(newClusterCollection(), this);
		weightCalculator = newWeightCalculator();

		/**
		 * Fors the set X by adding from gene trees and
		 * by adding using ASTRAL-II hueristics
		 */
		dataCollection.formSetX(this);

		
		if (options.isExactSolution()) {
	          Logging.log("calculating all possible bipartitions ...");
		    dataCollection.addAllPossibleSubClusters(this.dataCollection.clusters.getTopVertex().getCluster());
		}

	      
		if (extraTrees != null && extraTrees.size() > 0 && options.getAddExtra() != 3) {		
	        Logging.log("calculating extra bipartitions from extra input trees ...");
			dataCollection.addExtraBipartitionsByInput(extraTrees,options.isExtrarooted());
			int s = this.dataCollection.clusters.getClusterCount();
			/*
			 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
			 */
			Logging.log("Number of Clusters after additions from extra trees: "
					+ s);
		}
		
		if (toRemoveExtraTrees != null && toRemoveExtraTrees.size() > 0 && this.removeExtraTree) {		
	        Logging.log("Removing extra bipartitions from extra input trees ...");
			dataCollection.removeExtraBipartitionsByInput(toRemoveExtraTrees,true);
			int s = this.dataCollection.clusters.getClusterCount();
			/*
			 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
			 */
			Logging.log("Number of Clusters after deletion of extra tree bipartitions: "
					+ s);
		}
		
		if (this.options.isOutputSearchSpace()) {
			for (Set<Vertex> s: dataCollection.clusters.getSubClusters()) {
				for (Vertex v : s) {
					System.out.println(v.getCluster());
				}
			}
		}

		Logging.logTimeMessage("" );

		Logging.log("partitions formed in "
			+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");

		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		
		// Obsolete 
		weightCalculator.preCalculateWeights(trees, extraTrees);
		

		Logging.log("Dynamic Programming starting after "
				+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
		
	}
	public abstract void setupMisc();

	public List<Solution> inferSpeciesTree() {

		List<Solution> solutions;

		solutions = findTreesByDP(this.dataCollection.clusters);

		return (List<Solution>) solutions;
	}

	public abstract IClusterCollection newClusterCollection();
	
	public abstract AbstractWeightCalculator<T> newWeightCalculator();

	public abstract AbstractComputeMinCostTask<T> newComputeMinCostTask(AbstractInference<T> iInference,
			Vertex all, IClusterCollection clusters);
	
	public abstract Long getTotalCost(Vertex all);
	
	public double getDLbdWeigth() {
		return options.getDLbdWeigth();
	}

	
	public double getCS() {
		return options.getCS();
	}

	public double getCD() {
		return options.getCD();
	}

	public int getBranchAnnotation() {
		return this.options.getBranchannotation();
	}


	public void setDLbdWeigth(double d) {
		options.setDLbdWeigth(d);
	}

//	protected Object semiDeepCopy() {
//		try {
//			AbstractInference<T> clone =  (AbstractInference<T>) super.clone();
//			clone.dataCollection = (AbstractDataCollection<T>) this.dataCollection.clone();
//			clone.weightCalculator = (AbstractWeightCalculator<T>) this.weightCalculator.clone();
//			return clone;
//		} catch (CloneNotSupportedException e) {
//			e.printStackTrace();
//			throw new RuntimeException("unexpected error");
//		}
//	}
}
