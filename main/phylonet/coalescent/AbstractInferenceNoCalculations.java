package phylonet.coalescent;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.ConcurrentLinkedQueue;

import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.util.Collapse;

public abstract class AbstractInferenceNoCalculations<T> extends AbstractInference<T> {
	
	ConcurrentLinkedQueue<ICalculateWeightTask<Tripartition>> queue1;
	
	public AbstractInferenceNoCalculations(Options options, List<Tree> trees,
			List<Tree> extraTrees) {
		super(options, trees, extraTrees);
		df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		DecimalFormatSymbols dfs = DecimalFormatSymbols.getInstance();
		dfs.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(dfs);
		
	}
	
	public AbstractInferenceNoCalculations(AbstractInference in) {
		super(in.options, in.trees, in.extraTrees);
		this.dataCollection = in.dataCollection;
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
			Tree tr = sol._st;
			Collapse.expand(cd, (MutableTree) tr);
			for (TNode node : tr.postTraverse())
				if (((STINode) node).getData() == null)
					((STINode<Integer>) node).setData(Integer.valueOf(0));
		}
	}

	private int getResolutionsNumber(int nodeNumber) {
		int total = 1;
		for (int i = 3; i <= nodeNumber; i++) {
			total *= (2 * i - 3);
		}
		return total;
	}

	public void mapNames() {
		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}
        for (Tree tr : trees) {
            String[] leaves = tr.getLeaves();
            for (int i = 0; i < leaves.length; i++) {
                GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
            }
        }
        
        GlobalMaps.taxonNameMap.checkMapping(trees);

		System.err.println("Number of taxa: " + GlobalMaps.taxonIdentifier.taxonCount()+
		        " (" + GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesCount() +" species)"
		);
		System.err.println("Taxa: " + GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesNames());
	}

	public abstract double scoreGeneTree(Tree scorest, boolean initialize) ;

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
		
		
		//readScoreFromQueue scoretaker = new readScoreFromQueue(this, all, clusters, queue2);
		
		//TurnTaskToScores scorer = new TurnTaskToScores(weightCalculator, queue1, queue2);
		
		//writeTasktoQueue taskmaker = new writeTasktoQueue(this, all, clusters, queue1, scorer);

		//(new Thread(taskmaker)).start();
		try {
			//vertexStack.push(all);
			
			AbstractComputeMinCostTask<T> allTask = newComputeMinCostTask(this,all,clusters,true);
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
		
		System.err.println("Total Number of elements weighted: "+ weightCalculatorNoCalculations.getCalculatedWeightCount());

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
			Vertex pe = minVertices.pop();
			STITreeCluster stCluster = spm.
					getSTClusterForGeneCluster(pe.getCluster());
			//System.out.println(pe._min_rc);
			//System.out.println(pe._min_lc);
			minClusters.add(stCluster);
			//System.out.println(pe.getCluster().getClusterSize()+"\t"+pe._max_score);
			// int k = sigmaNs/(stTaxa.length-1);

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
			sol._st = Utils.buildTreeFromClusters(minClusters, spm.getSTTaxonIdentifier());
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
        System.err.println("Final optimization score: " + cost);

		return (List<Solution>) (List<Solution>) solutions;
	}
	
	public void setupSearchSpace() {
		long startTime = System.currentTimeMillis();

		mapNames();

		dataCollection = newCounter(newClusterCollection());
		weightCalculator = newWeightCalculator();
		weightCalculatorNoCalculations = newWeightCalculatorNoCalculations();
		dataCollection.computeTreePartitions(this);

		if (this.getAddExtra() != 0) {
		    System.err.println("calculating extra bipartitions to be added at level " + this.getAddExtra() +" ...");
		    dataCollection.addExtraBipartitionByExtension(this);
		}
		
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

		System.err.println("partitions formed in "
			+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
		
		weightCalculator.preCalculateWeights(trees, extraTrees);

		System.err.println("Dynamic Programming starting after "
				+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");
		
	}
	public List<Solution> inferSpeciesTree() {
		weightCalculatorNoCalculations = newWeightCalculatorNoCalculations();

		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		
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
	abstract AbstractWeightCalculatorNoCalculations<T> newWeightCalculatorNoCalculations();

	abstract AbstractComputeMinCostTask<T> newComputeMinCostTask(AbstractInference<T> dlInference,
			Vertex all, IClusterCollection clusters, boolean isWriteToQueue);
	
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

}
