package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.util.Collapse;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;

public abstract class Inference<T> {

	protected boolean rooted = true;
	protected boolean extrarooted = true;
	protected List<Tree> trees;
	protected List<Tree> extraTrees = null;
	protected boolean exactSolution;
	
	protected String[] gtTaxa;
	//protected String[] stTaxa;

	Collapse.CollapseDescriptor cd = null;
	private double DLbdWeigth;
	private double CS;
	private double CD;
	
	DataCollection<T> dataCollection;
	WeightCalculator<T> weightCalculator;

	public Inference(boolean rooted, boolean extrarooted, List<Tree> trees,
			List<Tree> extraTrees, boolean exactSolution) {
		super();
		this.rooted = rooted;
		this.extrarooted = extrarooted;
		this.trees = trees;
		this.extraTrees = extraTrees;
		this.exactSolution = exactSolution;
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

	public void mapNames() {
		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}
		if (GlobalMaps.taxonNameMap != null && GlobalMaps.taxonNameMap.taxonMap != null) {
			Map<String,String> taxonMap = GlobalMaps.taxonNameMap.taxonMap;
			String error = Trees.checkMapping(trees, taxonMap);
			if (error != null) {
				throw new RuntimeException("Gene trees have a leaf named "
						+ error
						+ " that hasn't been defined in the mapping file");
			}

			List<String> temp1 = new LinkedList<String>();
			List<String> temp2 = new LinkedList<String>();
			for (String s : taxonMap.keySet()) {
				temp1.add(s);
				if (!( temp2).contains(taxonMap.get(s))) {
					temp2.add(taxonMap.get(s));
				}
			}
			gtTaxa = (String[]) temp2.toArray();
			/*stTaxa = new String[temp2.size()];

			
			for (int i = 0; i < stTaxa.length; i++) {
				stTaxa[i] = ((String) ((List) temp2).get(i));
				GlobalMaps.taxonIdentifier.taxonId(stTaxa[i]);
				
			}*/
		} else if (GlobalMaps.taxonNameMap != null && GlobalMaps.taxonNameMap.taxonMap == null) {
			
			Set<String> taxalist = new HashSet<String>();
			Set<String> genelist = new HashSet<String>();
			for (Tree tr : trees) {
				String[] leaves = tr.getLeaves();
				for (int i = 0; i < leaves.length; i++) {
					String leaf = leaves[i];				
					genelist.add(leaf);
					taxalist.add(GlobalMaps.taxonNameMap.getTaxonName(leaf));
				}
			}			

			//stTaxa = new String[taxalist.size()];
			gtTaxa = new String[genelist.size()];

			int index = 0;
			/*for (String taxon : taxalist) {
				stTaxa[(index++)] = taxon;
				GlobalMaps.taxonIdentifier.taxonId(taxon);
			}*/
			index = 0;
			for (String gene : genelist) {
				gtTaxa[(index++)] = gene;
			}
		} else {
/*			cd = null;
			if (rooted & extraTrees == null & GlobalMaps.taxonNameMap == null && false) {
				cd = doCollapse(trees);
			}*/

			List<String> taxalist = new ArrayList<String>();
			for (Tree tr : trees) {
				for (TNode node : tr.postTraverse()) {
					if ((node.isLeaf()) && (!taxalist.contains(node.getName()))) {
						taxalist.add(node.getName());
					}
				}
			}

			gtTaxa = new String[taxalist.size()];

			int index = 0;
			for (String taxon : taxalist) {
				gtTaxa[(index++)] = taxon;
				GlobalMaps.taxonIdentifier.taxonId(taxon);
			}
		}

		System.err.println("Number of taxa: " + GlobalMaps.taxonIdentifier.taxonCount());
		System.err.println("Taxa: " + Arrays.toString(GlobalMaps.taxonIdentifier.getAllTaxonNames()));
	}
	public abstract void scoreGeneTree(STITree scorest) ;

	List<Solution> findTreesByDP(ClusterCollection clusters) {
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
			ComputeMinCostTask allTask = newComputeMinCostTask(this,all,clusters);
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

		List minClusters = new LinkedList();
		List coals = new LinkedList();
		Stack minVertices = new Stack();
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
		while (!minVertices.isEmpty()) {
			Vertex pe = (Vertex) minVertices.pop();
			//System.out.println(pe._min_rc);
			//System.out.println(pe._min_lc);
			minClusters.add(pe.getCluster());
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
			Object tr = new STITree();
			for (String s : GlobalMaps.taxonIdentifier.getAllTaxonNames()) {
				((MutableTree) tr).getRoot().createChild(s);
			}
			sol._st = ((Tree) tr);
		} else {
			sol._st = Utils.buildTreeFromClusters(minClusters);
		}

		Object map = new HashMap();
		for (TNode node : sol._st.postTraverse()) {
			BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
			if (node.isLeaf()) {
				bs.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
				((Map) map).put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = (BitSet) ((Map) map).get(child);
					bs.or(childCluster);
				}
				((Map) map).put(node, bs);
			}
//            System.err.println("Node: "+node);
			STITreeCluster c = new STITreeCluster();
			c.setCluster(bs);
//            System.err.println("m[0]: "+((STITreeCluster)minClusters.get(0)).toString2());
//            System.err.println("C: "+c.toString2());
//            System.err.println("Equals: "+((STITreeCluster)minClusters.get(0)).equals(c));
			if (c.getClusterSize() == GlobalMaps.taxonIdentifier.taxonCount()) {
				((STINode) node).setData(Double.valueOf(0));
			} else {
				int pos = minClusters.indexOf(c);                                
				((STINode) node).setData((Double) coals.get(pos));
			}
		}

		Long cost = getTotalCost(all);
		sol._totalCoals = cost;
		System.err.println("Total Number of elements weighted: "+ weightCalculator.getCalculatedWeightCount());
		solutions.add(sol);
        System.err.println("Final optimization score: " + cost);

		return (List<Solution>) (List<Solution>) solutions;
	}
	
	public List<Solution> inferSpeciesTree() {
		long startTime = System.currentTimeMillis();

		mapNames();

		ClusterCollection clusters = newClusterCollection();

		List<Solution> solutions;

		dataCollection = newCounter(clusters);
		weightCalculator = newWeightCalculator();

		dataCollection.computeTreePartitions(this);

		if (extraTrees != null) {		
			dataCollection.addExtraBipartitionsByInput(clusters, extraTrees,extrarooted);					
		}
		
		if (exactSolution) {
			dataCollection.addAllPossibleSubClusters(clusters.getTopVertex().getCluster());
		}

		//counter.addExtraBipartitionsByHeuristics(clusters);

		System.err.println("partitions formed in "
			+ (System.currentTimeMillis() - startTime) / 1000.0D + " secs");

		weightCalculator.preCalculateWeights(trees, extraTrees);
		

		System.err.println("Dynamic Programming starting after "
				+ (System.currentTimeMillis() - startTime) / 1000.0D
				+ " secs");

		solutions = findTreesByDP(clusters);

/*		if (GlobalMaps.taxonNameMap == null && rooted && extraTrees == null && false) {
			restoreCollapse(solutions, cd);
			}*/

		return (List<Solution>) solutions;
	}

	abstract ClusterCollection newClusterCollection();
	
	abstract DataCollection newCounter(ClusterCollection clusters);
	
	abstract WeightCalculator<T> newWeightCalculator();

	abstract ComputeMinCostTask<T> newComputeMinCostTask(Inference<T> dlInference,
			Vertex all, ClusterCollection clusters);
	
	abstract Long getTotalCost(Vertex all);
	
	public double getDLbdWeigth() {
		return DLbdWeigth;
	}

	public void setDLbdWeigth(double dLbdWeigth) {
		DLbdWeigth = dLbdWeigth;
	}

	public double getCS() {
		return CS;
	}

	public void setCS(double cS) {
		CS = cS;
	}

	public double getCD() {
		return CD;
	}

	public void setCD(double cD) {
		CD = cD;
	}

	protected int [] calc(Tree gtTree, SchieberVishkinLCA lcaLookup, Tree stTree) {
		int [] res = {0,0,0};
		Stack<TNode> stack = new Stack<TNode>();			
		for (TNode gtNode : gtTree.postTraverse()) {
			if (gtNode.isLeaf()) {
			    	TNode node = stTree.getNode(GlobalMaps.taxonNameMap !=null ? 
					GlobalMaps.taxonNameMap.getTaxonName(gtNode.getName()):
						gtNode.getName());
			    	if (node == null) {
					throw new RuntimeException("Leaf " + gtNode.getName() +
						" was not found in species tree; mapped as: "+
						GlobalMaps.taxonNameMap.getTaxonName(gtNode.getName())); 
			    	}
			    	stack.push(node);
				//System.out.println("stack: " +this.taxonNameMap.getTaxonName(gtNode.getName()));
			} else {
				TNode rightLCA = stack.pop();
				TNode leftLCA = stack.pop();
				// If gene trees are incomplete, we can have this case
				if (rightLCA == null || leftLCA == null) {
					stack.push(null);
					continue;
				}
				TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
				stack.push(lca);
				if (lca == leftLCA || lca == rightLCA) {
					// LCA in stTree dominates gtNode in gene tree
					res[0]++;
					if (lca == leftLCA && lca == rightLCA) {
						res[1] += 0;
					} else {
						res[1] += (lca == leftLCA) ?
									d(rightLCA,lca) + 1:
									d(leftLCA,lca) + 1;
					}
				} else {
					res[1] += (d(rightLCA,lca) + d(leftLCA,lca));
				}
			}
		}
		TNode rootLCA = stack.pop();
		res[2] = res[1];
		res[1] += d(rootLCA,stTree.getRoot()) + (rootLCA == stTree.getRoot()?0:1);
		return res;
	}

	private int d(TNode down, TNode upp) {
		int ret = 0;
		TNode t = down;
		//System.err.println("Down: "+down+"\nUPP: "+upp);
		while (t != upp) {ret++; t=t.getParent();}
		return Math.max(ret-1,0);
	}

	
}
