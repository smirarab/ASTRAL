package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.ForkJoinPool;


import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.util.Collapse;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;

public class MGDInference_DP {
	static boolean _print = true;
	int optimizeDuploss = 1;
	boolean rooted = true;
	boolean fast = false;
	boolean extrarooted = true;

	List<Tree> trees;
	private List<Tree> extraTrees = null;
	Map<STITreeCluster, Vertex> clusterToVertex;
	Map<Integer, Set<Vertex>> clusters;
	int sigmaNs;
	DuplicationWeightCounter counter;
	TaxonNameMap taxonNameMap = null;
	
	class TaxonNameMap {
		Map<String, String> taxonMap;
		String pattern = null;
		String rep = null;
		public TaxonNameMap (Map<String, String> taxonMap) {
			this.taxonMap = taxonMap;
		}
		public TaxonNameMap (String pattern, String rep) {
			this.pattern = pattern;
			this.rep = rep;
		}
		public String getTaxonName(String geneName) {
			if (geneName == null || "".equals(geneName)) {
				throw new RuntimeException("Empty name?");
			}
			if (pattern != null) {
				String s = geneName.replaceAll(pattern,rep);
				//System.err.println("Mapped " + geneName + " to " + s);
				return s;
			} else {
				return taxonMap.get(geneName);
			}
		}		
	}

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		this.taxonNameMap = null;
	}

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees,
			Map<String, String> taxonMap) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		this.taxonNameMap = new TaxonNameMap(taxonMap);
	}

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees,
			String pattern, String rep) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		this.taxonNameMap = new TaxonNameMap (pattern, rep);
	}
	
	public static void main(String[] args) {
		if ((args == null) || args.length == 0 || (args[0].equals("-h"))
				|| (args.length < 1)) {
			printUsage();
			return;
		}
		
		boolean optimizeDuploss = false;
		boolean rooted = true;
		boolean fast = false;
		boolean extrarooted = true;
		
		Map<String, String> taxonMap = null;
		String rep = null;
		String pattern = null;
		List<Tree> trees = null;
		List<Tree> extraTrees = null;
		String output = null;
		// boolean explore = false;
		// double proportion = 0.0D;
		// boolean exhaust = false;
		double bootstrap = 1.0D;
		double time = -1.0D;
		boolean unresolved = false;
		long startTime = System.currentTimeMillis();
		String line;
		BufferedReader treeBufferReader = null;
		BufferedReader extraTreebuffer = null;
		try {
			List<String[]> options = getOptions(args);
			for (String[] option : options) {
				if (option[0].equals("-i")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					treeBufferReader = new BufferedReader(new FileReader(
							option[1]));

					trees = new ArrayList();
					// String line;
				} else if (option[0].equals("-ex")) {
					if (option.length != 2) {
						printUsage();
						return;
					}					
					extraTreebuffer = new BufferedReader(new FileReader(
							option[1]));

					extraTrees = new ArrayList();					
				} else if (option[0].equals("-a")) {
					if ( (option.length != 2) && (option.length != 3)) {
						printUsage();
						return;
					}
					if (option.length == 2) {
						BufferedReader br = new BufferedReader(new FileReader(
								option[1]));
	
						taxonMap = new HashMap<String, String>();
						while ((line = br.readLine()) != null) {
							// String line;
							String[] mapString = line.trim().split(";");
							for (String s : mapString) {
								String species = s.substring(0, s.indexOf(":"))
										.trim();
								s = s.substring(s.indexOf(":") + 1);
								String[] alleles = s.split(",");
								for (String allele : alleles) {
									allele = allele.trim();
									if (taxonMap.containsKey(allele)) {
										System.err
												.println("The input file is not in correct format");
										System.err
												.println("An gene name can only map to one species");
										System.exit(-1);
									} else {
										taxonMap.put(allele, species);
									}
								}
							}
						}
						br.close();
					} else {
						pattern = option[1];
						rep = option [2];
						if (rep.equals("/delete/")) {
							rep = "";
						}
					}
				} else if (option[0].equals("-o")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					output = option[1];
					setPrinting (false);
				} else if (option[0].equals("-x")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					// exhaust = true;
				} else if (option[0].equals("-e")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					// explore = true;
					if (option.length != 2)
						continue;
					try {
						// proportion = Double.parseDouble(option[1]);
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-b")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					if (option.length != 2)
						continue;
					try {
						bootstrap = Double.parseDouble(option[1]);
						if ((bootstrap <= 1.0D) && (bootstrap > 0.0D))
							continue;
						printUsage();
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-ur")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					unresolved = true;
				} else if (option[0].equals("-f")) {
					if ((option.length != 1)) {
						printUsage();
						return;
					}
					fast = true;
				} else if (option[0].equals("-u")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					rooted = false;
				} else if (option[0].equals("-xu")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					extrarooted = false;
				} else if (option[0].equals("-t")) {
					if ((option.length != 2) || (unresolved)) {
						printUsage();
						return;
					}
					if (option.length != 2)
						continue;
					try {
						time = Double.parseDouble(option[1]);
						if (time > 0.0D)
							continue;
						printUsage();
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}
				} else if (option[0].equals("-dl")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					optimizeDuploss = true;
				} else {
					printUsage();
					return;
				}
			}

			if (treeBufferReader == null) {
				System.err.println("The input file has not been specified.");
				printUsage();
				return;
			}

			System.err.println("Gene trees are treated as "
					+ (rooted ? "rooted" : "unrooted"));
			int l = 0;
			try {
				while ((line = treeBufferReader.readLine()) != null) {
					l++;
					Set<String> previousTreeTaxa = new HashSet<String>();
					if (line.length() > 0) {
						NewickReader nr = new NewickReader(new StringReader(line));
						if (rooted) {
							STITree gt = new STITree(true);
							nr.readTree(gt);
							if (previousTreeTaxa.isEmpty()) {
								previousTreeTaxa.addAll(Arrays.asList(gt
										.getLeaves()));
							} else {
								if (!previousTreeTaxa.containsAll(Arrays.asList(gt
										.getLeaves()))) {
									throw new RuntimeException(
											"Not all trees are on the same set of taxa: "
													+ gt.getLeaves() + "\n"
													+ previousTreeTaxa);
								}
							}
							trees.add(gt);
						} else {						
							Tree tr = nr.readTree();
							trees.add(tr);
						}
					}
				}
				treeBufferReader.close();
			} catch (ParseException e) {
				treeBufferReader.close();
				throw new RuntimeException("Failed to Parse Tree number: " + l ,e);
			}			

			if (extraTreebuffer != null) {
				while ((line = extraTreebuffer.readLine()) != null) {
					if (line.length() > 0) {					
						NewickReader nr = new NewickReader(
								new StringReader(line));
						if (extrarooted) {
							STITree gt = new STITree(true);
							nr.readTree(gt);
							extraTrees.add(gt);
						} else {
							Tree tr = nr.readTree();
							extraTrees.add(tr);
						}
					}
				}				
				extraTreebuffer.close();
			}
		} catch (IOException e) {
			System.err.println("Error when reading trees. The function exits.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		} catch (ParseException e) {
			System.err
					.println("Error when parsing the Newick representation from input file.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		}

		if (trees.size() == 0) {
			System.err.println("Empty list of trees. The function exits.");
			return;
		}

		if (_print) {
			System.err.println("Reading trees in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		startTime = System.currentTimeMillis();
		MGDInference_DP inference;
		
		
		
		if (rep != null) {			
			inference = new MGDInference_DP(trees, extraTrees,
				pattern, rep);
		} else if (taxonMap != null) {
			inference = new MGDInference_DP(trees, extraTrees,
				taxonMap);
		} else {
			inference = new MGDInference_DP(trees, extraTrees);
		}
		
		inference.optimizeDuploss = optimizeDuploss ? 3 : 1;
		inference.rooted = rooted;
		inference.fast = fast;
		inference.extrarooted = extrarooted;

		List<Solution> solutions = inference.inferSpeciesTree();

		//if (_print) {
			System.err.println("Optimal tree inferred in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		//}

		if ((_print)) {
			for (Solution s : solutions)
				System.out.println(
					        s._st.toStringWD()
						+ " \n"
						+ s._totalCoals
						+ (optimizeDuploss ? " duplication+loss"
								: " duplicatins") + " in total");
		} else
			try {
				FileWriter fw = new FileWriter(output);
				for (Solution s : solutions) {
					fw.write(s._st.toString()+ " \n");
				}
				fw.close();
			} catch (IOException e) {
				System.err.println("Error when writing the species tree");
				System.err.println(e.getMessage());
				e.printStackTrace();
			}
	}

	protected static void printUsage() {
		System.out
				.println("This tool infers the species tree from rooted gene trees despite lineage sorting.");
		System.out.println("Usage is:");
		System.out
				.println("\tMGDInference_DP -i input [-a mapping] [-dl] [-u] [-ex extra_trees] [-xu] [-o output] [-f] ");
		System.out
				.println("\t-i gene tree file: The file containing gene trees. (required)");
		System.out
				.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled. Or, you can specify two reqular expressions for automatic name conversion (optional)");
		System.out
				.println("\t-o species tree file: The file to store the species tree. (optional)");
		System.out.println("\t-dl optimize duploss instead of duplications");
		System.out.println("\t-u treat input gene trees as unrooted");
		System.out.println("\t-ex provide extra trees to add to set of STBs searched");
		System.out.println("\t-xu treat extra trees input gene trees as unrooted");
		System.out.println("\t-f perform fast and less-accurate subtree-bipartition based search.");
		System.out.println();
	}

	public static void setPrinting(boolean print) {
		_print = print;
	}

	public List<Solution> inferSpeciesTree() {
		long startTime = System.currentTimeMillis();

		String[] gtTaxa;
		String[] stTaxa;
		Collapse.CollapseDescriptor cd = null;

		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}
		if (taxonNameMap != null && taxonNameMap.taxonMap != null) {
			Map<String,String> taxonMap = taxonNameMap.taxonMap;
			String error = Trees.checkMapping(trees, taxonMap);
			if (error != null) {
				throw new RuntimeException("Gene trees have leaf named "
						+ error
						+ "that hasn't been defined in the mapping file");
			}

			List temp1 = new LinkedList();
			List temp2 = new LinkedList();
			for (String s : taxonMap.keySet()) {
				temp1.add(s);
				if (!((List) temp2).contains(taxonMap.get(s))) {
					((List) temp2).add((String) taxonMap.get(s));
				}
			}
			gtTaxa = new String[temp1.size()];
			stTaxa = new String[temp2.size()];

			for (int i = 0; i < gtTaxa.length; i++) {
				gtTaxa[i] = ((String) temp1.get(i));
			}
			for (int i = 0; i < stTaxa.length; i++) {
				stTaxa[i] = ((String) ((List) temp2).get(i));
			}
		} else if (taxonNameMap != null && taxonNameMap.taxonMap == null) {
			
			Set<String> taxalist = new HashSet<String>();
			Set<String> genelist = new HashSet<String>();
			for (Tree tr : trees) {
				String[] leaves = tr.getLeaves();
				for (int i = 0; i < leaves.length; i++) {
					String leaf = leaves[i];				
					genelist.add(leaf);
					taxalist.add(taxonNameMap.getTaxonName(leaf));
				}
			}			

			stTaxa = new String[taxalist.size()];
			gtTaxa = new String[genelist.size()];

			int index = 0;
			for (String taxon : taxalist) {
				stTaxa[(index++)] = taxon;
			}
			index = 0;
			for (String gene : genelist) {
				gtTaxa[(index++)] = gene;
			}
		} else {
			cd = null;
			if (rooted & extraTrees == null & taxonNameMap == null && false) {
				cd = doCollapse(trees);
			}

			List<String> taxalist = new ArrayList<String>();
			for (Tree tr : trees) {
				for (TNode node : tr.postTraverse()) {
					if ((node.isLeaf()) && (!taxalist.contains(node.getName()))) {
						taxalist.add(node.getName());
					}
				}
			}

			stTaxa = new String[taxalist.size()];

			int index = 0;
			for (String taxon : taxalist) {
				stTaxa[(index++)] = taxon;
			}
			gtTaxa = stTaxa;
		}

		System.err.println("Number of taxa: " + stTaxa.length);
		System.err.println("Taxa: " + Arrays.toString(stTaxa));

		clusters = new HashMap<Integer, Set<Vertex>>(stTaxa.length);

		List<Solution> solutions;

		counter = new DuplicationWeightCounter(gtTaxa, stTaxa, rooted,taxonNameMap, clusters);

		int sigmaN = counter.computeTreeSTBipartitions(trees);

		if (extraTrees != null) {		
			counter.addExtraBipartitionsByInput(clusters, extraTrees,extrarooted);					
		}

		counter.addExtraBipartitionsByHeuristics(clusters);

		if (_print) {
			System.err.println("STBs formed in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		counter.preCalculateWeights(trees, extraTrees);
		
		if (_print) {
			System.err.println("Weights pre-caluclated after "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		sigmaNs = (optimizeDuploss == 3) ? 
				sigmaN + 2 * (stTaxa.length - 1) * trees.size(): 
				sigmaN;

		solutions = findTreesByDP(stTaxa, counter, trees, taxonNameMap);

		if (taxonNameMap == null && rooted && extraTrees == null && false) {
			restoreCollapse(solutions, cd);
		}

		return (List<Solution>) solutions;
	}

	private Collapse.CollapseDescriptor doCollapse(List<Tree> trees) {
		Collapse.CollapseDescriptor cd = Collapse.collapse(trees);
		return cd;
	}

	private void restoreCollapse(List<Solution> sols,
			Collapse.CollapseDescriptor cd) {
		for (Solution sol : sols) {
			Tree tr = sol._st;
			Collapse.expand(cd, (MutableTree) tr);
			for (TNode node : tr.postTraverse())
				if (((STINode) node).getData() == null)
					((STINode) node).setData(Integer.valueOf(0));
		}
	}

	  public static Tree buildTreeFromClusters(List<STITreeCluster> clusters)
	  {
	    if ((clusters == null) || (clusters.size() == 0)) {
	      System.err.println("Empty list of clusters. The function returns a null tree.");
	      return null;
	    }

	    MutableTree tree = new STITree();

	    String[] taxa = ((STITreeCluster)clusters.get(0)).getTaxa();
	    for (int i = 0; i < taxa.length; i++) {
	      tree.getRoot().createChild(taxa[i]);
	    }

	    for (STITreeCluster tc : clusters) {
	      if ((tc.getClusterSize() <= 1) || (tc.getClusterSize() == tc.getTaxa().length))
	      {
	        continue;
	      }

	      Set clusterLeaves = new HashSet();
	      TNode node;
	      for (String l : tc.getClusterLeaves()) {
	        node = tree.getNode(l);
	        clusterLeaves.add(node);
	      }

	      SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(tree);
	      TNode lca = lcaFinder.getLCA(clusterLeaves);

	      Object movedChildren = new LinkedList();
	      for (TNode child : lca.getChildren()) {
	        BitSet childCluster = new BitSet(taxa.length);
	        for (TNode cl : child.getLeaves()) {
	          for (int i = 0; i < taxa.length; i++) {
	            if (taxa[i].equals(cl.getName())) {
	              childCluster.set(i);
	              break;
	            }
	          }
	        }

	        BitSet temp = (BitSet)childCluster.clone();
	        temp.and(tc.getCluster());
	        if (temp.equals(childCluster)) {
	          ((List)movedChildren).add(child);
	        }

	      }

	      STINode newChild = ((STINode)lca).createChild();

	      while (!((List)movedChildren).isEmpty()) {
	        newChild.adoptChild((TMutableNode)((List)movedChildren).get(0));
	        ((List)movedChildren).remove(0);
	      }
	    }

	    return (Tree)tree;
	  }
	  
	private List<Solution> findTreesByDP(String[] stTaxa,
			DuplicationWeightCounter counter, List<Tree> trees,
			TaxonNameMap taxonNameMap) {
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
		clusterToVertex = new HashMap<STITreeCluster, Vertex>(sigmaNs);
		for (Set<Vertex> vs : clusters.values()) {
			for (Vertex vertex : vs) {
				vertex._max_score = -1;
				clusterToVertex.put(vertex._cluster, vertex);
			}
		}

		Vertex all = (Vertex) clusters.get(Integer.valueOf(stTaxa.length))
				.toArray()[0];
		System.err.println("Sigma N: " + sigmaNs);

		System.err.println("Size of largest cluster: " +all._cluster.getClusterSize());

		try {
			//vertexStack.push(all);
			ComputeMinCostTask allTask = new ComputeMinCostTask(this,all);
			ForkJoinPool pool = new ForkJoinPool(1);
			pool.invoke(allTask);
			Integer v = all._max_score;
			if (v == -2) {
				throw new CannotResolveException(all._cluster.toString());
			}
		} catch (CannotResolveException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new RuntimeException("Was not able to build a fully resolved tree ",e);
		}

		if (_print) {
			//System.err.println("Weights are: "
				//	+ counter.weights);
		}
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
			minClusters.add(pe._cluster);
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
				coals.add(0);
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
			for (String s : stTaxa) {
				((MutableTree) tr).getRoot().createChild(s);
			}
			sol._st = ((Tree) tr);
		} else {
			sol._st = buildTreeFromClusters(minClusters);
		}
		//System.err.println("SOL: " + sol._st);
		//System.err.println("coals: " + coals);
		//System.err.println("min cluster: " + minClusters);
		Object map = new HashMap();
		for (TNode node : sol._st.postTraverse()) {
			BitSet bs = new BitSet(stTaxa.length);
			if (node.isLeaf()) {
				for (int i = 0; i < stTaxa.length; i++) {
					if (node.getName().equals(stTaxa[i])) {
						bs.set(i);
						break;
					}
				}
				((Map) map).put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = (BitSet) ((Map) map).get(child);
					bs.or(childCluster);
				}
				((Map) map).put(node, bs);
			}
//            System.err.println("Node: "+node);
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.setCluster(bs);
//            System.err.println("m[0]: "+((STITreeCluster)minClusters.get(0)).toString2());
//            System.err.println("C: "+c.toString2());
//            System.err.println("Equals: "+((STITreeCluster)minClusters.get(0)).equals(c));
			if (c.getClusterSize() == stTaxa.length) {
				((STINode) node).setData(Integer.valueOf(0));
			} else {
				int pos = minClusters.indexOf(c);                                
				((STINode) node).setData((Integer) coals.get(pos));
			}
		}

		sol._totalCoals = all._min_cost;
		System.out.println("total cost: " + all._min_cost);
		solutions.add(sol);

		return (List<Solution>) (List<Solution>) solutions;
	}

	// TODO: fix this.	
	int maxEL = 1000000000;
	
/*	private int computeMinCost(Vertex v) throws CannotResolveException {
		
		if (v._max_score == -2) {
			throw new CannotResolveException(v._cluster.toString());
		}
		// Already calculated. Don't re-calculate.
		if (v._max_score != -1) {
			return v._max_score - maxEL;
		}
		// If in duploss mode, need to get MDC cost as well
		if (optimizeDuploss) {
			if (v._el_num == -1) {
				if (taxonNameMap == null) {
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, rooted);
				} else {
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, taxonNameMap, rooted);
				}
			}
		} else {
			v._el_num = 0;
		}
		
		int clusterSize = v._cluster.getClusterSize();
		
		// SIA: base case for singelton clusters.
		{
			if (clusterSize <= 1) {
				// SIA: TODO: this is 0, right?
				v._min_cost = 0;
				v._max_score = maxEL - v._el_num;
				v._min_lc = (v._min_rc = null);
				return v._max_score - maxEL;
			}
		}
		Set<STBipartition> clusterBiPartitions = counter
				.getClusterBiPartitions(v._cluster);

		// STBipartition bestSTB = null;
	 	if (fast) {	
			if (clusterBiPartitions == null) {
				if (v._cluster.getClusterSize() <= 3) {
					for (int j=0; j< v._cluster.getClusterSize(); j++){
						STITreeCluster c1 = new STITreeCluster(v._cluster.getTaxa());
						c1.addLeaf(v._cluster.getClusterLeaves()[j]);
						STITreeCluster c2 = new STITreeCluster(v._cluster.getTaxa());
						for (int i = 0; i < v._cluster.getClusterSize(); i++) {
							if (i != j){
								c2.addLeaf(v._cluster.getClusterLeaves()[i]);
							}
						}
						if (clusterToVertex.containsKey(c2)) {
							STBipartition stb = new STBipartition(c1, c2, v._cluster);
							if (clusterBiPartitions == null)
								clusterBiPartitions = new HashSet<STBipartition>(3);
							clusterBiPartitions.add(stb);
							System.err.println("Adding: " + stb);
						}
					}
				}
			}
	
			if (clusterBiPartitions == null) {
				System.err.println("Warn: the following cluster ( " + v._cluster.getClusterSize()+" taxa ) has no STBs:\n"
						+ v._cluster);
				v._max_score = -2;
				throw new CannotResolveException(v._cluster.toString());
			}
			for (STBipartition stb : clusterBiPartitions) {
	
				Vertex lv = clusterToVertex.get(stb.cluster1);
				Vertex rv = clusterToVertex.get(stb.cluster2);
	
				
				 * if (lv == null || rv == null) {
				 * //System.out.println("There is no STB for one half of : " + stb);
				 * continue; }
				 
	
				try {
	
					//vertexStack.push(lv);
					int lscore = computeMinCost(lv);
					//vertexStack.pop();
					//vertexStack.push(rv);
					int rscore = computeMinCost(rv);
					//vertexStack.pop();
	
					int w = counter.getBiPartitionDPWeight(lv._cluster,
							rv._cluster, v._cluster);
	
					int z = optimizeDuploss ? 3 : 1;
	
					int c = z * w - v._el_num;
	
					if ((v._max_score != -1)
							&& (lscore + rscore + c + maxEL <= v._max_score)) {
						continue;
					}
					v._max_score = (lscore + rscore + c) + maxEL;
					v._min_cost = sigmaNs
							- (c + lv._max_score + rv._max_score - 2 * maxEL);
					// stem.out.println(maxEL - (z*w + lv._max_score +
					// rv._max_score));
					v._min_lc = lv;
					v._min_rc = rv;
					v._c = c;				
				} catch (CannotResolveException c) {
					System.err.println("Warn: cannot resolve: " + c.getMessage());
				}
	
				// bestSTB = stb;
			}
		} else {

			for (int i = 1; i <= (clusterSize / 2); i++) {
				Set<Vertex> leftList = clusters.get(i);
				if (leftList != null) {
					for (Vertex lv : leftList) {
						if (v._cluster.containsCluster(lv._cluster)) {						
	
							Set<Vertex> rightList = clusters.get(clusterSize - i);
							if (rightList != null) {							
								for (Vertex rv : rightList) {
									if (! (lv._cluster.isDisjoint(rv._cluster) && v._cluster.containsCluster(rv._cluster)) ) {
										continue;
									}
									try {
										int lscore = computeMinCost(lv);
										int rscore = computeMinCost(rv);
	
										int w = counter.getBiPartitionDPWeight(lv._cluster,
												rv._cluster, v._cluster);
	
										int z = optimizeDuploss ? 3 : 1;
	
										int c = z * w - v._el_num;
	
										if ((v._max_score != -1)
												&& (lscore + rscore + c + maxEL <= v._max_score)) {
											continue;
										}
										v._max_score = (lscore + rscore + c) + maxEL;
										v._min_cost = sigmaNs
												- (c + lv._max_score + rv._max_score - 2 * maxEL);
										// stem.out.println(maxEL - (z*w + lv._max_score +
										// rv._max_score));
										v._min_lc = lv;
										v._min_rc = rv;
										v._c = c;
	
										break;	// Already found the only pair of clusters whose union is v's cluster.
									} catch (CannotResolveException c) {
										//System.err.println("Warn: cannot resolve: " + c.getMessage());											
									}
								}								
							}
						}
					}
				}
			}
	
		}
		
		if (v._min_lc == null) {
			System.err.println("WARN: No Resolution found for ( " + v._cluster.getClusterSize()+" taxa ):\n"
					+ v._cluster);
			v._max_score = -2;
			throw new CannotResolveException(v._cluster.toString());
		}

		
		 * if (clusterSize > 5){ counter.addGoodSTB(bestSTB, clusterSize); }
		 
		return v._max_score - maxEL;
	}
*/
	/*
	 * private int computeAllClusters(List<Tree> trees, String[] stTaxa,
	 * Map<String, String> taxonMap, Map<Integer, List<Vertex>> clusters) { int
	 * n = stTaxa.length; int maxEL = 0; if (n <= 0) {
	 * System.err.println("Empty list of taxa."); return -1; }
	 * 
	 * BitSet counter = new BitSet(n); boolean done = false;
	 * 
	 * while (!done) { int i = 0; while ((i < n) && (counter.get(i))) {
	 * counter.clear(i); i++; } if (i >= n) { done = true; } else {
	 * counter.set(i, true);
	 * 
	 * STITreeCluster tc = new STITreeCluster(stTaxa); tc.setCluster((BitSet)
	 * counter.clone()); Vertex v = new Vertex(); v._cluster = tc; if (taxonMap
	 * == null) { v._el_num = DuplicationWeightCounter.getClusterCoalNum( trees,
	 * tc, true); } else { v._el_num =
	 * DuplicationWeightCounter.getClusterCoalNum( trees, tc, taxonMap, true); }
	 * if (v._el_num > maxEL) { maxEL = v._el_num; }
	 * 
	 * int size = tc.getClusterSize(); List l = (List)
	 * clusters.get(Integer.valueOf(size)); if (l == null) { l = new
	 * LinkedList(); }
	 * 
	 * l.add(v); clusters.put(Integer.valueOf(size), l); } } maxEL++; return
	 * maxEL; }
	 */

	private int computeAllClusters(List<Tree> trees, String[] stTaxa,
			Map<Integer, List<Vertex>> clusters) {
		throw new RuntimeException("Not Implemented");
		// return computeAllClusters(trees, stTaxa, null, clusters);
	}

	protected static List<String[]> getOptions(String[] args) {
		LinkedList opts = new LinkedList();
		LinkedList arg_list = new LinkedList();

		int i = 0;
		while (i < args.length) {
			if (args[i].charAt(0) != '-') {
				printUsage();
				System.exit(-1);
			}

			arg_list.clear();
			arg_list.addFirst(args[i]);
			i++;

			while ((i < args.length) && (args[i].charAt(0) != '-')) {
				arg_list.addLast(args[i]);
				i++;
			}
			String[] arg_array = new String[arg_list.size()];
			arg_list.toArray(arg_array);

			opts.addLast(arg_array);
		}
		return opts;
	}

	private int getResolutionsNumber(int nodeNumber) {
		int total = 1;
		for (int i = 3; i <= nodeNumber; i++) {
			total *= (2 * i - 3);
		}
		return total;
	}

	/*
	 * private int tryBinaryResolutions(Tree tr, double time, String[] taxa,
	 * List<Tree> gts, Map<String, String> taxonMap) { List nodelist = new
	 * ArrayList(); List degreelist = new ArrayList(); int totalResolutions = 0;
	 * for (Iterator iterator = (new PostTraversal(tr.getRoot())).iterator();
	 * iterator .hasNext();) { TNode node = (TNode) iterator.next(); int
	 * childCount = node.getChildCount(); if (childCount > 2) {
	 * nodelist.add(node); int resolutionsNumber =
	 * getResolutionsNumber(childCount);
	 * degreelist.add(Integer.valueOf(resolutionsNumber)); totalResolutions +=
	 * resolutionsNumber; } } int addedxl = 0; for (int i = 0; i <
	 * nodelist.size(); i++) { TNode unresolvedNode = (TNode) nodelist.get(i);
	 * Map id2node = new HashMap(); for (TNode child :
	 * unresolvedNode.getChildren()) {
	 * id2node.put(Integer.valueOf(child.getID()), child); } Integer[] childIDs
	 * = (Integer[]) id2node.keySet().toArray( new Integer[0]); Object
	 * cluster2xl = new HashMap(); double endtime; if (time == -1.0D) { endtime
	 * = -1.0D; } else { endtime = time (((Integer)
	 * degreelist.get(i)).intValue() / totalResolutions) 1000.0D +
	 * System.currentTimeMillis(); } Solution sol = addMoreLeaves(null,
	 * childIDs, 0, id2node, taxa, gts, taxonMap, endtime, (Map) cluster2xl);
	 * TNode parent = unresolvedNode.getParent(); int xl = ((Integer) ((STINode)
	 * unresolvedNode).getData()) .intValue(); ((STINode)
	 * unresolvedNode).removeNode(); if (parent != null) { TNode newnode =
	 * ((STINode) parent).createChild(sol.getTree() .getRoot()); ((STINode)
	 * newnode).setData(Integer.valueOf(xl)); } else { for (TNode child :
	 * sol.getTree().getRoot().getChildren()) { ((STINode)
	 * tr.getRoot()).createChild(child); } } addedxl += sol.getCoalNum(); }
	 * return addedxl; }
	 */

	/*
	 * private Solution addMoreLeaves(STITree<Integer> preTree, Integer[]
	 * leavesid, int index, Map<Integer, TNode> id2node, String[] taxa,
	 * List<Tree> gts, Map<String, String> taxonMap, double endTime, Map<BitSet,
	 * Integer> cluster2xl) { if (preTree == null) { preTree = new
	 * STITree(false); STINode root = preTree.getRoot(); STINode newnode =
	 * root.createChild(); newnode.setData(leavesid[(index++)]); STINode innode
	 * = root.createChild(); newnode = innode.createChild();
	 * newnode.setData(leavesid[(index++)]); newnode = innode.createChild();
	 * newnode.setData(leavesid[(index++)]); } Solution sol = null; if (index ==
	 * leavesid.length) { sol = tryAllRootings(preTree, id2node, taxa, gts,
	 * taxonMap, endTime, cluster2xl); } else { int id =
	 * leavesid[(index++)].intValue(); for (Iterator iterator = (new
	 * PostTraversal(preTree.getRoot())) .iterator(); iterator.hasNext();) {
	 * TNode n = (TNode) iterator.next(); if (!n.isLeaf()) { if
	 * (n.getChildCount() != 2) { throw new RuntimeException("Not binary!"); }
	 * Iterator it = n.getChildren().iterator(); for (int i = 0; i < 2; i++) {
	 * STINode child = (STINode) it.next(); STITree newTree = new
	 * STITree(preTree); TNode peerChild = newTree.getNode(child.getID()); TNode
	 * peerParent = peerChild.getParent(); STINode newchild = ((STINode)
	 * peerParent).createChild(); newchild.adoptChild((TMutableNode) peerChild);
	 * STINode newnode = newchild.createChild();
	 * newnode.setData(Integer.valueOf(id)); Solution thissol =
	 * addMoreLeaves(newTree, leavesid, index, id2node, taxa, gts, taxonMap,
	 * endTime, cluster2xl); if ((sol == null) || (sol.getCoalNum() >
	 * thissol.getCoalNum())) { sol = thissol; } double now =
	 * System.currentTimeMillis(); if ((endTime != -1.0D) && (now > endTime)) {
	 * return sol; } if (n.isRoot()) { break; } } } } } return sol; }
	 */

	/*
	 * private Solution tryAllRootings(Tree subtree, Map<Integer, TNode>
	 * id2node, String[] taxa, List<Tree> gts, Map<String, String> taxonMap,
	 * double endTime, Map<BitSet, Integer> cluster2xl) { Solution sol = new
	 * Solution(); sol._totalCoals = -1; for (Tree rootedsubtree :
	 * subtree.getAllRootingTrees()) { TNode peerNode; for (Iterator iterator1 =
	 * (new PostTraversal(rootedsubtree .getRoot())).iterator();
	 * iterator1.hasNext();) { TNode replacingNode = (TNode) iterator1.next();
	 * if (replacingNode.isLeaf()) { int id = ((Integer) ((STINode)
	 * replacingNode).getData()) .intValue(); peerNode = (TNode)
	 * id2node.get(Integer.valueOf(id)); if (peerNode.isLeaf()) { ((STINode)
	 * replacingNode).setName(peerNode.getName()); } else { for (TNode child :
	 * peerNode.getChildren()) { ((STINode) replacingNode).createChild(child); }
	 * ((STINode) replacingNode).setName(""); }
	 * 
	 * ((STINode) replacingNode) .setData((Integer) ((STINode)
	 * peerNode).getData()); }
	 * 
	 * }
	 * 
	 * int xl = 0; Object map = new HashMap(); for (Iterator iterator2 = (new
	 * PostTraversal(rootedsubtree .getRoot())).iterator();
	 * iterator2.hasNext();) { TNode node = (TNode) iterator2.next(); BitSet bs
	 * = new BitSet(); if (node.isLeaf()) { for (int i = 0; i < taxa.length;
	 * i++) { if (node.getName().equals(taxa[i])) { bs.set(i); break; } } ((Map)
	 * map).put(node, bs); } else { for (TNode child : node.getChildren()) {
	 * BitSet childCluster = (BitSet) ((Map) map).get(child);
	 * bs.or(childCluster); } ((Map) map).put(node, bs); } if (!node.isRoot()) {
	 * Integer el = (Integer) ((STINode) node).getData(); if (el == null) { el =
	 * (Integer) cluster2xl.get(bs); if (el == null) { STITreeCluster tc = new
	 * STITreeCluster(taxa); tc.setCluster(bs); if (taxonMap == null) { el =
	 * Integer.valueOf(DuplicationWeightCounter .getClusterCoalNum(gts, tc,
	 * true)); } else { el = Integer.valueOf(DuplicationWeightCounter
	 * .getClusterCoalNum(gts, tc, taxonMap, true)); } cluster2xl.put(bs, el); }
	 * ((STINode) node).setData(el); xl += el.intValue(); } } }
	 * 
	 * if ((sol.getCoalNum() == -1) || (sol.getCoalNum() > xl)) {
	 * sol._totalCoals = xl; sol._st = rootedsubtree; } } return (Solution) sol;
	 * }
	 */

	static class Vertex {
		public STITreeCluster _cluster = null;
		public int _el_num = -1;
		public int _min_cost = -1;
		public int _max_score = -1;
		public int _c = 0;
		public Vertex _min_lc = this._min_rc = null;
		public Vertex _min_rc;
		public List<Vertex> _subcl = null;

		public Vertex() {
		}

		public String toString() {
			return this._cluster.toString() + "/" + this._max_score;
		}

		@Override
		public boolean equals(Object obj) {
			return ((Vertex) obj)._cluster.equals(this._cluster);
		}

		@Override
		public int hashCode() {
			return _cluster.hashCode();
		}

	}
}
