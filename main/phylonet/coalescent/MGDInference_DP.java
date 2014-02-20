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
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.util.Collapse;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;
import phylonet.coalescent.GlobalMaps.*;

public class MGDInference_DP {
	static boolean _print = true;
	int optimizeDuploss = 1; //one means dup, 3 means duploss
	boolean rooted = true;
	boolean fast = false;
	boolean extrarooted = true;
	double DLbdWeigth;
	double CS;
	double CD;

	List<Tree> trees;
	private List<Tree> extraTrees = null;
	//Map<STITreeCluster, Vertex> clusterToVertex;
	double sigmaNs;
	DuplicationWeightCounter counter;
	private boolean exactSolution;
	private STITree st = null;
	

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees,
			Map<String, String> taxonMap) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		if (taxonMap != null) {
			GlobalMaps.taxonNameMap = new TaxonNameMap(taxonMap);
		}
	}

	public MGDInference_DP(List<Tree> trees, List<Tree> extraTrees,
			String pattern, String rep) {
		super();
		this.trees = trees;
		this.extraTrees = extraTrees;
		GlobalMaps.taxonNameMap = new TaxonNameMap (pattern, rep);
	}
	
	public static void main(String[] args) {
		if ((args == null) || args.length == 0 || (args[0].equals("-h"))
				|| (args.length < 1)) {
			printUsage();
			return;
		}
		
		int optimizeDuploss = 0;
		boolean rooted = true;
		boolean fast = false;
		boolean extrarooted = true;
		boolean exactSolution = false;
		
		Map<String, String> taxonMap = null;
		String rep = null;
		String pattern = null;
		List<Tree> trees = null;
		List<Tree> extraTrees = null;
		String output = null;
		STITree scorest = null;
		// boolean explore = false;
		// double proportion = 0.0D;
		// boolean exhaust = false;
		double bootstrap = 1.0D;
		double cs = 1.0D;
		double cd = 1.0D;
		double time = -1.0D;
		double wd = 1.0D;
		double wh = 1.0D;
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
				} else if (option[0].equals("-st")) {
					if (option.length != 2) {
						printUsage();
						return;
					}					
					BufferedReader tmp = new BufferedReader(new FileReader(
							option[1]));
					line = tmp.readLine();
					NewickReader nr = new NewickReader(new StringReader(line));
					scorest = new STITree(true);
					nr.readTree(scorest);				
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
												.println("Any gene name can only map to one species");
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

				} else if (option[0].equals("-cs")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					try {
						cs = Double.parseDouble(option[1]);
						if ((cs <= 1.0D) && (cs >= 0.0D))
							continue;
						printUsage();
						return;
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter");
						printUsage();
						return;
					}

				} else if (option[0].equals("-cd")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					try {
						cd = Double.parseDouble(option[1]);
						if ((cd <= 1.0D) && (cd >= 0.0D))
							continue;
						printUsage();
						return;
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
				} else if (option[0].equals("-xt")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					exactSolution = true;
				} else if (option[0].equals("-dl")) {
					optimizeDuploss = 1;
					if (option.length != 2) {
						printUsage();
						return;
					}					
					try {
						if (option[1].equals("auto")) {
							wh = -1;
							continue;
						} else {
							wh = Double.parseDouble(option[1]);
							if (wh >= 0.0D)
								continue;
						}
						printUsage();
						return;
					} catch (NumberFormatException e) {
						System.err.println("Error in reading parameter wd");
						printUsage();
						return;
					}
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
			inference = new MGDInference_DP(trees, extraTrees, null);
		}
		
		inference.optimizeDuploss = optimizeDuploss > 0 ? 3 : 1;
		inference.DLbdWeigth = wh; 
		inference.rooted = rooted;
		inference.fast = fast;
		inference.extrarooted = extrarooted;
		inference.CS = cs;
		inference.CD = cd;
		inference.exactSolution = exactSolution;
		inference.st = scorest;

		if (scorest != null) {
			inference.scoreGeneTree();
			System.exit(0);
		}
		
		List<Solution> solutions = inference.inferSpeciesTree();

		//if (_print) {
			System.err.println("Optimal tree inferred in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		//}

		if ((_print)) {
			String metric;
			if (optimizeDuploss == 0) {
				metric = "duplications";
			} else if (optimizeDuploss == 1) {
				metric = "duplication+loss (homomorphic)";
			} else {
				metric = "duplication+loss (original)";
			}
			for (Solution s : solutions)
				System.out.println(
					        s._st.toStringWD()
						+ " \n"
						+ s._totalCoals
						+ " " + metric + " in total");
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

	private int [] calc(Tree gtTree, SchieberVishkinLCA lcaLookup, Tree stTree) {
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
	
	private void scoreGeneTree() {
		// first calculated duplication cost by looking at gene trees. 
		
		SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(this.st);
		Integer duplications = 0;
		Integer losses = 0;
		Integer lossesstd = 0;
		
		for (Tree gtTree : this.trees) {
			int[] res = calc(gtTree,lcaLookup, this.st);
			duplications += res[0];
			losses += res[1];
			
			STITree hmst = new STITree(this.st);
			//hmst.constrainByLeaves(Arrays.asList(gtTree.getLeaves()));
			SchieberVishkinLCA hmlcaLookup = new SchieberVishkinLCA(hmst);
			int[] res2 = calc(gtTree,hmlcaLookup, hmst);
			
			lossesstd += res2[2];
		}
		System.out.println("Total number of duplications is: "+duplications);
		System.out.println("Total number of losses (bd) is: "+losses);
		System.out.println("Total number of losses (std) is: "+lossesstd);
		System.out.println("Total number of duploss (bd) is: " + (losses+duplications));
		System.out.println("Total number of duploss (st) is: " + (lossesstd+duplications));
		System.out.println("Total weighted (wd = "+this.DLbdWeigth+") loss is: " + (lossesstd + this.DLbdWeigth*(losses-lossesstd)));
	}

	private int d (TNode down, TNode upp) {
		int ret = 0;
		TNode t = down;
		//System.err.println("Down: "+down+"\nUPP: "+upp);
		while (t != upp) {ret++; t=t.getParent();}
		return Math.max(ret-1,0);
	}
	protected static void printUsage() {
		System.out
				.println("This tool infers the species tree from rooted gene trees despite lineage sorting.");
		System.out.println("Usage is:");
		System.out
				.println("\tMGDInference_DP -i input [-a mapping] [-dl | -dll] [-ex extra_trees] [-o output] [-cs number] [-cd number] [-xt] [-s species tree] [-wd duplication weight]");
		System.out
				.println("\t-i gene tree file: The file containing gene trees. (required)");
		System.out
				.println("\t-st species tree file: The file containing a species tree to be scored.\n" +
						 "\t                       If this option is provided the software only scores the species tree.");
		System.out
				.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled.\n" +
						 "\t                 Alternatively, two reqular expressions for automatic name conversion (optional)");
		System.out
				.println("\t-o species tree file: The file to store the species tree. (optional)");
		System.out.println("\t-dl N: optimize duplications and losses. Use -dl 0 for standard (homomorphic) definition, and -dl 1 for ``bd'' definition. Any value in between weights the impact of missing taxa on the tree.");
		System.out.println("\t-xt find the exact solution by looking at all clusters.");
		//System.out.println("\t-u treat input gene trees as unrooted (Not implemented!)");
		System.out.println("\t-ex provide extra trees to add to set of STBs searched");
		//System.out.println("\t-xu treat extra trees input gene trees as unrooted (Not implemented!)");
		System.out.println("\t-cs and " +
						   "-cd thes two options set two parameters (cs and cd) to a value between 0 and 1. \n" +
						   "\t    For any cluster C if |C| >= cs*|taxa|, we add complementary clusters (with respect to C) of all subclusters of C\n" +
						   "\t    if size of the subcluster is >= cd*|C|.\n" +
						   "\t    By default cs = cd = 1; so no extra clusters are added. Lower cs and cd values could result in better scores\n" +
						   "\t    (especially when gene trees have low taxon occupancy) but can also increase the running time dramatically.");
		
		//System.out.println("\t-f perform fast and less-accurate subtree-bipartition based search (Not implemented!).");
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
		if (GlobalMaps.taxonNameMap != null && GlobalMaps.taxonNameMap.taxonMap != null) {
			Map<String,String> taxonMap = GlobalMaps.taxonNameMap.taxonMap;
			String error = Trees.checkMapping(trees, taxonMap);
			if (error != null) {
				throw new RuntimeException("Gene trees have a leaf named "
						+ error
						+ " that hasn't been defined in the mapping file");
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
				GlobalMaps.taxonIdentifier.taxonId(stTaxa[i]);
				
			}
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

			stTaxa = new String[taxalist.size()];
			gtTaxa = new String[genelist.size()];

			int index = 0;
			for (String taxon : taxalist) {
				stTaxa[(index++)] = taxon;
				GlobalMaps.taxonIdentifier.taxonId(taxon);
			}
			index = 0;
			for (String gene : genelist) {
				gtTaxa[(index++)] = gene;
			}
		} else {
			cd = null;
			if (rooted & extraTrees == null & GlobalMaps.taxonNameMap == null && false) {
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
				GlobalMaps.taxonIdentifier.taxonId(taxon);
			}
			gtTaxa = stTaxa;
		}

		System.err.println("Number of taxa: " + stTaxa.length);
		System.err.println("Taxa: " + Arrays.toString(stTaxa));

		 ClusterCollection clusters = new BasicClusterCollection(stTaxa.length);

		List<Solution> solutions;

		counter = new DuplicationWeightCounter(gtTaxa, stTaxa, rooted, clusters);

		double sigmaN = counter.computeTreeSTBipartitions(this);

		if (extraTrees != null) {		
			counter.addExtraBipartitionsByInput(clusters, extraTrees,extrarooted);					
		}
		
		if (exactSolution) {
			counter.addAllPossibleSubClusters(clusters.getTopVertex().getCluster(),  stTaxa.length);
		}

		//counter.addExtraBipartitionsByHeuristics(clusters);

		if (_print) {
			System.err.println("STBs formed in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		counter.preCalculateWeights(trees, extraTrees);
		
		if (_print) {
			System.err.println("DP starting after "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		sigmaNs = sigmaN ;

		solutions = findTreesByDP(stTaxa, counter, trees, GlobalMaps.taxonNameMap,clusters);

		if (GlobalMaps.taxonNameMap == null && rooted && extraTrees == null && false) {
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

	    //String[] taxa = ((STITreeCluster)clusters.get(0)).getTaxa();
	    for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++) {
	      tree.getRoot().createChild(GlobalMaps.taxonIdentifier.getTaxonName(i));
	    }

	    for (STITreeCluster tc : clusters) {
	      if ((tc.getClusterSize() <= 1) || (tc.getClusterSize() == GlobalMaps.taxonIdentifier.taxonCount()))
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
	        BitSet childCluster = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
	        for (TNode cl : child.getLeaves()) {
	          int i = GlobalMaps.taxonIdentifier.taxonId(cl.getName());
	          childCluster.set(i);
	        }
	        

	        BitSet temp = (BitSet)childCluster.clone();
	        temp.and(tc.getBitSet());
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
			TaxonNameMap taxonNameMap, ClusterCollection clusters) {
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
		System.err.println("Sigma N: " + sigmaNs);

		System.err.println("Size of largest cluster: " +all.getCluster().getClusterSize());

		try {
			//vertexStack.push(all);
			ComputeMinCostTask allTask = new ComputeMinCostTask(this,all,clusters);
			//ForkJoinPool pool = new ForkJoinPool(1);
			allTask.compute();
			double v = all._max_score;
			if (v == Integer.MIN_VALUE) {
				throw new CannotResolveException(all.getCluster().toString());
			}
		} catch (CannotResolveException e) {
			System.err.println("Was not able to build a fully resolved tree. Not" +
					"enough STBs present in input gene trees ");
			e.printStackTrace();
			System.exit(1);
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
			minClusters.add(pe.getCluster());
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
			STITreeCluster c = new STITreeCluster();
			c.setCluster(bs);
//            System.err.println("m[0]: "+((STITreeCluster)minClusters.get(0)).toString2());
//            System.err.println("C: "+c.toString2());
//            System.err.println("Equals: "+((STITreeCluster)minClusters.get(0)).equals(c));
			if (c.getClusterSize() == stTaxa.length) {
				((STINode) node).setData(Double.valueOf(0));
			} else {
				int pos = minClusters.indexOf(c);                                
				((STINode) node).setData((Double) coals.get(pos));
			}
		}

		sol._totalCoals = (int) (sigmaNs - all._max_score);
		System.out.println("total cost: " + (sigmaNs - all._max_score));
		solutions.add(sol);

		return (List<Solution>) (List<Solution>) solutions;
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

}
