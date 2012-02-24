package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import phylonet.coalescent.DuplicationWeightCounter.STBipartition;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.util.Collapse;
import phylonet.tree.util.Trees;


public class MGDInference_DP {
	private static boolean _print = true;
	private static boolean optimizeDuploss = false;
	private static boolean rooted = true;

	public static void main(String[] args) {
		if ((args == null) || args.length == 0 || (args[0].equals("-h")) || (args.length < 1)) {
			printUsage();
			return;
		}
		Map taxonMap = null;
		List<Tree> trees = null;
		String output = null;
		boolean explore = false;
		double proportion = 0.0D;
		boolean exhaust = false;
		double bootstrap = 1.0D;
		double time = -1.0D;
		boolean unresolved = false;
		long startTime = System.currentTimeMillis();
		String line;
		BufferedReader treesbr = null;
		try {
			List<String[]> options = getOptions(args);
			for (String[] option : options) {
				if (option[0].equals("-i")) {
					if (option.length != 2) {
						printUsage();
						return;
					}					
					treesbr = new BufferedReader(new FileReader(
							option[1]));

					trees = new ArrayList();
					// String line;
				} else if (option[0].equals("-a")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					BufferedReader br = new BufferedReader(new FileReader(
							option[1]));

					taxonMap = new HashMap();
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
											.println("An allele can only map to one species");
									System.exit(-1);
								} else {
									taxonMap.put(allele, species);
								}
							}
						}
					}
					br.close();
				} else if (option[0].equals("-o")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					output = option[1];
				} else if (option[0].equals("-x")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					exhaust = true;
				} else if (option[0].equals("-e")) {
					if (option.length > 2) {
						printUsage();
						return;
					}
					explore = true;
					if (option.length != 2)
						continue;
					try {
						proportion = Double.parseDouble(option[1]);
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
				} else if (option[0].equals("-u")) {
					if ((option.length != 1) || (time != -1.0D)) {
						printUsage();
						return;
					}
					rooted = false;
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

			if (treesbr == null) {
				System.err.println("The input file has not been specified.");
				printUsage();
				return;
			}
			
			System.out.println("Gene trees are treated as " + (rooted?"rooted":"unrooted"));
			while ((line = treesbr.readLine()) != null) {
				Set<String> previousTreeTaxa = new HashSet<String>();
				if (line.length() > 0) {
					NewickReader nr = new NewickReader(
							new StringReader(line));
					if (rooted) {
						STITree gt = new STITree(true);
						nr.readTree(gt);							
						if (previousTreeTaxa.isEmpty()) {
							previousTreeTaxa.addAll(Arrays.asList(gt.getLeaves()));
						} else {
							if (! previousTreeTaxa.containsAll(Arrays.asList(gt.getLeaves()))) {
								throw new RuntimeException("Not all trees are on the same set of taxa: "
										+gt.getLeaves() + "\n"+ previousTreeTaxa);
							}
						}
						trees.add(gt);
					} else {
						Tree tr = nr.readTree();
						trees.add(tr);
					}
				}
			}
			treesbr.close();
		
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
			System.out.println("Reading trees in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}

		startTime = System.currentTimeMillis();
		MGDInference_DP inference = new MGDInference_DP();
		List<Solution> solutions;
		if (taxonMap == null) {
			solutions = inference.inferSpeciesTree(trees, explore, proportion,
					exhaust, bootstrap, unresolved, time);
		} else {
			//throw new RuntimeException("Not Implemented");
			solutions = inference.inferSpeciesTree(trees, taxonMap, explore,
					proportion, exhaust, bootstrap, unresolved, time);
		}

		if (_print) {
			System.out.println("Optimal tree inferred in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");			
		}

		if ((output == null) && (_print)) {
			for (Solution s : solutions)
				System.out.println(s._st.toStringWD() + " \n" + s._totalCoals
						+ (optimizeDuploss? " duplication+loss":" duplicatins") + " in total");
		} else
			try {
				FileWriter fw = new FileWriter(output);
				for (Solution s : solutions) {
					if (_print) {
						fw.write(s._st.toStringWD() + " \n" + s._totalCoals
								+ " duplications in total\n");
					} else {
						fw.write(s._st.toString() + " \n" + s._totalCoals
								+ " duplications in total\n");
					}
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
				.println("\t-i input [-e proportion] [-a mapping] [-dl] [-o output]");
		System.out
				.println("\t-i gene tree file: The file containing gene trees. (required)");
		System.out
				.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled. (optional)");
		System.out
				.println("\t-o species tree file: The file to store the species tree. (optional)");
		System.out
				.println("\t-dl optimize duploss instead of duplications");
		System.out
				.println("\t-u treat input gene trees as unrooted");				
		System.out.println();
	}

	public static void setPrinting(boolean print) {
		_print = print;
	}

	private Map<STITreeCluster,Vertex> clusterToVertex;

	public List<Solution> inferSpeciesTree(List<Tree> trees, boolean explore,
			double proportion, boolean exhaust, double bootstrap,
			boolean unresolved, double time) {
		
		long startTime = System.currentTimeMillis();
		if ((trees == null) || (trees.size() == 0)) {
			throw new IllegalArgumentException("empty or null list of trees");
		}

		if (bootstrap < 1.0D) {
			for (Tree tr : trees) {
				if (Trees.handleBootStrapInTree(tr, bootstrap) == -1) {
					throw new IllegalArgumentException(
							"Input gene trees have nodes that don't have bootstrap value");
				}
			}

		}

		Collapse.CollapseDescriptor cd = null;
		if (rooted) {			
			cd = doCollapse(trees);
		}

		List<String> taxalist = new ArrayList<String>();
		for (Tree tr : trees) {
			for (TNode node : tr.postTraverse()) {
				if ((node.isLeaf())
						&& (!taxalist.contains(node.getName()))) {
					taxalist.add(node.getName());
				}
			}
		}

		String[] taxa = new String[taxalist.size()];

		int index = 0;
		for (String taxon : taxalist) {
			taxa[(index++)] = taxon;
		}
		System.out.println("taxa size: " + index);

		Map<Integer,Set<Vertex>> clusters = new HashMap<Integer, Set<Vertex>>();
		/*if (!exhaust) {
			//SIA: we are here
			computeTreeClusters(trees, taxa, clusters);
			//computeTreeSTBipartitions(trees, taxa, clusters);
		} else {
			throw new RuntimeException("Not Implemented");
			//maxEL = computeAllClusters(trees, taxa, clusters);
		}*/
		List<Solution> solutions;
		
		DuplicationWeightCounter counter = new DuplicationWeightCounter(taxa,rooted);
		
		counter.computeTreeSTBipartitions(trees, null, clusters);
				
	
		counter.addExtraBipartitions(clusters, taxa);
		
		if (_print) {
			System.out.println("STBs formed in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}		


		counter.calculateWeights(taxa);

		if (false && _print) {
			System.out.println("Trees weights calculated in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}		

		if (explore) {
			solutions = null;//findTreesByClique(clusters, taxa, proportion);
			throw new RuntimeException("Not Implemented");
		} else {
			int maxEL;
			//SIA: we are here
			maxEL = (taxa.length - 1) * trees.size();  
			maxEL = optimizeDuploss ? maxEL + 2 * (taxa.length - 1) * trees.size() : maxEL;
			solutions = findTreesByDP(clusters, taxa, maxEL,counter,trees, null);
		}

		if (!unresolved) {
			time *= 60.0D;
			for (Solution sol : solutions) {
				if (!Trees.isBinary(sol._st)) {
					throw new RuntimeException("Where did we get unresolved trees from?");
/*					sol._totalCoals = (tryBinaryResolutions(sol._st, time,
							taxa, trees, null) + sol._totalCoals);
*/				}
			}
		}

		if (rooted) {
			restoreCollapse(solutions, cd);
		}
		return (List<Solution>) solutions;
	}

	public List<Solution> inferSpeciesTree(List<Tree> trees,
			Map<String, String> taxonMap, boolean explore, double proportion,
			boolean exhaust, double bootstrap, boolean unresolved, double time) {
		long startTime = System.currentTimeMillis();
		if ((trees == null) || (trees.size() == 0)) {
			System.err
					.println("Empty list of trees. The function returns a null tree.");
			return null;
		}

		String error = Trees.checkMapping(trees, taxonMap);
		if (error != null) {
			throw new RuntimeException("Gene trees have leaf named " + error
					+ "that hasn't been defined in the mapping file");
		}

		if (bootstrap < 1.0D) {
			for (Tree tr : trees) {
				if (Trees.handleBootStrapInTree(tr, bootstrap) == -1) {
					throw new IllegalArgumentException(
							"Input gene trees have nodes that don't have bootstrap value");
				}
			}

		}

		List temp1 = new LinkedList();
		Object temp2 = new LinkedList();
		for (String s : taxonMap.keySet()) {
			temp1.add(s);
			if (!((List) temp2).contains(taxonMap.get(s))) {
				((List) temp2).add((String) taxonMap.get(s));
			}
		}

		String[] gtTaxa = new String[temp1.size()];
		String[] stTaxa = new String[((List) temp2).size()];

		for (int i = 0; i < gtTaxa.length; i++) {
			gtTaxa[i] = ((String) temp1.get(i));
		}
		for (int i = 0; i < stTaxa.length; i++) {
			stTaxa[i] = ((String) ((List) temp2).get(i));
		}

		Map<Integer,Set<Vertex>> clusters = new HashMap<Integer, Set<Vertex>>();
		int maxEL;
		
/*		if (!exhaust) {
			maxEL = computeTreeClusters(trees, stTaxa, gtTaxa, taxonMap,
					clusters);
		} else{
			//maxEL = computeAllClusters(trees, stTaxa, taxonMap, clusters);
			throw new RuntimeException("Not Implemented");
		}
*/				
		
		List<Solution> solutions;
		
		DuplicationWeightCounter counter = new DuplicationWeightCounter(gtTaxa,stTaxa,rooted);
		
		int sigmaN = counter.computeTreeSTBipartitions(trees, taxonMap, clusters);
				
		counter.calculateWeights(stTaxa);
		
		counter.addExtraBipartitions(clusters, stTaxa);
		
		if (_print) {
			System.out.println("STBs formed in "
					+ (System.currentTimeMillis() - startTime) / 1000.0D
					+ " secs");
		}		
		
		maxEL = optimizeDuploss ? sigmaN + 2 * (stTaxa.length - 1) * trees.size() : sigmaN;
		
		if (explore) {
			solutions = null;//findTreesByClique(clusters, taxa, proportion);
			throw new RuntimeException("Not Implemented");
		} else {
			//SIA: we are here
/*			if (DUPLOSS){
				maxEL = 3*(stTaxa.length - 1) * trees.size();  // 3 is for duploss 
			} */
			// System.out.println("maxEL: " + maxEL + taxa.length + trees.size());
			solutions = findTreesByDP(clusters, stTaxa, maxEL, counter,trees, taxonMap);
		}

		if (!unresolved) {
			time *= 60.0D;
			for (Solution sol : solutions) {
				if (!Trees.isBinary(sol._st)) {
					throw new RuntimeException("Where did we get unresolved trees from?");
/*					sol._totalCoals = (tryBinaryResolutions(sol._st, time,
							stTaxa, trees, taxonMap) + sol._totalCoals);
*/				}
			}
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

	private List<Solution> findTreesByDP(Map<Integer, Set<Vertex>> clusters,
			String[] stTaxa, int sigmaN, DuplicationWeightCounter counter, List<Tree> trees, Map<String, String> taxonMap) {
		List<Solution> solutions = new ArrayList<Solution>();

		
/*		clusterToVertex = new HashMap<STITreeCluster, Vertex>();
		for (Set<Vertex> vs: clusters.values()) {
			for (Vertex vertex : vs) {
				clusterToVertex.put(vertex._cluster,vertex);
			}			
		}				
		Vertex all = (Vertex) clusters.get(Integer
				.valueOf(stTaxa.length)).toArray()[0];			
		computeMinCost(clusters, all, sigmaN, counter,trees, taxonMap);
		
		System.out.println("first round finished, adding new STBs");
		counter.addExtraBipartitions(clusters, stTaxa);
*/		
		clusterToVertex = new HashMap<STITreeCluster, Vertex>();
		for (Set<Vertex> vs: clusters.values()) {
			for (Vertex vertex : vs) {
				vertex._max_score = -1;
				clusterToVertex.put(vertex._cluster,vertex);
			}			
		}		
		
		Vertex all = (Vertex) clusters.get(Integer
				.valueOf(stTaxa.length)).toArray()[0];
		System.out.println("maxEL: " + sigmaN);
		
		computeMinCost(clusters, all, sigmaN, counter,trees, taxonMap);

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

			minClusters.add(pe._cluster);
			int k = sigmaN/(stTaxa.length-1);
			
			if (pe._min_rc != null) {
				minVertices.push(pe._min_rc);
			}
			if (pe._min_lc != null) {
				minVertices.push(pe._min_lc);
			}
			if (pe._min_lc != null && pe._min_rc != null) {
				coals.add(pe._min_cost);
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
			Object tr = new STITree();
			for (String s : stTaxa) {
				((MutableTree) tr).getRoot().createChild(s);
			}
			sol._st = ((Tree) tr);
		} else {
			sol._st = Trees.buildTreeFromClusters(minClusters);
		}

		Object map = new HashMap();
		for (TNode node : sol._st.postTraverse()) {
			BitSet bs = new BitSet();
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
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.setCluster(bs);
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

/*	private List<Solution> findTreesByClique(Map<Integer, List<Vertex>> cmap,
			String[] stTaxa, double proportion) {
		List solutions = new LinkedList();

		List clusters = new ArrayList();
		int addEL = 0;

		for (Map.Entry entry : cmap.entrySet()) {
			if (((Integer) entry.getKey()).intValue() == 1) {
				List<Vertex> l = (List) entry.getValue();
				for (Vertex v : l) {
					addEL += v._el_num;
				}
			} else if (((Integer) entry.getKey()).intValue() < stTaxa.length) {
				List<Vertex> l = (List) entry.getValue();
				for (Vertex v : l) {
					STITreeClusterWD c = new STITreeClusterWD(v._cluster);
					c.setData(Integer.valueOf(v._el_num));
					clusters.add(c);
				}
			}
		}

		double[][] compatibilityMatrix = new double[clusters.size()][clusters
				.size()];
		for (int i = 0; i < clusters.size() - 1; i++) {
			STITreeCluster cl1 = (STITreeCluster) clusters.get(i);
			for (int j = i + 1; j < clusters.size(); j++) {
				STITreeCluster cl2 = (STITreeCluster) clusters.get(j);

				if (cl1.isCompatible(cl2)) {
					compatibilityMatrix[i][j] = 1.0D;
				} else {
					compatibilityMatrix[i][j] = 0.0D;
				}
			}

		}

		MaxClique mc = new MaxClique(compatibilityMatrix);
		List<int[]> nodeCliques = mc.calculateGroups(false, 0.0D);
		List<Solution> maxCliques = new LinkedList();

		for (int max = ((STITreeClusterWD) clusters.get(0)).getTaxa().length - 2; max > 1; max--) {
			for (int[] nodes : nodeCliques) {
				if (nodes.length != max) {
					continue;
				}
				int sum = 0;
				for (int id : nodes) {
					sum += ((Integer) ((STITreeClusterWD) clusters.get(id))
							.getData()).intValue();
				}
				Solution s = new Solution();
				s._clusterIDs = nodes;
				s._totalCoals = (sum + addEL);
				maxCliques.add(s);
			}
			if (maxCliques.size() > 0)
				break;
		}
		Solution s2;
		for (int i = 1; i < maxCliques.size(); i++) {
			Solution s1 = (Solution) maxCliques.get(i);
			for (int j = 0; j < i; j++) {
				s2 = (Solution) maxCliques.get(j);
				if (s1._totalCoals < s2._totalCoals) {
					maxCliques.remove(s1);
					maxCliques.add(j, s1);
					break;
				}
			}

		}

		int minCoal = ((Solution) maxCliques.get(0))._totalCoals;
		int maxCoal = (int) ((1.0D + proportion / 100.0D) * minCoal);

		Solution s;
		for (Iterator iterator4 = maxCliques.iterator(); iterator4.hasNext(); solutions
				.add(s)) {
			s = (Solution) iterator4.next();
			if (s._totalCoals > maxCoal)
				break;
			List minClusters = new ArrayList();
			List coals = new ArrayList();
			int ai1[];
			int k1 = (ai1 = s._clusterIDs).length;
			for (int j1 = 0; j1 < k1; j1++) {
				int id = ai1[j1];
				STITreeCluster c = (STITreeCluster) clusters.get(id);
				minClusters.add(c);
				coals.add((Integer) ((STITreeClusterWD) c).getData());
			}

			Vertex v;
			for (Iterator iterator5 = ((List) cmap.get(Integer.valueOf(1)))
					.iterator(); iterator5.hasNext(); coals.add(Integer
					.valueOf(v._el_num))) {
				v = (Vertex) iterator5.next();
				minClusters.add(v._cluster);
			}

			Tree st = Trees.buildTreeFromClusters(minClusters);
			Map map = new HashMap();
			for (Iterator iterator6 = st.postTraverse().iterator(); iterator6
					.hasNext();) {
				TNode node = (TNode) iterator6.next();
				BitSet bs = new BitSet();
				if (node.isLeaf()) {
					for (int i = 0; i < stTaxa.length; i++) {
						if (!node.getName().equals(stTaxa[i]))
							continue;
						bs.set(i);
						break;
					}

					map.put(node, bs);
				} else {
					BitSet childCluster;
					for (Iterator iterator7 = node.getChildren().iterator(); iterator7
							.hasNext(); bs.or(childCluster)) {
						TNode child = (TNode) iterator7.next();
						childCluster = (BitSet) map.get(child);
					}

					map.put(node, bs);
				}
				STITreeCluster c = new STITreeCluster(stTaxa);
				c.setCluster(bs);
				if (c.getClusterSize() == stTaxa.length) {
					((STINode) node).setData(Integer.valueOf(0));
				} else {
					int pos = minClusters.indexOf(c);
					((STINode) node).setData((Integer) coals.get(pos));
				}
			}

			s._st = st;
		}

		return solutions;
	}
*/
/*	private int computeTreeClusters(List<Tree> trees, String[] stTaxa,
			String[] gtTaxa, Map<String, String> taxonMap,
			Map<Integer, Set<Vertex>> clusters) {
		int maxEL = 0;
		for (Tree tr : trees) {
			for (Iterator localIterator2 = tr.getClusters(gtTaxa, false)
					.iterator(); localIterator2.hasNext();) {
				STITreeCluster tc = (STITreeCluster) localIterator2.next();

				STITreeCluster stCluster = new STITreeCluster(stTaxa);
				for (String s : tc.getClusterLeaves()) {
					stCluster.addLeaf((String) taxonMap.get(s));
				}

				int csize = stCluster.getClusterSize();

				if (csize <= 1) {
					continue;
				}
				//me
			    if (csize > 1) {
					List<Vertex> l = new LinkedList();
					Vertex v = new Vertex();

					((Vertex) v)._cluster = stCluster;
					((Vertex) v)._el_num = DeepCoalescencesCounter
							.getClusterCoalNum(trees, stCluster, taxonMap, true);
					if (((Vertex) v)._el_num > maxEL) {
						maxEL = ((Vertex) v)._el_num;
					}
					((Vertex) v)._min_cost = -1;

					((List) l).add(v);
					clusters.put(Integer.valueOf(csize), l);
				}
				//me
				if (! clusters.containsKey(Integer.valueOf(csize))) {
					List<Vertex> l = new LinkedList();
					clusters.put(Integer.valueOf(csize), l);

				}
				List<Vertex> l = (List) clusters
						.get(Integer.valueOf(csize));
				boolean found = false;

				for (Vertex v : l) {
					if (v._cluster.equals(stCluster)) {
						found = true;
						break;
					}
				}

				if (!found) {
					Vertex nv = new Vertex();
					nv._cluster = stCluster;
					//me
					nv._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, stCluster, taxonMap, true);
					if (nv._el_num > maxEL) {
						maxEL = nv._el_num;
					}
					//me
					nv._min_cost = -1;

					l.add(nv);
				}
			} 

		}

		STITreeCluster all = new STITreeCluster(stTaxa);
		String as[];
		int j = (as = stTaxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i];
			all.addLeaf(t);
		}

		Vertex v = new Vertex();
		v._cluster = all;
		v._el_num = 0;  //me
		v._min_cost = -1;

		List la = new LinkedList();
		la.add(v);

		clusters.put(Integer.valueOf(all.getClusterSize()), la);

		List<Vertex> l1 = new LinkedList<Vertex>();
		for (String t : stTaxa) {
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.addLeaf(t);

			v = new Vertex();
			v._cluster = c;
			//me
			v._el_num = DeepCoalescencesCounter.getClusterCoalNum(trees, c,
					taxonMap, true);
			if (v._el_num > maxEL) {
				maxEL = v._el_num;
			}
			//me
			v._min_cost = -1;

			l1.add(v);
		}
		maxEL++;
		clusters.put(Integer.valueOf(1), l1);
		return maxEL;
		
		int maxEL = 0;
		// SIA: here they simply find internal vertices, and add
		// them (with their weights) to *clusters* hashmap
		for (Tree tr : trees) {
			maxEL += tr.getLeafCount() - 1;
			for (Iterator clustersIt = tr.getClusters(gtTaxa, false)
					.iterator(); clustersIt.hasNext();) {
				STITreeCluster tc = (STITreeCluster) clustersIt.next();
				
				STITreeCluster stCluster = new STITreeCluster(stTaxa);
				for (String s : tc.getClusterLeaves()) {
					stCluster.addLeaf((String) taxonMap.get(s));
				}

				int csize = stCluster.getClusterSize();
											

				if (! clusters.containsKey(Integer.valueOf(csize))) {
					clusters.put(Integer.valueOf(csize), new HashSet<Vertex>());
				}
				Set<Vertex> l = clusters.get(Integer
						.valueOf(csize));

				Vertex nv = new Vertex();
				
				nv._cluster = stCluster;
				
				//STITreeCluster tc = new STITreeCluster(taxa);
				
				if (!l.contains(nv)) {
					//me
					nv._el_num = -1;
					nv._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, tc, true);
					if (nv._el_num > maxEL) {
						maxEL = nv._el_num;
					}
					//me
					nv._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, tc, true);
					if (nv._el_num > maxEL) {
						maxEL = nv._el_num;
					}
					nv._min_cost = -1;

					l.add(nv);	
				}
			}

		}

		// SIA: and add a node for the root. It weight is 0
		STITreeCluster all = new STITreeCluster(stTaxa);	
		String as[];
		int j = (as = stTaxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i];
			all.addLeaf(t);
		}
		Vertex v = new Vertex();
		v._cluster = all;
		v._el_num = -1; 
		v._min_cost = -1;

		HashSet la = new HashSet();
		la.add(v);

		clusters.put(Integer.valueOf(all.getClusterSize()), la);
		
		// SIA: And finally add nodes for singelton clades
		HashSet<Vertex> l1 = new HashSet<Vertex>();
		for (String t : stTaxa) {
			STITreeCluster c = new STITreeCluster(stTaxa);
			c.addLeaf(t);
			v = new Vertex();
			v._cluster = c;
			//me
			v._el_num = -1;			
			//me
			v._min_cost = -1;

			l1.add(v);
		}
		
		clusters.put(Integer.valueOf(1), l1);
		int s = 0;
		for (Integer c: clusters.keySet()){
			s += clusters.get(c).size();
		}
		System.out.println("Number of Clusters considered: " +s);
		
		return maxEL;
	}
*/
	/*protected int computeTreeClusters(List<Tree> trees, String[] taxa,
			Map<Integer, Set<Vertex>> clusters) {
		
		int maxEL = 0;
		
		// SIA: here they simply find internal vertices, and add
		// them (with their weights) to *clusters* hashmap
		for (Tree tr : trees) {
			for (Iterator clustersIt = (
					rooted ? 
					tr.getClusters(taxa, false):tr.getBipartitionClusters(taxa,false))
					.iterator(); clustersIt.hasNext();) {
				STITreeCluster tc = (STITreeCluster) clustersIt.next();
				int tc_size = tc.getClusterSize();			

				if (! clusters.containsKey(Integer.valueOf(tc_size))) {
					clusters.put(Integer.valueOf(tc_size), new HashSet());
				}
				Set<Vertex> l = clusters.get(Integer
						.valueOf(tc_size));

				Vertex nv = new Vertex();
				nv._cluster = tc;
				
				if (!l.contains(nv)) {
					nv._el_num = -1;
					nv._min_cost = -1;

					l.add(nv);	
				}
			}

		}

		// SIA: and add a node for the root. It weight is 0
		STITreeCluster all = new STITreeCluster(taxa);
		String as[];
		int j = (as = taxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i];
			all.addLeaf(t);
		}
		Vertex v = new Vertex();
		v._cluster = all;
		v._el_num = -1; //me
		v._min_cost = -1;
		HashSet la = new HashSet();
		la.add(v);
		clusters.put(Integer.valueOf(all.getClusterSize()), la);

		// SIA: And finally add nodes for singelton clades
		HashSet<Vertex> l1 = new HashSet();
		String as1[];
		int i1 = (as1 = taxa).length;
		for (int k = 0; k < i1; k++) {
			String t = as1[k];
			STITreeCluster c = new STITreeCluster(taxa);
			c.addLeaf(t);
			v = new Vertex();
			v._cluster = c;
			v._el_num = -1; //me
			v._min_cost = -1;
			l1.add(v);
		}
		maxEL++;
		clusters.put(Integer.valueOf(1), l1);
		int s = 0;
		for (Integer c: clusters.keySet()){
			s += clusters.get(c).size();
		}
		System.out.println("Number of Clusters considered: " +s);
		
		//System.out.println(clusters);
		return maxEL;
	}
*/
	private int computeMinCost(Map<Integer, Set<Vertex>> clusters, Vertex v,
			int sigmaNs, DuplicationWeightCounter counter, List<Tree> trees, Map<String, String> taxonMap) {
		int maxEL = 1000000000;
		// SIA: Already calculated. Don't re-calculate.
		if (v._max_score != -1) {
			return v._max_score - maxEL;
		}
		// SIA: If in duploss mode, need to get MDC cost as well
		if (optimizeDuploss) {			
			if (v._el_num == -1) {
				if (taxonMap == null) {
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, rooted);
				} else {
					v._el_num = DeepCoalescencesCounter.getClusterCoalNum(
							trees, v._cluster, taxonMap, rooted);
				}
			}
		} else {
			v._el_num = 0;
		}
		// SIA: base case for singelton clusters.
		int clusterSize = v._cluster.getClusterSize();
		if (clusterSize <= 1) {
			// SIA: TODO: this is 0, right?
			v._min_cost = 0;
			v._max_score = maxEL - v._el_num;
			v._min_lc = (v._min_rc = null);			
			return v._max_score - maxEL;
		}
		
		Set<STBipartition> clusterBiPartitions = counter.getClusterBiPartitions(v._cluster);
		
		//STBipartition bestSTB = null;
				
		for (STBipartition stb: clusterBiPartitions) {

			Vertex lv = clusterToVertex.get(stb.cluster1);
			Vertex rv = clusterToVertex.get(stb.cluster2);
			
			if (lv == null || rv == null) {
				//System.out.println("There is no STB for one half of : " + stb);
				continue;
			}

			int lscore = computeMinCost(clusters, lv, sigmaNs, counter,trees,taxonMap);
			int rscore = computeMinCost(clusters, rv, sigmaNs, counter,trees,taxonMap);

			int w = counter.getBiPartitionDPWeight(lv._cluster, rv._cluster, v._cluster);

			int z = optimizeDuploss? 3 : 1;

			int c = z*w - v._el_num;

			if ((v._max_score != -1)
					&& (lscore + rscore + c + maxEL
							<= v._max_score)){
				continue;
			}							
			v._max_score = (lscore + rscore + c) + maxEL;
			v._min_cost = sigmaNs - (c + lv._max_score + rv._max_score - 2*maxEL);
			//stem.out.println(maxEL - (z*w + lv._max_score + rv._max_score));
			v._min_lc = lv;
			v._min_rc = rv;				
			//bestSTB = stb;
		}

/*		if (clusterSize > 5){
			counter.addGoodSTB(bestSTB, clusterSize);
		}
*/	
		return v._max_score - maxEL;
	}

	
	
/*	private int computeAllClusters(List<Tree> trees, String[] stTaxa,
			Map<String, String> taxonMap, Map<Integer, List<Vertex>> clusters) {
		int n = stTaxa.length;
		int maxEL = 0;
		if (n <= 0) {
			System.err.println("Empty list of taxa.");
			return -1;
		}

		BitSet counter = new BitSet(n);
		boolean done = false;

		while (!done) {
			int i = 0;
			while ((i < n) && (counter.get(i))) {
				counter.clear(i);
				i++;
			}
			if (i >= n) {
				done = true;
			} else {
				counter.set(i, true);

				STITreeCluster tc = new STITreeCluster(stTaxa);
				tc.setCluster((BitSet) counter.clone());
				Vertex v = new Vertex();
				v._cluster = tc;
				if (taxonMap == null) {
					v._el_num = DuplicationWeightCounter.getClusterCoalNum(
							trees, tc, true);
				} else {
					v._el_num = DuplicationWeightCounter.getClusterCoalNum(
							trees, tc, taxonMap, true);
				}
				if (v._el_num > maxEL) {
					maxEL = v._el_num;
				}

				int size = tc.getClusterSize();
				List l = (List) clusters.get(Integer.valueOf(size));
				if (l == null) {
					l = new LinkedList();
				}

				l.add(v);
				clusters.put(Integer.valueOf(size), l);
			}
		}
		maxEL++;
		return maxEL;
	}*/

	private int computeAllClusters(List<Tree> trees, String[] stTaxa,
			Map<Integer, List<Vertex>> clusters) {
		throw new RuntimeException("Not Implemented");
		//return computeAllClusters(trees, stTaxa, null, clusters);
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

/*	private int tryBinaryResolutions(Tree tr, double time, String[] taxa,
			List<Tree> gts, Map<String, String> taxonMap) {
		List nodelist = new ArrayList();
		List degreelist = new ArrayList();
		int totalResolutions = 0;
		for (Iterator iterator = (new PostTraversal(tr.getRoot())).iterator(); iterator
				.hasNext();) {
			TNode node = (TNode) iterator.next();
			int childCount = node.getChildCount();
			if (childCount > 2) {
				nodelist.add(node);
				int resolutionsNumber = getResolutionsNumber(childCount);
				degreelist.add(Integer.valueOf(resolutionsNumber));
				totalResolutions += resolutionsNumber;
			}
		}
		int addedxl = 0;
		for (int i = 0; i < nodelist.size(); i++) {
			TNode unresolvedNode = (TNode) nodelist.get(i);
			Map id2node = new HashMap();
			for (TNode child : unresolvedNode.getChildren()) {
				id2node.put(Integer.valueOf(child.getID()), child);
			}
			Integer[] childIDs = (Integer[]) id2node.keySet().toArray(
					new Integer[0]);
			Object cluster2xl = new HashMap();
			double endtime;
			if (time == -1.0D) {
				endtime = -1.0D;
			} else {
				endtime = time
						* (((Integer) degreelist.get(i)).intValue() / totalResolutions)
						* 1000.0D + System.currentTimeMillis();
			}
			Solution sol = addMoreLeaves(null, childIDs, 0, id2node, taxa, gts,
					taxonMap, endtime, (Map) cluster2xl);
			TNode parent = unresolvedNode.getParent();
			int xl = ((Integer) ((STINode) unresolvedNode).getData())
					.intValue();
			((STINode) unresolvedNode).removeNode();
			if (parent != null) {
				TNode newnode = ((STINode) parent).createChild(sol.getTree()
						.getRoot());
				((STINode) newnode).setData(Integer.valueOf(xl));
			} else {
				for (TNode child : sol.getTree().getRoot().getChildren()) {
					((STINode) tr.getRoot()).createChild(child);
				}
			}
			addedxl += sol.getCoalNum();
		}
		return addedxl;
	}*/

/*	private Solution addMoreLeaves(STITree<Integer> preTree,
			Integer[] leavesid, int index, Map<Integer, TNode> id2node,
			String[] taxa, List<Tree> gts, Map<String, String> taxonMap,
			double endTime, Map<BitSet, Integer> cluster2xl) {
		if (preTree == null) {
			preTree = new STITree(false);
			STINode root = preTree.getRoot();
			STINode newnode = root.createChild();
			newnode.setData(leavesid[(index++)]);
			STINode innode = root.createChild();
			newnode = innode.createChild();
			newnode.setData(leavesid[(index++)]);
			newnode = innode.createChild();
			newnode.setData(leavesid[(index++)]);
		}
		Solution sol = null;
		if (index == leavesid.length) {
			sol = tryAllRootings(preTree, id2node, taxa, gts, taxonMap,
					endTime, cluster2xl);
		} else {
			int id = leavesid[(index++)].intValue();
			for (Iterator iterator = (new PostTraversal(preTree.getRoot()))
					.iterator(); iterator.hasNext();) {
				TNode n = (TNode) iterator.next();
				if (!n.isLeaf()) {
					if (n.getChildCount() != 2) {
						throw new RuntimeException("Not binary!");
					}
					Iterator it = n.getChildren().iterator();
					for (int i = 0; i < 2; i++) {
						STINode child = (STINode) it.next();
						STITree newTree = new STITree(preTree);
						TNode peerChild = newTree.getNode(child.getID());
						TNode peerParent = peerChild.getParent();
						STINode newchild = ((STINode) peerParent).createChild();
						newchild.adoptChild((TMutableNode) peerChild);
						STINode newnode = newchild.createChild();
						newnode.setData(Integer.valueOf(id));
						Solution thissol = addMoreLeaves(newTree, leavesid,
								index, id2node, taxa, gts, taxonMap, endTime,
								cluster2xl);
						if ((sol == null)
								|| (sol.getCoalNum() > thissol.getCoalNum())) {
							sol = thissol;
						}
						double now = System.currentTimeMillis();
						if ((endTime != -1.0D) && (now > endTime)) {
							return sol;
						}
						if (n.isRoot()) {
							break;
						}
					}
				}
			}
		}
		return sol;
	}*/

/*	private Solution tryAllRootings(Tree subtree, Map<Integer, TNode> id2node,
			String[] taxa, List<Tree> gts, Map<String, String> taxonMap,
			double endTime, Map<BitSet, Integer> cluster2xl) {
		Solution sol = new Solution();
		sol._totalCoals = -1;
		for (Tree rootedsubtree : subtree.getAllRootingTrees()) {
			TNode peerNode;
			for (Iterator iterator1 = (new PostTraversal(rootedsubtree
					.getRoot())).iterator(); iterator1.hasNext();) {
				TNode replacingNode = (TNode) iterator1.next();
				if (replacingNode.isLeaf()) {
					int id = ((Integer) ((STINode) replacingNode).getData())
							.intValue();
					peerNode = (TNode) id2node.get(Integer.valueOf(id));
					if (peerNode.isLeaf()) {
						((STINode) replacingNode).setName(peerNode.getName());
					} else {
						for (TNode child : peerNode.getChildren()) {
							((STINode) replacingNode).createChild(child);
						}
						((STINode) replacingNode).setName("");
					}

					((STINode) replacingNode)
							.setData((Integer) ((STINode) peerNode).getData());
				}

			}

			int xl = 0;
			Object map = new HashMap();
			for (Iterator iterator2 = (new PostTraversal(rootedsubtree
					.getRoot())).iterator(); iterator2.hasNext();) {
				TNode node = (TNode) iterator2.next();
				BitSet bs = new BitSet();
				if (node.isLeaf()) {
					for (int i = 0; i < taxa.length; i++) {
						if (node.getName().equals(taxa[i])) {
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
				if (!node.isRoot()) {
					Integer el = (Integer) ((STINode) node).getData();
					if (el == null) {
						el = (Integer) cluster2xl.get(bs);
						if (el == null) {
							STITreeCluster tc = new STITreeCluster(taxa);
							tc.setCluster(bs);
							if (taxonMap == null) {
								el = Integer.valueOf(DuplicationWeightCounter
										.getClusterCoalNum(gts, tc, true));
							} else {
								el = Integer.valueOf(DuplicationWeightCounter
										.getClusterCoalNum(gts, tc, taxonMap,
												true));
							}
							cluster2xl.put(bs, el);
						}
						((STINode) node).setData(el);
						xl += el.intValue();
					}
				}
			}

			if ((sol.getCoalNum() == -1) || (sol.getCoalNum() > xl)) {
				sol._totalCoals = xl;
				sol._st = rootedsubtree;
			}
		}
		return (Solution) sol;
	}*/

	static class Vertex {
		public STITreeCluster _cluster = null;
		public  int _el_num = -1;
		public int _min_cost = -1;
		public int _max_score = -1;
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
			return ((Vertex)obj)._cluster.equals(this._cluster);
		}

		@Override
		public int hashCode() {
			return _cluster.hashCode();
		}

	}
}