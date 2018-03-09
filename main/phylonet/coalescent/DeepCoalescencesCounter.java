package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import phylonet.network.io.ExNewickException;
import phylonet.network.io.ExNewickReader;
import phylonet.network.model.NetNode;
import phylonet.network.model.Network;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;

public class DeepCoalescencesCounter {
	public static void main(String[] args) {
		if ((args == null) || (args[0].equals("-h")) || (args.length < 2)) {
			printUsage();
			return;
		}
		List<Tree> speciesTrees = null;
		Network net = null;
		List geneTrees = new ArrayList();
		Map taxonMap = null;
		boolean rooted = true;
		double bootstrap = 1.0D;
		boolean isTree = true;
		String line;
		try {
			List<String []> options = getOptions(args);
			BufferedReader br = new BufferedReader(new FileReader(args[1]));
			STITree gt;
			while ((line = br.readLine()) != null) {
				//String line;
				line.trim();
				if (line.length() > 0) {
					NewickReader nr = new NewickReader(new StringReader(line));
					gt = new STITree(true);
					nr.readTree(gt);
					geneTrees.add(gt);
				}
			}
			br.close();
			for (String[] option : options) {
				if (option[0].equals("-n")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					isTree = false;
				} else if (option[0].equals("-a")) {
					if (option.length != 2) {
						printUsage();
						return;
					}
					taxonMap = new HashMap();
					br = new BufferedReader(new FileReader(option[1]));
					while ((line = br.readLine()) != null) {
						String[] mapString = line.split(";");
						for (String s : mapString) {
							String species = s.substring(0, s.indexOf(":"))
									.trim();
							s = s.substring(s.indexOf(":") + 1);
							String[] alleles = s.split(",");
							for (String allele : alleles) {
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

				} else if (option[0].equals("-u")) {
					if (option.length != 1) {
						printUsage();
						return;
					}
					rooted = false;
				} else {
					printUsage();
					return;
				}
			}

			if (isTree) {
				speciesTrees = new ArrayList();
				br = new BufferedReader(new FileReader(args[0]));
				while ((line = br.readLine()) != null) {
					line.trim();
					if (line.length() > 0) {
						NewickReader nr = new NewickReader(new StringReader(
								line));
						speciesTrees.add(nr.readTree());
					}
				}
				br.close();
			} else {
				ExNewickReader enr = new ExNewickReader(new FileReader(args[0]));
				net = enr.readNetwork();
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
		} catch (ExNewickException e) {
			System.err
					.println("Error when parsing the Newick representation from input file.");
			System.err.println(e.getMessage());
			e.printStackTrace();
			return;
		}

		if (isTree) {
			if (taxonMap == null) {
				int index = 1;
				for (Tree st : speciesTrees) {
					int coalNum = countExtraCoal(geneTrees, st, rooted,
							bootstrap);
					System.out.println("Species_Tree#" + index++ + " = "
							+ st.toStringWD());
					System.out.println("Total number of extra lineages: "
							+ coalNum);
				}
			} else {
				int index = 1;
				for (Tree st : speciesTrees) {
					int coalNum = countExtraCoal(geneTrees, st, taxonMap,
							rooted, bootstrap);
					System.out.println("Species_Tree#" + index++ + " = "
							+ st.toStringWD());
					System.out.println("Total number of extra lineages: "
							+ coalNum);
				}
			}
		} else
			countExtraCoal(geneTrees, net, taxonMap);
	}

	public static void printUsage() {
		System.out.println();
		System.out.println("This tool counts the number of deep coalescences");
		System.out.println("Usage is:");
		System.out
				.println("\tdeep_coal_count species-tree-file gene-trees-file [-u] [-b threshold] [-a mapping]");
		System.out
				.println("\tspecies-tree-file: File that contains species trees(required)");
		System.out
				.println("\tgene-trees-file: File that contains gene trees(required)");
		System.out
				.println("\t-a mapping file: The file containing the mapping from alleles to speceis if multiple alleles sampled (optional)");
		System.out
				.println("\t-u unrooted: Specify the gene trees and species tree to be treated as unrooted(optional)");
		System.out.println("\t-b bootstrap threshold(optional)");
	}

	public static int countExtraCoal(List<Tree> gts, Tree st, boolean rooted,
			double bootstrap) {
		int sum = 0;
		String[] taxa = st.getLeaves();

		if (bootstrap < 1.0D) {
			for (Tree tr : gts) {
				if (Trees.handleBootStrapInTree(tr, bootstrap) == -1) {
					throw new IllegalArgumentException(
							"Input gene trees have nodes that don't have bootstrap value");
				}
			}
		}

		Map map = new HashMap();
		for (TNode node : st.postTraverse()) {
			BitSet bs = new BitSet();
			if (node.isLeaf()) {
				for (int i = 0; i < taxa.length; i++) {
					if (node.getName().equals(taxa[i])) {
						bs.set(i);
						break;
					}
				}
				map.put(node, bs);
				((STINode) node).setData(Integer.valueOf(0));
			} else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = (BitSet) map.get(child);
					bs.or(childCluster);
				}
				map.put(node, bs);
				STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
				c.setCluster(bs);
				if (c.getClusterSize() == taxa.length) {
					((STINode) node).setData(Integer.valueOf(0));
				} else {
					int el = getClusterCoalNum(gts, c, rooted);
					((STINode) node).setData(Integer.valueOf(el));
					sum += el;
				}
			}
		}
		return sum;
	}

	public static int countExtraCoal(List<Tree> gts, Tree st,
			Map<String, String> taxonMap, boolean rooted, double bootstrap) {
		String error = Trees.checkMapping(gts, taxonMap);
		if (error != null) {
			throw new RuntimeException("Gene trees have leaf named " + error
					+ "that hasn't been defined in the mapping file");
		}

		int sum = 0;
		String[] stTaxa = st.getLeaves();

		if (bootstrap < 1.0D) {
			for (Tree tr : gts) {
				if (Trees.handleBootStrapInTree(tr, bootstrap) == -1) {
					throw new IllegalArgumentException(
							"Input gene trees have nodes that don't have bootstrap value");
				}
			}

		}

		Map map = new HashMap();
		for (TNode node : st.postTraverse()) {
			BitSet bs = new BitSet();
			if (node.isLeaf()) {
				for (int i = 0; i < stTaxa.length; i++) {
					if (node.getName().equals(stTaxa[i])) {
						bs.set(i);
						break;
					}
				}
				map.put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet childCluster = (BitSet) map.get(child);
					bs.or(childCluster);
				}
				map.put(node, bs);
			}
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.setCluster(bs);
			if (c.getClusterSize() == stTaxa.length) {
				((STINode) node).setData(Integer.valueOf(0));
			} else {
				throw new RuntimeException("Not implemented");
//				int el = getClusterCoalNum(gts, c, null, rooted);
//				((STINode) node).setData(Integer.valueOf(el));
//				sum += el;
			}

		}

		return sum;
	}

	public static int countExtraCoal(List<Tree> gts, Network net,
			Map<String, String> taxonMap) {
		int coal_sum = 0;
		Map nname2tamount = new HashMap();
		Tree superst = networkToTree(net, nname2tamount);
		System.out.println(superst);
		for (Tree gt : gts) {
			List<String> gt_taxa = Arrays.asList(gt.getLeaves());
			if (taxonMap == null) {
				taxonMap = new HashMap();
				for (String taxon : gt_taxa) {
					taxonMap.put(taxon, taxon);
				}
			}
			List<Map> allmappings = new ArrayList();
			Map firstmap = new HashMap();
			((Map) firstmap).putAll(taxonMap);
			allmappings.add(firstmap);
			List<Map> temp;
			for (String gtleaf : gt_taxa) {
				String nleaf = (String) taxonMap.get(gtleaf);
				temp = new ArrayList();
				temp.addAll(allmappings);
				allmappings.clear();
				for (int j = 1; j <= ((Integer) nname2tamount.get(nleaf))
						.intValue(); j++) {
					String st_leaf = nleaf + "_" + j;
					for (Map mapping : temp) {
						Map new_mapping = new HashMap();
						new_mapping.putAll(mapping);
						new_mapping.put(gtleaf, st_leaf);
						allmappings.add(new_mapping);
					}
				}
			}
			int min_coal = 2147483647;
			Object gtlist = new ArrayList();
			for (Map mapping : allmappings) {
				((List) gtlist).clear();
				((List) gtlist).add(gt);
				int coal = countExtraCoal((List) gtlist, superst, mapping,
						true, 1.0D);
				if (min_coal > coal) {
					min_coal = coal;
				}
			}
			System.out.println(gt);
			System.out.println(min_coal);
		}
		return coal_sum;
	}
	static int counter = 0;
	public static int getClusterCoalNum(List<Tree> trees,
			STITreeCluster cluster, boolean rooted) {
		int weight = 0;
		counter ++;
		for (Tree tr : trees) {
			if (rooted) {
				weight += getClusterCoalNum_rooted(tr, cluster);
			} else {
				weight += getClusterCoalNum_unrooted(tr, cluster);
			}
		}

		return weight;
	}

	public static int getClusterCoalNumMap(List<Tree> trees,
			STITreeCluster cluster, boolean rooted) {
		int weight = 0;

		for (Tree tr : trees) {
			if (rooted) {
				weight += getClusterCoalNum_rootedMap(tr, cluster);
			} else {
				weight += getClusterCoalNum_unrootedMap(tr, cluster);
			}
		}

		return weight;
	}

	public static int getClusterCoalNum_rooted(Tree tr, STITreeCluster cluster) {
		Map map = new HashMap();
		
		int count = 0;
		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				int index = GlobalMaps.taxonIdentifier.taxonId(node.getName());
				BitSet bs = new BitSet();

				bs.set(index);
				if (cluster.containsCluster(bs)) {
					count++;
				}

				map.put(node, bs);
			} else {
				BitSet bs = new BitSet();
				int intersect = 0;
				int childCount = node.getChildCount();
				for (TNode child : node.getChildren()) {
					BitSet v = (BitSet) map.get(child);
					bs.or(v);
					if ((!cluster.containsCluster(v)))
						continue;
					intersect++; //0,1, or two 
				}

				if (cluster.containsCluster(bs)) {
					count -= node.getChildCount();
					count++;
				} else if (intersect > 1) {
					count -= intersect;
					count++;
				}

				map.put(node, bs);
			}
		}
		return Math.max(count - 1, 0);
	}

	public static int getClusterCoalNum_unrooted(Tree tr, STITreeCluster cluster) {
		Map map = new HashMap();
		List taxalist = new ArrayList();
		String[] taxa = tr.getLeaves();
		int ntaxa = taxa.length;
		for (String leaf : taxa) {
			taxalist.add(leaf);
		}
		STITreeCluster concluster = new STITreeCluster(GlobalMaps.taxonIdentifier);
		for (Integer leaf : cluster) {
			if (taxalist.contains(GlobalMaps.taxonIdentifier.getTaxonName(leaf))) {
				concluster.addLeaf(leaf);
			}
		}
		if (concluster.getClusterSize() == ntaxa) {
			return 0;
		}
		Object coveragelist = new ArrayList();
		for (int i = concluster.getBitSet().nextSetBit(0); i >= 0; i = concluster
				.getBitSet().nextSetBit(i + 1)) {
			BitSet bs = new BitSet(ntaxa);
			((BitSet) bs).set(i);
			((List) coveragelist).add(bs);
		}
		int intersect;
		BitSet virtualbs;
		for (Object tn = tr.postTraverse().iterator(); ((Iterator) tn)
				.hasNext();) {
			TNode node = (TNode) ((Iterator) tn).next();
			if (((List) coveragelist).size() <= 1) {
				break;
			}
			BitSet bs = new BitSet(ntaxa);
			intersect = 0;
			virtualbs = new BitSet(ntaxa);
			if (node.isLeaf()) {
				int index = taxalist.indexOf(node.getName());
				bs.set(index);
				map.put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet v = (BitSet) map.get(child);
					bs.or(v);
					if (concluster.containsCluster(v)) {
						intersect++;
						virtualbs.or(v);
					}
				}

				if ((concluster.containsCluster(bs)) || (intersect > 1)) {
					if (concluster.containsCluster(bs)) {
						virtualbs = bs;
					}
					for (int i = 0; i < ((List) coveragelist).size(); i++) {
						BitSet exbs = (BitSet) ((List) coveragelist).get(i);
						BitSet temp = (BitSet) virtualbs.clone();
						temp.and(exbs);
						if (temp.equals(exbs)) {
							((List) coveragelist).remove(i);
							i--;
						}
					}
					((List) coveragelist).add(virtualbs);
				}

				map.put(node, bs);
			}

			if (!node.isRoot()) {
				BitSet complementbs = (BitSet) bs.clone();
				complementbs.flip(0, taxa.length);
				if (intersect > 0) {
					complementbs.or(virtualbs);
				}
				if (concluster.containsCluster(complementbs)) {
					for (int i = 0; i < ((List) coveragelist).size(); i++) {
						BitSet exbs = (BitSet) ((List) coveragelist).get(i);
						BitSet temp = (BitSet) complementbs.clone();
						temp.and(exbs);
						if (temp.equals(exbs)) {
							((List) coveragelist).remove(i);
							i--;
						}
					}
					((List) coveragelist).add(complementbs);
					break;
				}

			}

		}

		return Math.max(0, ((List) coveragelist).size() - 1);
	}

	public static int getClusterCoalNum_rootedMap(Tree tr, STITreeCluster cluster) {
		Map map = new HashMap();
		int count = 0;
		for (TNode node : tr.postTraverse()) {
			if (node.isLeaf()) {
				String stTaxon = GlobalMaps.taxonNameMap.getTaxonName(node.getName());
				int index = GlobalMaps.taxonIdentifier.taxonId(stTaxon);
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				bs.set(index);
				if (cluster.containsCluster(bs)) {
					count++;
				}

				map.put(node, bs);
			} else {
				BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
				int intersect = 0;
				int childCount = node.getChildCount();
				for (TNode child : node.getChildren()) {
					BitSet v = (BitSet) map.get(child);
					bs.or(v);
					if ((childCount <= 2) || (!cluster.containsCluster(v)))
						continue;
					intersect++;
				}

				if (cluster.containsCluster(bs)) {
					count -= node.getChildCount();
					count++;
				} else if (intersect > 1) {
					count -= intersect;
					count++;
				}

				map.put(node, bs);
			}
		}
		return Math.max(count - 1, 0);
	}

	public static int getClusterCoalNum_unrootedMap(Tree tr,
			STITreeCluster cluster) {
		Map map = new HashMap();
		List gtTaxalist = new ArrayList();
		String[] gtTaxa = tr.getLeaves();
		int ngtTaxa = gtTaxa.length;
		for (String leaf : gtTaxa) {
			gtTaxalist.add(leaf);
		}
		STITreeCluster concluster = new STITreeCluster(GlobalMaps.taxonIdentifier);
		for (TNode n : tr.getNodes()) {
			if ((!n.isLeaf())
					|| (!cluster.containsLeaf(GlobalMaps.taxonNameMap.getTaxonName(n.getName()))))
				continue;
			concluster.addLeaf(GlobalMaps.taxonIdentifier.taxonId(n.getName()));
		}

		Object coveragelist = new ArrayList();
		for (int i = concluster.getBitSet().nextSetBit(0); i >= 0; i = concluster
				.getBitSet().nextSetBit(i + 1)) {
			BitSet bs = new BitSet(ngtTaxa);
			((BitSet) bs).set(i);
			((List) coveragelist).add(bs);
		}
		int intersect;
		BitSet virtualbs;
		for (Object tn = tr.postTraverse().iterator(); ((Iterator) tn)
				.hasNext();) {
			TNode node = (TNode) ((Iterator) tn).next();
			if (((List) coveragelist).size() <= 1) {
				break;
			}
			BitSet bs = new BitSet(ngtTaxa);
			intersect = 0;
			virtualbs = new BitSet(ngtTaxa);
			if (node.isLeaf()) {
				int index = gtTaxalist.indexOf(node.getName());
				bs.set(index);
				map.put(node, bs);
			} else {
				for (TNode child : node.getChildren()) {
					BitSet v = (BitSet) map.get(child);
					bs.or(v);
					if (concluster.containsCluster(v)) {
						intersect++;
						virtualbs.or(v);
					}
				}

				if ((concluster.containsCluster(bs)) || (intersect > 1)) {
					if (concluster.containsCluster(bs)) {
						virtualbs = bs;
					}
					for (int i = 0; i < ((List) coveragelist).size(); i++) {
						BitSet exbs = (BitSet) ((List) coveragelist).get(i);
						BitSet temp = (BitSet) virtualbs.clone();
						temp.and(exbs);
						if (temp.equals(exbs)) {
							((List) coveragelist).remove(i);
							i--;
						}
					}
					((List) coveragelist).add(virtualbs);
				}

				map.put(node, bs);
			}

			if (!node.isRoot()) {
				BitSet complementbs = (BitSet) bs.clone();
				complementbs.flip(0, gtTaxa.length);
				if (intersect > 0) {
					complementbs.or(virtualbs);
				}
				if (concluster.containsCluster(complementbs)) {
					for (int i = 0; i < ((List) coveragelist).size(); i++) {
						BitSet exbs = (BitSet) ((List) coveragelist).get(i);
						BitSet temp = (BitSet) complementbs.clone();
						temp.and(exbs);
						if (temp.equals(exbs)) {
							((List) coveragelist).remove(i);
							i--;
						}
					}
					((List) coveragelist).add(complementbs);
					break;
				}
			}
		}
		return Math.max(0, ((List) coveragelist).size() - 1);
	}

	private static Tree networkToTree(Network net,
			Map<String, Integer> nname2tamount) {
		MutableTree tree = new STITree();
		Queue source = new LinkedList();
		Queue dest = new LinkedList();
		source.offer(net.getRoot());
		dest.offer(tree.getRoot());
		while (!source.isEmpty()) {
			NetNode parent = (NetNode) source.poll();
			TMutableNode peer = (TMutableNode) dest.poll();

			int index = 0;
			for (NetNode child : (Iterable<NetNode>)parent.getChildren()) {
				TMutableNode copy;
				if (child.getName() == "") {
					copy = peer.createChild("");
				} else {
					Integer amount = (Integer) nname2tamount.get(child
							.getName());
					if (amount == null) {
						amount = Integer.valueOf(0);
					}
					nname2tamount.put(child.getName(), amount = Integer
							.valueOf(amount.intValue() + 1));
					String newname = child.getName() + "_" + amount;
					copy = peer.createChild(newname);
				}

				double distance = child.getParentDistance(parent);
				if (distance == (-1.0D / 0.0D)) {
					copy.setParentDistance(0.0D);
				} else {
					copy.setParentDistance(distance);
				}

				source.offer(child);
				dest.offer(copy);
				index++;
			}
		}
		Trees.removeBinaryNodes(tree);
		return tree;
	}

	private static List<String[]> getOptions(String[] args) {
		LinkedList opts = new LinkedList();
		LinkedList arg_list = new LinkedList();

		int i = 2;
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
}