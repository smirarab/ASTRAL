package phylonet.coalescent;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.CharBuffer;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;

import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.util.Trees;
import phylonet.util.BitSet;

public class WQDataCollection extends AbstractDataCollection<Tripartition> implements Cloneable {

	public List<STITreeCluster> treeAllClusters = new ArrayList<STITreeCluster>();

	Tripartition[] finalTripartitions = null;
	int[] finalCounts = null;

	public Integer[] geneTreesAsInts;
	public int maxHeight;
	// private float[][] similarityMatrix;
	// private Integer[][] orderedTaxonBySimilarity;

	SimilarityMatrix similarityMatrix;
	SimilarityMatrix speciesSimilarityMatrix;

	private int n;
	private int algorithm;
	private SpeciesMapper spm;

	private boolean SLOW = false;
	private final double[] GREEDY_ADDITION_THRESHOLDS = new double[] { 0, 1 / 100., 1 / 50., 1 / 20., 1 / 10., 1 / 5.,
			1 / 3. };
	private final int GREEDY_DIST_ADDITTION_LAST_THRESHOLD_INDX = 3;
	private final int GREEDY_ADDITION_MAX_POLYTOMY_MIN = 50;
	private final int GREEDY_ADDITION_MAX_POLYTOMY_MULT = 10;
	private final int GREEDY_ADDITION_DEFAULT_RUNS = 10;
	private final double GREEDY_ADDITION_MIN_FREQ = 0.01;
	private final int GREEDY_ADDITION_IMPROVEMENT_REWARD = 2;
	private final int POLYTOMY_RESOLUTIONS = 2;
	private List<Tree> geneTrees;
	private List<Tree> completedGeeneTrees;
	private boolean outputCompleted;

	public WQDataCollection(WQClusterCollection clusters, int alg, AbstractInference<Tripartition> inference) {
		this.clusters = clusters;
		this.algorithm = alg;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
		this.SLOW = inference.getAddExtra() >= 2;
		this.geneTrees = inference.trees;
		this.completedGeeneTrees = new ArrayList<Tree>();
		this.outputCompleted = inference.shouldOutputCompleted();
	}

	private STITreeCluster getClusterForNodeName(String nodeName) {
		STITreeCluster cluster = new STITreeCluster();
		Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
		cluster.addLeaf(taxonID);
		return cluster;
	}

	private void addTripartition(STITreeCluster l_cluster, STITreeCluster r_cluster, STITreeCluster remaining,
			TNode node, Map<Tripartition, Integer> geneTreeTripartitonCount) {

		Tripartition trip = new Tripartition(l_cluster, r_cluster, remaining);
		geneTreeTripartitonCount.put(trip,
				geneTreeTripartitonCount.containsKey(trip) ? geneTreeTripartitonCount.get(trip) + 1 : 1);
	}

	void findGenetreeTripartitions(Map<Tripartition, Integer> geneTreeTripartitonCount,
			List<STITreeCluster> treeCompteleClusters) {

		System.err.println("Calculating tripartitions from gene trees ");
		int t = 0;
		for (Tree tr : this.geneTrees) {
			// System.err.print(".");
			Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
			STITreeCluster gtAll = treeCompteleClusters.get(t++);
			BitSet gtAllBS = gtAll.getBitSet();

			for (TNode node : tr.postTraverse()) {
				if (node.isLeaf()) {
					STITreeCluster cluster = getClusterForNodeName(node.getName());
					stack.add(cluster);
				} else {

					ArrayList<STITreeCluster> childbslist = new ArrayList<STITreeCluster>();
					BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					for (TNode child : node.getChildren()) {
						STITreeCluster pop = stack.pop();
						childbslist.add(pop);
						bs.or(pop.getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster((BitSet) bs.clone());
					stack.add(cluster);

					STITreeCluster remaining = cluster.complementaryCluster();
					remaining.getBitSet().and(gtAllBS);
					if (remaining.getClusterSize() != 0) {
						childbslist.add(remaining);
					}
					// System.err.println(childbslist.size());
					for (int i = 0; i < childbslist.size(); i++) {
						for (int j = i + 1; j < childbslist.size(); j++) {
							for (int k = j + 1; k < childbslist.size(); k++) {
								addTripartition(childbslist.get(i), childbslist.get(j), childbslist.get(k), node,
										geneTreeTripartitonCount);
							}
						}
					}

				}
			}

		}
	}

	void addTreeBipartitionsToX(List<Tree> trees) {

		Tree[] greedies = new Tree[POLYTOMY_RESOLUTIONS];
		for (int i = 0; i < greedies.length; i++) {
			greedies[i] = Utils.greedyConsensus(trees, true, GlobalMaps.taxonIdentifier);
			resolveByUPGMA((MutableTree) greedies[i]);
		}

		for (Tree tr : trees) {

			Stack<STITreeCluster> stack = new Stack<STITreeCluster>();

			for (TNode node : tr.postTraverse()) {
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					STITreeCluster cluster = getClusterForNodeName(node.getName());
					stack.add(cluster);
					STITreeCluster remaining = cluster.complementaryCluster();
					addARawBipartitionToX(cluster, remaining);
				} else {
					ArrayList<BitSet> childbslist = new ArrayList<BitSet>();
					BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					for (TNode child : node.getChildren()) {
						STITreeCluster pop = stack.pop();
						childbslist.add(pop.getBitSet());
						bs.or(pop.getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster(bs);
					stack.add(cluster);
					STITreeCluster remaining = cluster.complementaryCluster();

					if (addARawBipartitionToX(cluster, remaining)) {
						// if (! addTripartition) {
						// System.err.println(t+ " Extra bipartition added: " +
						// spm.getSTClusterForGeneCluster(cluster) +" |
						// "+spm.getSTClusterForGeneCluster(remaining));
						// }
					}

					/*
					 * while (childbslist.size() > 2) { STITreeCluster c1 =
					 * childbslist.remove(GlobalMaps.random.nextInt(childbslist.
					 * size())); STITreeCluster c2 =
					 * childbslist.remove(GlobalMaps.random.nextInt(childbslist.
					 * size()));
					 * 
					 * STITreeCluster newcl = new STITreeCluster(c1);
					 * newcl.getBitSet().or(c2.getBitSet()); STITreeCluster remm
					 * = newcl.complementaryCluster(); addBipartitionToX(newcl,
					 * remm);
					 * 
					 * childbslist.add(newcl); }
					 */

					if (childbslist.size() > 2) {
						// sampleAndResolve(childbslist.toArray(new BitSet[]{}),
						// false);
						boolean isRoot = remaining.getClusterSize() == 0;
						int d = childbslist.size() + (isRoot ? 0 : 1);
						BitSet[] polytomy = new BitSet[d];
						int i = 0;
						for (BitSet child : childbslist) {
							polytomy[i++] = child;
						}
						if (!isRoot) {
							polytomy[i] = remaining.getBitSet();
						}

						// TODO: do multiple samples
						HashMap<String, Integer> randomSample = this.randomSampleAroundPolytomy(polytomy);

						// STITree<Boolean> restrictedTree = new
						// STITree(greedy);
						// restrictedTree.constrainByLeaves(randomSample.keySet());
						// Utils.randomlyResolve(restrictedTree);
						for (int j = 0; j < greedies.length; j++) {
							for (BitSet restrictedBitSet : Utils.getBitsets(randomSample, greedies[j])) {
								this.addSubSampledBitSetToX(polytomy, restrictedBitSet);
							}
						}

						// System.err.print(".");
					}
				}
			}
			// System.err.println("+");

		}

	}

	Tree getCompleteTree(Tree tr, BitSet gtAllBS) {

		if (gtAllBS.cardinality() < 3) {
			throw new RuntimeException("Tree " + tr.toNewick() + " has less than 3 taxa; it cannot be completed");
		}
		STITree trc = new STITree(tr);

		Trees.removeBinaryNodes(trc);

		for (int missingId = gtAllBS.nextClearBit(0); missingId < n; missingId = gtAllBS.nextClearBit(missingId + 1)) {

			int closestId = similarityMatrix.getClosestPresentTaxonId(gtAllBS, missingId);

			STINode closestNode = trc.getNode(GlobalMaps.taxonIdentifier.getTaxonName(closestId));

			trc.rerootTreeAtNode(closestNode);
			Trees.removeBinaryNodes(trc);

			Iterator cit = trc.getRoot().getChildren().iterator();
			STINode c1 = (STINode) cit.next();
			STINode c2 = (STINode) cit.next();
			STINode start = closestNode == c1 ? c2 : c1;

			int c1random = -1;
			int c2random = -1;
			while (true) {
				if (start.isLeaf()) {
					break;
				}

				cit = start.getChildren().iterator();
				c1 = (STINode) cit.next();
				c2 = (STINode) cit.next();

				// TODO: what if c1 or c2 never appears in the same tree as
				// missing and closestId .
				if (c1random == -1) {
					c1random = GlobalMaps.taxonIdentifier.taxonId(Utils.getLeftmostLeaf(c1));
				}
				if (c2random == -1) {
					c2random = GlobalMaps.taxonIdentifier.taxonId(Utils.getLeftmostLeaf(c2));
				}
				int betterSide = similarityMatrix.getBetterSideByFourPoint(missingId, closestId, c1random, c2random);
				if (betterSide == closestId) {
					break;
				} else if (betterSide == c1random) {
					start = c1;
					// Currently, c1random is always under left side of c1
					c2random = -1;
				} else if (betterSide == c2random) {
					start = c2;
					// Currently, c2random is always under left side of c2
					c1random = c2random;
					c2random = -1;
				}

			}
			if (start.isLeaf()) {
				STINode newnode = start.getParent().createChild(GlobalMaps.taxonIdentifier.getTaxonName(missingId));
				STINode newinternalnode = start.getParent().createChild();
				newinternalnode.adoptChild(start);
				newinternalnode.adoptChild(newnode);
			} else {
				STINode newnode = start.createChild(GlobalMaps.taxonIdentifier.getTaxonName(missingId));
				STINode newinternalnode = start.createChild();
				newinternalnode.adoptChild(c1);
				newinternalnode.adoptChild(c2);
			}
		}

		return trc;
	}

	private boolean addARawBitSetToX(BitSet bs) {
		STITreeCluster cluster = new STITreeCluster();
		cluster.setCluster(bs);
		return this.addARawBipartitionToX(cluster, cluster.complementaryCluster());
	}

	/*
	 * Adds bipartitions to X, ensuring that individuals from the same species
	 * are not on both sides of any bipartitions Uses heuristics to move
	 * individuals around if they are on both sides.
	 */
	private boolean addARawBipartitionToX(STITreeCluster c1, STITreeCluster c2) {

		boolean added = false;

		STITreeCluster c1copy = new STITreeCluster(c1);
		STITreeCluster c2copy = new STITreeCluster(c2);
		BitSet b1copy = c1copy.getBitSet();
		BitSet b2copy = c2copy.getBitSet();

		/*
		 * Find out for each species whether they are more frequent in left or
		 * right
		 */
		int[] countsC1c = new int[spm.getSpeciesCount()], countsC2c = new int[spm.getSpeciesCount()];
		for (int i = b1copy.nextSetBit(0); i >= 0; i = b1copy.nextSetBit(i + 1)) {
			int sID = spm.getSpeciesIdForTaxon(i);
			countsC1c[sID] += 10;
			if (spm.getLowestIndexIndividual(sID) == i) {
				countsC1c[sID]++;
			}
		}
		for (int i = b2copy.nextSetBit(0); i >= 0; i = b2copy.nextSetBit(i + 1)) {
			int sID = spm.getSpeciesIdForTaxon(i);
			countsC2c[sID] += 10;
			if (spm.getLowestIndexIndividual(sID) == i) {
				countsC2c[sID]++;
			}
		}
		/*
		 * Add a bipartition where every individual is moved to the side where
		 * it is more common
		 */
		BitSet bs1Voted = new BitSet(spm.getSpeciesCount());
		for (int i = 0; i < countsC2c.length; i++) {
			if (countsC1c[i] > countsC2c[i]) {
				bs1Voted.set(i);
			}
		}
		STITreeCluster c1Voted = spm.getGeneClusterForSTCluster(bs1Voted);
		added |= this.addCompletedSpeciesFixedBipartionToX(c1Voted, c1Voted.complementaryCluster());

		/*
		 * Add two more bipartitions by adding all individuals from each species
		 * to each side they appear at least once
		 */
		// spm.addMissingIndividuals(c1copy.getBitSet());
		// spm.addMissingIndividuals(c2copy.getBitSet());
		// STITreeCluster c1cComp = c1copy.complementaryCluster();
		// STITreeCluster c2cComp = c2copy.complementaryCluster();
		// added |= this.addCompletedSpeciesFixedBipartionToX(c1copy, c1cComp);
		// added |= this.addCompletedSpeciesFixedBipartionToX(c2copy, c2cComp);

		return added;

	}

	public void addExtraBipartitionsByInput(List<Tree> extraTrees, boolean extraTreeRooted) {

		List<Tree> completedExtraGeeneTrees = new ArrayList<Tree>();
		for (Tree tr : extraTrees) {
			String[] gtLeaves = tr.getLeaves();
			STITreeCluster gtAll = new STITreeCluster();
			for (int i = 0; i < gtLeaves.length; i++) {
				gtAll.addLeaf(GlobalMaps.taxonIdentifier.taxonId(gtLeaves[i]));
			}
			Tree trc = getCompleteTree(tr, gtAll.getBitSet());

			completedExtraGeeneTrees.add(trc);
		}
		addTreeBipartitionsToX(completedExtraGeeneTrees);
	}

	private boolean addCompletedSpeciesFixedBipartionToX(STITreeCluster c1, STITreeCluster c2) {
		boolean added = false;
		int size = c1.getClusterSize();
		if (size == n || size == 0) {
			return false;
		}
		added |= addToClusters(c1, size);
		size = c2.getClusterSize();
		added |= addToClusters(c2, size);
		return added;
	}

	public void computeTreePartitions(AbstractInference<Tripartition> inference) {

		int haveMissing = preProcess(inference);

		calculateDistances();

		if (haveMissing > 0) {
			completeGeneTrees();
		} else {
			this.completedGeeneTrees = this.geneTrees;
		}

		/*
		 * Calculate gene tree clusters and bipartitions for X
		 */
		STITreeCluster all = new STITreeCluster();
		all.getBitSet().set(0, n);
		addToClusters(all, GlobalMaps.taxonIdentifier.taxonCount());
		System.err.println("Building set of clusters (X) from gene trees ");
		addTreeBipartitionsToX(this.completedGeeneTrees);

		initializeWeightCalculator(inference);

		long maxpos;

		if (inference instanceof AbstractInferenceNoCalculations) {
			((WQInferenceNoCalculations) inference).maxpossible = maxpos = this.calculateMaxPossible();
		} else
			((WQInference) inference).maxpossible = maxpos = this.calculateMaxPossible();
		System.err.println("Number of quartet trees in the gene trees: " + maxpos);

	}

	void initializeWeightCalculator(AbstractInference<Tripartition> inference) {
		/*
		 * If needed calculate gene tree tripartitions
		 */
		int k = this.geneTrees.size();
		Map<Tripartition, Integer> geneTreeTripartitonCount = new HashMap<Tripartition, Integer>(k * n);
		if (this.algorithm == 2 || this.algorithm == -1) {
			findGenetreeTripartitions(geneTreeTripartitonCount, treeAllClusters);
		}
		if (this.algorithm == -1) {
			this.setAlgorithm(geneTreeTripartitonCount.size(), k);
		}

		if (this.algorithm == 2) {
			System.err.println("Using tripartition-based weight calculation.");

			finalTripartitions = new Tripartition[geneTreeTripartitonCount.size()];
			finalCounts = new int[geneTreeTripartitonCount.size()];
			int i = 0;
			for (Entry<Tripartition, Integer> entry : geneTreeTripartitonCount.entrySet()) {
				finalTripartitions[i] = entry.getKey();
				finalCounts[i] = entry.getValue();
				i++;
			}
		} else {
			System.err.println("Using tree-based weight calculation.");
			List<Integer> temp = new ArrayList<Integer>();
			Stack<Integer> stack = new Stack<Integer>();
			maxHeight = 0;
			for (Tree tr : this.geneTrees) {
				for (TNode node : tr.postTraverse()) {
					if (node.isLeaf()) {
						temp.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
						stack.push(0);
					} else {
						temp.add(-node.getChildCount());
						int h = 0;
						for (int i = 0; i < node.getChildCount(); i++) {
							int childheight = stack.pop();
							if(childheight > h)
								h = childheight;
						}
						h++;
						stack.push(h);
					}
					if (node.isRoot()) {
						temp.add(Integer.MIN_VALUE);
						stack.clear();
					}
					if(stack.size()>maxHeight) {
						maxHeight = stack.size();
					}
				}
			}

			geneTreesAsInts = temp.toArray(new Integer[] {});
			if (WQWeightCalculator.WRITE_OR_DEBUG) {
				FileWriter writer = null;
				try {
					writer = new FileWriter("geneTreesAsInts.txt");

					for (int i = 0; i < geneTreesAsInts.length; i++) {
						writer.write(geneTreesAsInts[i].toString() + " ");
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} finally {
					if (writer != null) {
						try {
							writer.close();
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}
			}
		}

		if (geneTreeTripartitonCount.size() > 0) {
			long s = 0;
			for (Integer c : geneTreeTripartitonCount.values()) {
				s += c;
			}
			System.err.println("Tripartitons in gene trees (count): " + geneTreeTripartitonCount.size());
			System.err.println("Tripartitons in gene trees (sum): " + s);
		}
		System.err.println("Number of Default Clusters: " + clusters.getClusterCount());

		inference.weightCalculator.initializeWeightContainer(geneTreeTripartitonCount.size() * 2);
	}

	long calculateMaxPossible() {
		long weight = 0;
		Integer allsides = null;
		Iterator<STITreeCluster> tit = this.treeAllClusters.iterator();
		boolean newTree = true;

		Deque<Integer> stack = new ArrayDeque<Integer>();
		for (Integer gtb : this.geneTreesAsInts) {
			if (newTree) {
				allsides = tit.next().getBitSet().cardinality();
				newTree = false;
			}
			if (gtb >= 0) {
				stack.push(1);
			} else if (gtb == Integer.MIN_VALUE) {
				stack.clear();
				newTree = true;
			} else {
				ArrayList<Integer> children = new ArrayList<Integer>();
				Integer newSide = 0;
				for (int i = gtb; i < 0; i++) {
					Integer pop = stack.pop();
					children.add(pop);
					newSide += pop;
				}
				stack.push(newSide);
				Integer sideRemaining = allsides - newSide;
				if (sideRemaining != 0) {
					children.add(sideRemaining);
				}
				for (int i = 0; i < children.size(); i++) {
					Long a = children.get(i) + 0l;

					for (int j = i + 1; j < children.size(); j++) {
						Long b = children.get(j) + 0l;
						/*
						 * if (children.size() > 5) { if ((side1.s0+side2.s0 ==
						 * 0? 1 :0) + (side1.s1+side2.s1 == 0? 1 :0) +
						 * (side1.s2+side2.s2 == 0? 1:0) > 1) continue; }
						 */
						for (int k = j + 1; k < children.size(); k++) {
							Long c = children.get(k) + 0l;
							weight += (a + b + c - 3) * a * b * c;
						}
					}
				}
			}
		}
		return weight / 4l;
	}

	private void calculateDistances() {
		System.err.println("Calculating quartet distance matrix (for completion of X)");

		this.similarityMatrix = new SimilarityMatrix(n);
		this.similarityMatrix.populateByQuartetDistance(treeAllClusters, this.geneTrees);
		this.speciesSimilarityMatrix = this.similarityMatrix.convertToSpeciesDistance(spm);
	}

	int preProcess(AbstractInference<Tripartition> inference) {
		System.err.println("Number of gene trees: " + this.geneTrees.size());
		n = GlobalMaps.taxonIdentifier.taxonCount();

		int haveMissing = 0;
		for (Tree tree : this.geneTrees) {
			if (tree.getLeafCount() != n) {
				haveMissing++;
			}
			String[] gtLeaves = tree.getLeaves();
			STITreeCluster gtAll = new STITreeCluster();
			long ni = gtLeaves.length;
			for (int i = 0; i < ni; i++) {
				gtAll.addLeaf(GlobalMaps.taxonIdentifier.taxonId(gtLeaves[i]));
			}
			treeAllClusters.add(gtAll);
		}
		System.err.println(haveMissing + " trees have missing taxa");

		return haveMissing;
	}

	long maxPossibleScore(Tripartition trip) {

		long weight = 0;

		for (STITreeCluster all : this.treeAllClusters) {
			long a = trip.cluster1.getBitSet().intersectionSize(all.getBitSet()),
					b = trip.cluster2.getBitSet().intersectionSize(all.getBitSet()),
					c = trip.cluster3.getBitSet().intersectionSize(all.getBitSet());

			weight += (a + b + c - 3) * a * b * c;
		}
		return weight;
	}

	private void completeGeneTrees() {
		System.err.println("Will attempt to complete bipartitions from X before adding using a distance matrix.");
		int t = 0;
		BufferedWriter completedFile = null;
		if (this.outputCompleted) {
			String fn = GlobalMaps.outputfilename + ".completed_gene_trees";
			System.err.println("Outputting completed gene trees to " + fn);
			try {
				completedFile = new BufferedWriter(new FileWriter(fn));
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		for (Tree tr : this.geneTrees) {
			Tree trc = getCompleteTree(tr, this.treeAllClusters.get(t++).getBitSet());
			this.completedGeeneTrees.add(trc);
			if (completedFile != null) {
				try {
					completedFile.write(trc.toStringWD() + " \n");
					completedFile.flush();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		}
		if (completedFile != null) {
			try {
				completedFile.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}

	private void printoutdistmatrix(double[][] distSTMatrix) {
		for (String s : spm.getSTTaxonIdentifier().getAllTaxonNames()) {
			System.err.print(String.format("%1$8s", s));
		}
		System.err.println();
		for (int i = 0; i < spm.getSpeciesCount(); i++) {
			for (int j = 0; j < spm.getSpeciesCount(); j++) {
				System.err.print(String.format("%1$8.3f", distSTMatrix[i][j]));
			}
			System.err.println();
		}
	}

	private void setAlgorithm(int geneTreeTripartitonCountSize, int k) {
		if (this.algorithm != -1) {
			return;
		}
		if (k <= 0 || geneTreeTripartitonCountSize <= 0) {
			throw new RuntimeException("gene tree tripartition size or k not set properly");
		}
		if (this.algorithm == -1) {
			this.algorithm = (n <= 32 || (geneTreeTripartitonCountSize < k * 6)) ? 2 : 1;
		} else {
			throw new RuntimeException("Algorithm already set");
		}
	}

	public void addExtraBipartitionByDistance() {

		for (BitSet bs : speciesSimilarityMatrix.UPGMA()) {
			STITreeCluster g = spm.getGeneClusterForSTCluster(bs);
			this.addCompletedSpeciesFixedBipartionToX(g, g.complementaryCluster());
			// upgmac.add(g);
		}
		;
		if (SLOW) {
			for (BitSet bs : speciesSimilarityMatrix.getQuadraticBitsets()) {
				STITreeCluster g = spm.getGeneClusterForSTCluster(bs);
				this.addCompletedSpeciesFixedBipartionToX(g, g.complementaryCluster());
			}
			;
		}

		System.err.println("Number of Clusters after addition by distance: " + clusters.getClusterCount());
	}

	@Override
	public void addExtraBipartitionByExtension(AbstractInference<Tripartition> inference) {

		this.addExtraBipartitionByDistance();

		System.err.println("Adding to X using resolutions of greedy consensus ...");
		for (Tree tree : this.completedGeeneTrees) {
			tree.rerootTreeAtEdge(GlobalMaps.taxonIdentifier.getTaxonName(0));
			Trees.removeBinaryNodes((MutableTree) tree);
		}

		/*
		 * if (completeTrees.size() < 2) { System.err.println("Only "
		 * +completeTrees.size() +
		 * " complete trees found. Greedy-based completion not applicable.");
		 * return; }
		 */

		Collection<Tree> allGreedies = Utils.greedyConsensus(this.completedGeeneTrees, GREEDY_ADDITION_THRESHOLDS, true,
				1, GlobalMaps.taxonIdentifier);

		int th = 0;
		for (Tree cons : allGreedies) {
			double thresh = GREEDY_ADDITION_THRESHOLDS[th];
			System.err.println("Threshold " + thresh + ":");
			th = (th + 1) % GREEDY_ADDITION_THRESHOLDS.length;
			// System.err.println(cons);
			Stack<BitSet> greedyNodeStack = new Stack<BitSet>();
			for (TNode greedyNode : cons.postTraverse()) {
				if (greedyNode.isLeaf()) {
					BitSet nbs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					nbs.set(GlobalMaps.taxonIdentifier.taxonId(greedyNode.getName()));
					greedyNodeStack.push(nbs);
					continue;
				}

				BitSet[] childbs = new BitSet[greedyNode.getChildCount() + 1];

				BitSet greedyBS = new BitSet();
				for (int i = 0; i < greedyNode.getChildCount(); i++) {
					BitSet pop = greedyNodeStack.pop();
					greedyBS.or(pop);
					childbs[i] = pop;
				}
				greedyNodeStack.push(greedyBS);

				if (greedyNode.getChildCount() > 2 && greedyNode.getChildCount() < this.getRoundCount()) {

					BitSet comp = (BitSet) greedyBS.clone();
					comp.flip(0, n);
					childbs[greedyNode.getChildCount()] = comp;

					System.err.print("polytomy of size " + greedyNode.getChildCount());

					int k = 0;

					this.resolveByUPGMA(childbs);

					for (int j = 0; j < GREEDY_ADDITION_DEFAULT_RUNS + k; j++) {

						boolean quadratic = SLOW || (th <= GREEDY_DIST_ADDITTION_LAST_THRESHOLD_INDX
								&& j < GREEDY_ADDITION_DEFAULT_RUNS);

						if (sampleAndResolve(childbs, quadratic)) {
							k += GREEDY_ADDITION_IMPROVEMENT_REWARD;
						}
					}
					System.err.println("; rounds with additions with at least " + GREEDY_ADDITION_MIN_FREQ
							+ " support: " + k / GREEDY_ADDITION_IMPROVEMENT_REWARD + "; clusters: "
							+ clusters.getClusterCount());
				}
			}
		}
		System.err.println("Number of Clusters after addition by greedy: " + clusters.getClusterCount());
	}

	private long getRoundCount() {
		return this.GREEDY_ADDITION_MAX_POLYTOMY_MIN
				+ Math.round(Math.sqrt(n * this.GREEDY_ADDITION_MAX_POLYTOMY_MULT));
	}

	private boolean resolveByUPGMA(BitSet[] polytomyBSList) {
		boolean added = false;

		for (BitSet bs : this.similarityMatrix.resolveByUPGMA(Arrays.asList(polytomyBSList), true)) {
			added |= this.addARawBitSetToX(bs);
		}
		return added;
	}

	private HashMap<BitSet, Integer> returnBitSetCounts(List<Tree> genetrees, HashMap<String, Integer> randomSample) {

		HashMap<BitSet, Integer> counts = new HashMap<BitSet, Integer>();

		for (Tree gt : genetrees) {
			List<BitSet> bsList = Utils.getBitsets(randomSample, gt);

			for (BitSet bs : bsList) {
				if (counts.containsKey(bs)) {
					counts.put(bs, counts.get(bs) + 1);
					continue;
				}
				BitSet bs2 = (BitSet) bs.clone();
				bs2.flip(0, randomSample.size());
				if (counts.containsKey(bs2)) {
					counts.put(bs2, counts.get(bs2) + 1);
					continue;
				}
				counts.put(bs, 1);
			}
		}
		return counts;
	}

	private boolean sampleAndResolve(BitSet[] polytomyBSList, boolean addQuadratic) {

		boolean addedHighFreq = false;
		// random sample taxa
		HashMap<String, Integer> randomSample = randomSampleAroundPolytomy(polytomyBSList);

		addedHighFreq = resolveLinearly(polytomyBSList, randomSample);

		resolveByDistance(polytomyBSList, randomSample, addQuadratic);

		return addedHighFreq;
	}

	private boolean resolveLinearly(BitSet[] polytomyBSList, HashMap<String, Integer> randomSample) {
		int sampleSize = randomSample.size();
		// get bipartition counts in the induced trees
		HashMap<BitSet, Integer> counts = returnBitSetCounts(this.completedGeeneTrees, randomSample);

		// sort bipartitions
		TreeSet<Entry<BitSet, Integer>> countSorted = new TreeSet<Entry<BitSet, Integer>>(
				new Utils.BSComparator(true, sampleSize));
		countSorted.addAll(counts.entrySet());

		// build the greedy tree
		MutableTree greedyTree = new STITree<BitSet>();
		TNode[] tmpnodes = new TNode[sampleSize];
		for (int i = 0; i < sampleSize; i++) {
			tmpnodes[i] = greedyTree.getRoot().createChild(i + "");
			BitSet bs = new BitSet(sampleSize);
			bs.set(i);
			((STINode<BitSet>) tmpnodes[i]).setData(bs);
		}

		boolean added = false;
		boolean addedHighFreq = false;
		// System.err.print("^");
		for (Entry<BitSet, Integer> entry : countSorted) {

			BitSet newbs = entry.getKey();

			SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(greedyTree);
			Set<TNode> clusterLeaves = new HashSet<TNode>();
			TNode node;
			for (int i = newbs.nextSetBit(0); i >= 0; i = newbs.nextSetBit(i + 1)) {
				node = tmpnodes[i];
				clusterLeaves.add(node);
			}
			TNode lca = lcaFinder.getLCA(clusterLeaves);
			LinkedList<TNode> movedChildren = new LinkedList<TNode>();
			int nodes = clusterLeaves.size();
			for (TNode child : lca.getChildren()) {
				BitSet childCluster = ((STINode<BitSet>) child).getData();

				BitSet temp = (BitSet) childCluster.clone();
				temp.and(newbs);
				if (temp.equals(childCluster)) {
					movedChildren.add(child);
					nodes -= temp.cardinality();
				}

			}

			// boolean isPartOfGreedy = false;
			if (movedChildren.size() != 0 && nodes == 0) {

				STINode<BitSet> newChild = ((STINode<BitSet>) lca).createChild();
				newChild.setData(newbs);

				while (!movedChildren.isEmpty()) {
					newChild.adoptChild((TMutableNode) movedChildren.get(0));
					movedChildren.remove(0);
				}

				if (addSubSampledBitSetToX(polytomyBSList, newbs)) {
					if (GREEDY_ADDITION_MIN_FREQ <= (entry.getValue() + 0.0) / this.completedGeeneTrees.size()) {
						addedHighFreq = true;
					}
					added = true;
				}
			}

		}

		if (added) {
			for (TNode node : greedyTree.postTraverse()) {
				if (node.getChildCount() < 3) {
					continue;
				}
				ArrayList<BitSet> children = new ArrayList<BitSet>(node.getChildCount() + 1);
				BitSet rest = new BitSet(sampleSize);
				for (TNode child : node.getChildren()) {
					children.add(((STINode<BitSet>) child).getData());
					rest.or(((STINode<BitSet>) child).getData());
				}
				rest.flip(0, sampleSize);
				if (rest.cardinality() != 0)
					children.add(rest);

				for (BitSet bs : this.similarityMatrix.resolveByUPGMA(children, true)) {
					addSubSampledBitSetToX(polytomyBSList, bs);
				}

				while (children.size() > 2) {
					BitSet c1 = children.remove(GlobalMaps.random.nextInt(children.size()));
					BitSet c2 = children.remove(GlobalMaps.random.nextInt(children.size()));

					BitSet newbs = (BitSet) c1.clone();
					newbs.or(c2);
					addSubSampledBitSetToX(polytomyBSList, newbs);
					children.add(newbs);
				}
			}
		}
		return addedHighFreq;
	}

	private boolean resolveByDistance(BitSet[] polytomyBSList, HashMap<String, Integer> randomSample,
			boolean quartetAddition) {
		boolean added = false;

		SimilarityMatrix sampleSimMatrix = this.similarityMatrix.getInducedMatrix(randomSample);

		for (BitSet bs : sampleSimMatrix.UPGMA()) {
			added |= this.addSubSampledBitSetToX(polytomyBSList, bs);
		}
		if (quartetAddition) {
			for (BitSet bs : sampleSimMatrix.getQuadraticBitsets()) {
				added |= this.addSubSampledBitSetToX(polytomyBSList, bs);
			}
		}
		return added;
	}

	private HashMap<String, Integer> randomSampleAroundPolytomy(BitSet[] polyTomy) {
		HashMap<String, Integer> randomSample = new HashMap<String, Integer>();
		int ind = 0;
		for (BitSet child : polyTomy) {
			int sample = GlobalMaps.random.nextInt(child.cardinality());
			int p = child.nextSetBit(0);
			for (int i = 0; i < sample; i++) {
				p = child.nextSetBit(p + 1);
			}
			randomSample.put(GlobalMaps.taxonIdentifier.getTaxonName(p), ind);
			ind++;
		}
		return randomSample;
	}

	private boolean addSubSampledBitSetToX(BitSet[] childbs, BitSet restrictedBitSet) {
		BitSet stNewbs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
		for (int j = restrictedBitSet.nextSetBit(0); j >= 0; j = restrictedBitSet.nextSetBit(j + 1)) {
			stNewbs.or(childbs[j]);
		}

		return this.addARawBitSetToX(stNewbs);
	}

	private void resolveByUPGMA(MutableTree tree) {
		Stack<BitSet> stack = new Stack<BitSet>();
		for (TNode node : tree.postTraverse()) {
			BitSet bitset = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
			if (node.isLeaf()) {
				bitset.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
			} else {
				List<TMutableNode> children = new ArrayList<TMutableNode>();
				ArrayList<BitSet> poly = new ArrayList<BitSet>();
				for (TNode child : node.getChildren()) {
					BitSet cbs = stack.pop();
					poly.add(cbs);
					children.add((TMutableNode) child);
					bitset.or(cbs);
				}
				if (children.size() > 2) {

					for (BitSet bs : this.similarityMatrix.resolveByUPGMA(poly, false)) {
						TMutableNode newChild = ((TMutableNode) node).createChild();
						for (int i = bs.nextSetBit(0); i >= 0; i = bs.nextSetBit(i + 1)) {
							TMutableNode child = children.get(i);
							if (child.getParent() == node) {
								newChild.adoptChild(child);
							}
							children.set(i, newChild);
						}
					}
				}
			}
			stack.push(bitset);
		}
	}

	public Object clone() throws CloneNotSupportedException {
		WQDataCollection clone = (WQDataCollection) super.clone();
		clone.clusters = (WQClusterCollection) ((AbstractClusterCollection) this.clusters).clone();
		return clone;
	}
}
