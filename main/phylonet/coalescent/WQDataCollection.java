package phylonet.coalescent;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.Stack;

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

public class WQDataCollection extends AbstractDataCollection<Tripartition> {


	List<STITreeCluster> treeAllClusters = new ArrayList<STITreeCluster>();

	Tripartition [] finalTripartitions = null;
	int [] finalCounts = null;	

	Integer [] geneTreesAsInts;

	private float[][] similarityMatrix;
	private Integer[][] orderedTaxonBySimilarity;
	
	private int n;
	private int algorithm;
	private SpeciesMapper spm;

	private boolean DISTANCE_ADDITION = false;
	private final double [] GREEDY_ADDITION_THRESHOLDS = new double [] {0, 1/100., 1/50., 1/20., 1/10., 1/5., 1/3.} ;
	private final double GREEDY_ADDITION_MIN_FREQ = 0.01;
	private final int GREEDY_ADDITION_NOIMPROVEMENT_LIMIT = 20;
	//private final int P = 100;
	//private final int M = 1000;
	private List<Tree> geneTrees;
	private List<Tree> completedGeeneTrees;
	private boolean outputCompleted;

	public WQDataCollection( WQClusterCollection clusters, int alg, 
			AbstractInference<Tripartition> inference) {
		this.clusters = clusters;
		this.algorithm = alg;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
		this.DISTANCE_ADDITION = inference.getAddExtra() >= 2;
		this.geneTrees = inference.trees;
		this.completedGeeneTrees = new ArrayList<Tree>();
		this.outputCompleted = inference.outputCompleted;
	}


	private STITreeCluster getClusterForNodeName(String nodeName) {
		STITreeCluster cluster = new STITreeCluster();
		Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
		cluster.addLeaf(taxonID);
		return cluster;
	}


	private void addTripartition(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster remaining, TNode node,
			Map<Tripartition, Integer> geneTreeTripartitonCount) {
	
		Tripartition trip = new Tripartition(l_cluster, r_cluster, remaining);
		geneTreeTripartitonCount.put(trip,
				geneTreeTripartitonCount.containsKey(trip) ? 
						geneTreeTripartitonCount.get(trip) + 1 : 1);
	}


	void findGenetreeTripartitions( Map<Tripartition, Integer> geneTreeTripartitonCount,
			List<STITreeCluster> treeCompteleClusters) {
	
		System.err.println("Calculating tripartitions from gene trees ");
		int t = 0;
		for (Tree tr : this.geneTrees) {
			//System.err.print(".");
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
					for (TNode child: node.getChildren()) {
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
					//System.err.println(childbslist.size());
					for (int i = 0; i < childbslist.size(); i++) {
						for (int j = i+1; j < childbslist.size(); j++) {
							for (int k = j+1; k < childbslist.size(); k++) {
								addTripartition( childbslist.get(i),  childbslist.get(j), 
										childbslist.get(k), node, geneTreeTripartitonCount);
							}
						}					       
					}
					
				}
			}
	
		}
	}


	void addTreeBipartitionsToX(List<Tree> trees) {

		Tree[] greedies = new Tree[1];
		for (int i = 0; i < greedies.length; i++ ) {
			greedies[i] = Utils.greedyConsensus(trees, true);
			Utils.randomlyResolve((MutableTree) greedies[i]);
		}
		
		for (Tree tr : trees) {
			
			Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
			
			for (TNode node : tr.postTraverse()) {				
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					STITreeCluster cluster = getClusterForNodeName(node.getName());
					stack.add(cluster);
					STITreeCluster remaining = cluster.complementaryCluster();
					addBipartitionToX(cluster,remaining);
				} else {
					ArrayList<BitSet> childbslist = new ArrayList<BitSet>();
					BitSet bs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					for (TNode child: node.getChildren()) {
						STITreeCluster pop = stack.pop();
						childbslist.add(pop.getBitSet());
						bs.or(pop.getBitSet());
					}

					STITreeCluster cluster = new STITreeCluster();
					cluster.setCluster(bs);
					stack.add(cluster);
					STITreeCluster remaining = cluster.complementaryCluster();

					if (addBipartitionToX(cluster, remaining)) {
						//if (! addTripartition) {
							//System.err.println(t+ " Extra bipartition added: " + spm.getSTClusterForGeneCluster(cluster) +" | "+spm.getSTClusterForGeneCluster(remaining));
						//}
					}
					
					/*while (childbslist.size() > 2) {						
						STITreeCluster c1 = childbslist.remove(GlobalMaps.random.nextInt(childbslist.size()));
						STITreeCluster c2 = childbslist.remove(GlobalMaps.random.nextInt(childbslist.size()));
						
						STITreeCluster newcl = new STITreeCluster(c1);
						newcl.getBitSet().or(c2.getBitSet());						
						STITreeCluster remm = newcl.complementaryCluster();
						addBipartitionToX(newcl, remm);
						
						childbslist.add(newcl);
					}*/
				
					if (childbslist.size() > 2) {	
						//sampleAndResolve(childbslist.toArray(new BitSet[]{}), false);
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
						
						HashMap<String, Integer> randomSample = this.randomSampleAroundPolytomy(polytomy);
						
						//STITree<Boolean> restrictedTree = new STITree(greedy);
						//restrictedTree.constrainByLeaves(randomSample.keySet());
						//Utils.randomlyResolve(restrictedTree);
						for (int j = 0; j < greedies.length; j++) {
							for (BitSet restrictedBitSet :  Utils.getBitsets(randomSample, greedies[j])) {
								this.addSubSampledBitSetToX(polytomy, restrictedBitSet);
							}
						}
						
						//System.err.print(".");
					}
				}
			}
			//System.err.println("+");

		}

	}


	Tree getCompleteTree(Tree tr, BitSet gtAllBS) {
		STITree trc = new STITree(tr);		
		
		for (int missingId = gtAllBS.nextClearBit(0); missingId < n ; missingId = gtAllBS.nextClearBit(missingId+1)) {
			
			int closestId = -1;
			for (int j = 0; ; j++){
				if ( missingId > this.orderedTaxonBySimilarity[missingId][j] // already added
						|| gtAllBS.get(this.orderedTaxonBySimilarity[missingId][j]) // was in original tree
						) {
					closestId = this.orderedTaxonBySimilarity[missingId][j];
					break;
				}
				if (j > GlobalMaps.taxonIdentifier.taxonCount()) {
					throw new RuntimeException("Bug: this should not be reached");
				}
			}
			
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
				
				if (c1random == -1) {
					c1random = GlobalMaps.taxonIdentifier.taxonId(Utils.getLeftmostLeaf(c1));
				}
				if (c2random == -1) {
					c2random = GlobalMaps.taxonIdentifier.taxonId(Utils.getLeftmostLeaf(c2));
				}
				int betterSide = getBetterSideByFourPoint(missingId, closestId, c1random, c2random);
				if (betterSide == closestId) {
					break;
				} else if (betterSide == c1random) { 
					start = c1;
					//Currently, c1random is always under left side of c1
					c2random = -1;
				} else if (betterSide == c2random) {
					start = c2;
					//Currently, c2random is always under left side of c2
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
		
		//System.err.println(trc);
		
		return trc;
	}

	int getBetterSideByFourPoint(int x, int a, int b, int c) {
		double xa = this.similarityMatrix[x][a];
		double xb = this.similarityMatrix[x][b];
		double xc = this.similarityMatrix[x][c];
		double ab = this.similarityMatrix[a][b];
		double ac = this.similarityMatrix[a][c];
		double bc = this.similarityMatrix[b][c];
		double ascore = xa + bc  - (xb + ac); // Note this is similartiy, not distance
		double bscore = xb + ac  - (xa + bc); 
		double cscore = xc + ab - (xb + ac); 
		return ascore >= bscore ?
				ascore >= cscore ? a : c :
					bscore >= cscore ? b : c;	
	}

	private boolean addBipartitionToX(STITreeCluster c1, STITreeCluster c2) {

		boolean added = false;

		// TODO: should this be treated differently?
		if (c1.getClusterSize() == 1) {
			added |= addToClusters(c1, c1.getClusterSize());
		}

		STITreeCluster c1c = new STITreeCluster (c1);
		STITreeCluster c2c = new STITreeCluster (c2);
		
		BitSet b1c = c1c.getBitSet();
		BitSet b2c = c2c.getBitSet();
				
		int [] countsC1c = new int [spm.getSpeciesCount()], countsC2c = new int [spm.getSpeciesCount()];
		int s1 = 0, s2 = 0;
		for (int i = b1c.nextSetBit(0); i >=0 ; i = b1c.nextSetBit(i+1)) {
			countsC1c[spm.getSpeciesIdForTaxon(i)]++;  
			s1++;
		}
		for (int i = b2c.nextSetBit(0); i >=0 ; i = b2c.nextSetBit(i+1)) {
			countsC2c[spm.getSpeciesIdForTaxon(i)]++;   
			s2++;
		}  
		BitSet bs1 = new BitSet(spm.getSpeciesCount()); 
		for (int i = 0; i < countsC2c.length; i++) {
			if (countsC1c[i] > countsC2c[i] || ((countsC1c[i] == countsC2c[i]) && (s1 < s2))) {
				bs1.set(i);
			} 
		}
		STITreeCluster c1s = spm.getGeneClusterForSTCluster(bs1);
		added |=  this.addCompletedBipartionToX(c1s, c1s.complementaryCluster());


		spm.addMissingIndividuals(c1c.getBitSet());
		spm.addMissingIndividuals(c2c.getBitSet());

		added |= this.addCompletedBipartionToX(c1c, c2c);

		c1c = c1c.complementaryCluster();
		c2c = c2c.complementaryCluster();

		added |= this.addCompletedBipartionToX(c1c, c2c);

		return added;

	}


	public void addExtraBipartitionsByInput(
			List<Tree> extraTrees, boolean extraTreeRooted) {
		
		List<Tree> completedExtraGeeneTrees = new ArrayList<Tree>();
		for (Tree tr: extraTrees) {
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


	private boolean addCompletedBipartionToX(STITreeCluster c1, STITreeCluster c2) {
		boolean added = false;
		int size = c1.getClusterSize();
		if (size == n || size == 0) {
			return false;
		}
		added |= addToClusters(c1, size);  
		size  = c2.getClusterSize();
		added |= addToClusters(c2, size);
		return added;
	}


	public void computeTreePartitions(AbstractInference<Tripartition> inference) {

		int k = this.geneTrees.size();
		System.err.println("Number of gene trees: " + k);
		
		n = GlobalMaps.taxonIdentifier.taxonCount();

		int haveMissing = 0;
		
		for (Tree tree :  this.geneTrees) {
			if (tree.getLeafCount() != n) {
				haveMissing++;
			}
			String[] gtLeaves = tree.getLeaves();
			STITreeCluster gtAll = new STITreeCluster();
			long ni = gtLeaves.length;
			for (int i = 0; i < ni ; i++) {
				gtAll.addLeaf(GlobalMaps.taxonIdentifier.taxonId(gtLeaves[i]));
			}
			treeAllClusters.add(gtAll);
			((WQInference)inference).maxpossible += ni * (ni -1) * (ni - 2) * (ni -3) / 24l;
		}
		
		System.err.println( haveMissing + " trees have missing taxa");
		
		if ( (haveMissing > 0) || DISTANCE_ADDITION) {
			System.err.println("Calculating quartet distance matrix (for completion of X)");
			calculateDistances(treeAllClusters);
		}

		if (haveMissing > 0 ) {
			this.orderedTaxonBySimilarity = sortByDistance();
			System.err.println("Will attempt to complete bipartitions from X before adding using a distance matrix.");
			int t = 0;
			BufferedWriter completedFile = null;
			if (this.outputCompleted) {
				String fn = GlobalMaps.outputfilename + ".completed_gene_trees";
				System.err.println("Ouputting completed gene trees to " + fn);
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
						completedFile.write(trc.toStringWD()+ " \n");
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
		} else  {
			this.completedGeeneTrees = this.geneTrees;
		}

		/*
		 * Calculate gene tree clusters and bipartitions for X
		 */
		STITreeCluster all = new STITreeCluster();
		all.getBitSet().set(0, n);
		addToClusters(all, GlobalMaps.taxonIdentifier.taxonCount());	
		System.err.println("Building set of clusters (X) from gene trees ");
		
		addTreeBipartitionsToX( this.completedGeeneTrees);
		
		/*
		 * If needed calculate gene tree tripartitions
		 */
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
			for (Entry<Tripartition, Integer> entry : geneTreeTripartitonCount.entrySet()){
				finalTripartitions[i] = entry.getKey();
				finalCounts[i] = entry.getValue();
				i++;
			}
		} else { 
			System.err.println("Using tree-based weight calculation.");
			List<Integer> temp = new ArrayList<Integer>(); 

			for (Tree tr :  this.geneTrees) {
				for (TNode node : tr.postTraverse()) {
					if (node.isLeaf()) {                        
						temp.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
					} else {
						temp.add(-node.getChildCount());
					}
					if (node.isRoot()) {
						temp.add(Integer.MIN_VALUE);
					}
				}
			}
			geneTreesAsInts = temp.toArray(new Integer[]{});
		}


		if (geneTreeTripartitonCount.size() > 0) {
			long s = 0;
			for (Integer c : geneTreeTripartitonCount.values()) {
				s += c;
			}
			System.err.println("Tripartitons in gene trees (count): "
					+ geneTreeTripartitonCount.size());
			System.err.println("Tripartitons in gene trees (sum): " + s);
		}
		System.err.println("Number of Default Clusters: " + clusters.getClusterCount());

		inference.weightCalculator.initializeWeightContainer(
				geneTreeTripartitonCount.size() * 2);
		// System.err.println("sigma n is "+sigmaN);

	}
	
	private void printoutdistmatrix(double [][] distSTMatrix) {
		for (String s :spm.getSTTaxonIdentifier().getAllTaxonNames()) {
	        System.err.print(String.format("%1$8s",s));
	    }
	    System.err.println();
	    for (int i = 0; i < spm.getSpeciesCount(); i++) {
	        for (int j = 0; j< spm.getSpeciesCount(); j++) {
	            System.err.print(String.format("%1$8.3f",distSTMatrix[i][j]));
	        }
	        System.err.println();
	    }
	}

	private void setAlgorithm(int geneTreeTripartitonCountSize, int k){
		if (this.algorithm != -1) {
			return;
		}
		if (k <= 0 || geneTreeTripartitonCountSize <= 0) {
			throw new RuntimeException("gene tree tripartition size or k not set properly");
		}
		if (this.algorithm == -1) {
			this.algorithm = (n <= 32 || (geneTreeTripartitonCountSize < k*6)) ? 2 : 1;
		} else {
			throw new RuntimeException("Algorithm already set");
		}
	}

	private Integer[][] sortByDistance() {
		Integer [][] ret = new Integer[n][n];
		for (int i = 0; i < n; i++) {
			Integer [] indices = new Integer[n];
			for (int j = 0; j < n; j++) {
					indices[j] = j;
			}
			final float[] js = this.similarityMatrix[i];
			Arrays.sort(indices, new Comparator<Integer>() {

				@Override
				public int compare(Integer o1, Integer o2) {
					int comp = Float.compare(js[o1], js[o2]) ;
					return  comp == 0 ? - o1.compareTo(o2) : - comp;
				}
			});
			ret[i] = indices;
		}
		return ret;
	}


	private void updateDistanceForTwoNodes (Integer treeall, BitSet left,
			BitSet right, float[][] matrix) {
		int c = treeall - left.cardinality() - right.cardinality();
		c = c*(c-1)/2;
		for (int l = left.nextSetBit(0); l >= 0; l=left.nextSetBit(l+1)) {
			for (int r = right.nextSetBit(0); r >= 0; r=right.nextSetBit(r+1)) {
				matrix[l][r] += c;
				matrix[r][l] = matrix[l][r];
			}
		}
	}

	private void calculateDistances(List<STITreeCluster> treeAllClusters) {
		Deque<BitSet> stack = new ArrayDeque<BitSet>();
		this.similarityMatrix = new float[n][n];
		int [][] denom = new int [n][n];

		int k = 0;
		for (Tree tree :  this.geneTrees) {
			STITreeCluster treeallCL = treeAllClusters.get(k++);
			
			Integer treeall = treeallCL.getClusterSize();
			
			for (TNode node : tree.postTraverse()) {
				if (node.isLeaf()) {
					BitSet tmp = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
					stack.push(tmp);
				} else if (node.isRoot() && node.getChildCount() == 3){
					BitSet left = stack.pop();
					BitSet middle = stack.pop();
					BitSet right = stack.pop();
					updateDistanceForTwoNodes(treeall, left, right, similarityMatrix);
					updateDistanceForTwoNodes(treeall, left, middle, similarityMatrix);
					updateDistanceForTwoNodes(treeall, middle, right, similarityMatrix);
				} else {
					BitSet left = stack.pop();
					BitSet right = stack.pop();
					BitSet both = new BitSet();
					both.or(left);
					both.or(right);
					BitSet middle = new BitSet();
					middle.or(treeallCL.getBitSet());
					middle.andNot(both); 
					updateDistanceForTwoNodes(treeall, left, right, similarityMatrix);
					updateDistanceForTwoNodes(treeall, left, middle, similarityMatrix);
					updateDistanceForTwoNodes(treeall, middle, right, similarityMatrix);
					stack.push(both);
				}
			}

			BitSet all = treeallCL.getBitSet();
			int c = all.cardinality() - 2;
			for (int l = all.nextSetBit(0); l >= 0; l=all.nextSetBit(l+1)) {
				for (int r = all.nextSetBit(0); r >= 0; r=all.nextSetBit(r+1)) {
					denom[l][r] += c*(c-1)/2;
					denom[r][l] = denom[l][r];
				}
			}
			
		}
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (denom[i][j] == 0)
					similarityMatrix[i][j] = 0;
				else
					similarityMatrix[i][j] = similarityMatrix[i][j] / (denom[i][j]/2);
				similarityMatrix[j][i] = similarityMatrix[i][j];
			}
		}
	}

	public void addExtraBipartitionByDistance() {
		
		float [][] distSTMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
		float[][] denum = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				int stI =  spm.getSpeciesIdForTaxon(i);
				int stJ =  spm.getSpeciesIdForTaxon(j);
				distSTMatrix[stI][stJ] += this.similarityMatrix[i][j]; 
				distSTMatrix[stJ][stI] = distSTMatrix[stI][stJ];
				denum[stI][stJ] ++;
				denum[stJ][stI] ++;
			}
		}
		for (int i = 0; i < spm.getSpeciesCount(); i++) {
			for (int j = 0; j < spm.getSpeciesCount(); j++) {
				distSTMatrix[i][j] = denum[i][j] == 0 ? 0 : 
					distSTMatrix[i][j] / denum[i][j];
			}
			distSTMatrix[i][i] = 1;
			//System.err.println(Arrays.toString(this.distSTMatrix[i]));
		}
		System.err.println("Species tree distances calculated ...");

		ArrayList<Integer> inds = new ArrayList<Integer> (distSTMatrix.length);
		for (int i = 0; i < distSTMatrix.length; i++) {
			inds.add(i);
		}
		for (final float[] fs : distSTMatrix) {
			Collections.sort(inds, new Comparator<Integer>() {

				@Override
				public int compare(Integer i1, Integer i2) {
					return -(Float.compare(fs[i1],fs[i2]));
				}
			});
			
			BitSet stBS = new BitSet(spm.getSpeciesCount());
			//float previous = fs[inds.get(1)];
			//float lastStep = 0;
			for (int sp : inds) {
				stBS.set(sp);
				/*if (previous - fs[sp] < 0) {
					continue;
				}*/
				STITreeCluster g = spm.getGeneClusterForSTCluster(stBS);
				this.addCompletedBipartionToX(g, g.complementaryCluster());
				//lastStep = previous - fs[sp];
				//previous = fs[sp];
			}
			//System.err.println(this.clusters.getClusterCount());
		}

		System.err.println("Number of Clusters after addition by distance: " + clusters.getClusterCount());
	}

	@Override
	public void addExtraBipartitionByExtension(AbstractInference<Tripartition> inference) {
		if (DISTANCE_ADDITION) {
			this.addExtraBipartitionByDistance(); 
		}
	
		System.err.println("Adding to X using resolutions of greedy consensus ...");
		for (Tree tree:  this.completedGeeneTrees) {
			tree.rerootTreeAtNode(tree.getNode(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSpeciesName(0)));
			Trees.removeBinaryNodes((MutableTree) tree);
		}

		/*if (completeTrees.size() < 2) {
			System.err.println("Only "+completeTrees.size() + " complete trees found. Greedy-based completion not applicable.");
			return;
		}*/

		List<Tree> allGreedies = new ArrayList<Tree>();
		
		for (int j = 0; j < 1; j++) {
			allGreedies.addAll(Utils.greedyConsensus(this.completedGeeneTrees,
					GREEDY_ADDITION_THRESHOLDS, true));	
		}
		
		int th = 0;
		for (Tree cons: allGreedies) {
			System.err.println("Threshold " +GREEDY_ADDITION_THRESHOLDS[th++]+":");
			if (th == GREEDY_ADDITION_THRESHOLDS.length) 
				th = 0;
			System.err.println(cons);
			Stack<BitSet> greedyNodeStack = new Stack<BitSet>();
			for (TNode greedyNode :  cons.postTraverse()) {
				if (greedyNode.isLeaf()) {
					BitSet nbs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
					nbs.set(GlobalMaps.taxonIdentifier.taxonId(greedyNode.getName()));
					greedyNodeStack.push(nbs);
					continue;
				} 

				BitSet[] childbs = new BitSet [greedyNode.getChildCount()+1];

				BitSet greedyBS = new BitSet();
				for (int i = 0; i <  greedyNode.getChildCount(); i++) {
					BitSet pop = greedyNodeStack.pop();
					greedyBS.or(pop);
					childbs[i] = pop;
				}
				greedyNodeStack.push(greedyBS);

				if ( greedyNode.getChildCount() > 2 ) { // && greedyNode.getChildCount() < P) {

					BitSet comp = (BitSet) greedyBS.clone();
					comp.flip(0,n);
					childbs[greedyNode.getChildCount()] = comp;

					System.err.print("polytomy of size " + greedyNode.getChildCount());
					
					int k = 0;
					for (int j = 0; j < GREEDY_ADDITION_NOIMPROVEMENT_LIMIT;) {						

						if (!sampleAndResolve(childbs)) {
							j++;
							//System.err.println("+");
						} else {
							k++;
							j--;
						}
					}
					System.err.println("; rounds with additions with at least "+GREEDY_ADDITION_MIN_FREQ+ " support: " + k + "; clusters: "+clusters.getClusterCount());
				}
			}
		}
		System.err.println("Number of Clusters after addition by greedy: " + clusters.getClusterCount());
	}

	private HashMap<BitSet, Integer> returnBitSetCounts(List<Tree> genetrees,
			HashMap<String, Integer> randomSample) {
		
		HashMap<BitSet, Integer> counts = new HashMap<BitSet, Integer>();
		
		for (Tree gt : genetrees) {
			//STITree<Boolean> restrictedTree = new STITree(gt);
			//restrictedTree.constrainByLeaves(randomSample.keySet());	
			List<BitSet> bsList = Utils.getBitsets(randomSample, gt);
			
			for (BitSet bs : bsList) {
				if (counts.containsKey(bs)) {
					counts.put(bs, counts.get(bs) + 1);
					continue;
				}
				BitSet bs2 = (BitSet)bs.clone();
				bs2.flip(0,randomSample.size());
				if (counts.containsKey(bs2)) {
					counts.put(bs2, counts.get(bs2) + 1);
					continue;
				}
				counts.put(bs2, 1);
			}
		}
		return counts;
	}


	
	
	private boolean sampleAndResolve(BitSet[] polytomyBSList) {
		
		boolean addedHighFreq = false;
		// random sample taxa
		HashMap<String, Integer> randomSample = randomSampleAroundPolytomy(polytomyBSList);

		//System.err.print(".");
		// get bipartition counts in the induced trees
		HashMap<BitSet, Integer> counts = returnBitSetCounts(this.completedGeeneTrees, randomSample);
		
		// sort bipartitions
		TreeSet<Entry<BitSet,Integer>> countSorted = new 
				TreeSet<Entry<BitSet,Integer>>(new Utils.BSComparator(true,randomSample.size()));
		countSorted.addAll(counts.entrySet()); 
		
		// build the greedy tree
		MutableTree greedyTree = new STITree<BitSet>();
		TNode [] tmpnodes = new TNode[randomSample.size()];
		for (int i = 0; i < randomSample.size(); i++) {
		  tmpnodes[i] = greedyTree.getRoot().createChild(i+"");
		  BitSet bs = new BitSet(randomSample.size());
		  bs.set(i);
		  ((STINode<BitSet>)tmpnodes[i]).setData(bs);
		}				       
		    
		boolean added = false;
		//System.err.print("^");
		for (Entry<BitSet, Integer> entry : countSorted) {

			BitSet newbs = entry.getKey();

			SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(greedyTree);
			Set<TNode> clusterLeaves = new HashSet<TNode>();
			TNode node;
			for (int i = newbs.nextSetBit(0); i >= 0; i = newbs.nextSetBit(i+1)) {
				node = tmpnodes[i];
				clusterLeaves.add(node);
			}
			TNode lca = lcaFinder.getLCA(clusterLeaves);
			LinkedList<TNode> movedChildren = new LinkedList<TNode>();
			int nodes = clusterLeaves.size();
			for (TNode child : lca.getChildren()) {
				BitSet childCluster = ((STINode<BitSet>)child).getData();

				BitSet temp = (BitSet)childCluster.clone();
				temp.and(newbs);
				if (temp.equals(childCluster)) {
					movedChildren.add(child);
					nodes -= temp.cardinality();
				}

			}

			//boolean isPartOfGreedy = false;
			if ( movedChildren.size() != 0 && nodes == 0) {
				
				STINode<BitSet> newChild = ((STINode<BitSet>)lca).createChild();
				newChild.setData(newbs);
				
				while (!movedChildren.isEmpty()) {
					newChild.adoptChild((TMutableNode)movedChildren.get(0));
					movedChildren.remove(0);
				}

				if (addSubSampledBitSetToX(polytomyBSList, newbs)) {
					if (GREEDY_ADDITION_MIN_FREQ <= (entry.getValue()+.0d)/this.completedGeeneTrees.size()) {
						addedHighFreq = true;
						added  = true;
					}
				}
			}

/*			if (beyondGreedy &&
					GREEDY_ADDITION_MIN_FREQ <= (entry.getValue()+.0d)/this.completedGeeneTrees.size()) {			        		
				if (addSubSampledBitSetToX(polytomyBSList, newbs)) {
					//System.err.print("*");
					added = true;
				}
			} else if (isPartOfGreedy) {
				addSubSampledBitSetToX(polytomyBSList, newbs);
			}*/
		}
		
		
		if (added) {
			for (TNode node : greedyTree.postTraverse()) {
				if (node.getChildCount() > 2) {
					ArrayList<BitSet> children = new ArrayList<BitSet> ();
					for (TNode child : node.getChildren()) {
						children.add(((STINode<BitSet>)child).getData());
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
		}
		return addedHighFreq;
	}


	private HashMap<String,Integer>  randomSampleAroundPolytomy(BitSet[] polyTomy) {
		HashMap<String,Integer>  randomSample = new HashMap<String, Integer>();
		int ind = 0;
		for (BitSet child : polyTomy) {
			int sample = GlobalMaps.random.nextInt(child.cardinality());
			int p = child.nextSetBit(0);
			for (int i = 0; i< sample; i++) {
				p = child.nextSetBit(p+1);
			}
			randomSample.put(GlobalMaps.taxonIdentifier.getTaxonName(p),ind);
			ind++;
		}
		return randomSample;
	}
	
	private boolean sampleAndResolve2(List<Tree> genetrees, BitSet[] childbs) {
		
		boolean added = false;
		HashMap<String,Integer> randomSample = randomSampleAroundPolytomy(childbs);

		// get bipartition counts in the induced trees
		for (Tree gt : genetrees) {
			STITree<Boolean> rgt = new STITree(gt);
			rgt.constrainByLeaves(randomSample.keySet());

			Stack<BitSet> stack = new Stack<BitSet>();
			for (TNode rgtn : rgt.postTraverse()) {

				if (rgtn.isRoot() && rgtn.getChildCount() == 2) {
					continue;
				}
				BitSet bs = new BitSet(randomSample.size());
				if (rgtn.isLeaf()) {
					// Find the index of this leaf.
					int i =  randomSample.get(rgtn.getName());               
					bs.set(i); 
				}
				else {
					for (int i = 0; i <  rgtn.getChildCount(); i++) {
						bs.or(stack.pop());
					}
				}
				stack.push(bs);

				if (bs.cardinality() >= randomSample.size() - 1) {
					continue;
				}
				
				addSubSampledBitSetToX(childbs, bs);    

			}
		}
		
		
		return added;
	}


	private boolean addSubSampledBitSetToX(BitSet[] childbs, BitSet restrictedBitSet) {
		BitSet stNewbs = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
		for (int j = restrictedBitSet.nextSetBit(0); j >= 0; j = restrictedBitSet.nextSetBit(j+1)) {
			stNewbs.or(childbs[j]);
		}
		STITreeCluster g = new STITreeCluster();
		g.setCluster(stNewbs);

		return this.addCompletedBipartionToX(g, g.complementaryCluster());
	}
}
