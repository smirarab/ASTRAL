package phylonet.coalescent;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Stack;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class WQDataCollection extends AbstractDataCollection<Tripartition> {


	List<STITreeCluster> treeAllClusters = new ArrayList<STITreeCluster>();

	Tripartition [] finalTripartitions = null;
	int [] finalCounts = null;	

	Integer [] geneTreesAsInts;

	Integer[][] distMatrix;
	float [][] distSTMatrix;
	private int n;
	private int algorithm;
	private SpeciesMapper spm;


	public WQDataCollection( WQClusterCollection clusters, int alg) {
		this.clusters = clusters;
		this.algorithm = alg;
		this.spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper();
	}


	void traverseTrees(List<Tree> trees, boolean addTripartition,
			Map<Tripartition, Integer> geneTreeTripartitonCount, boolean complete) {

		System.err.println("Building clusters (and more) from gene trees ");
		for (Tree tr : trees) {
			Stack<STITreeCluster> stack = new Stack<STITreeCluster>();
			STITreeCluster gtAll = new STITreeCluster();
			String[] gtLeaves = tr.getLeaves();
			for (int i = 0; i < gtLeaves.length; i++) {
				gtAll.addLeaf(GlobalMaps.taxonIdentifier.taxonId(gtLeaves[i]));
			}
			treeAllClusters.add(gtAll);
			BitSet gtAllBS = gtAll.getBitSet();
			
			int [] neighbor = null;
			
			if (complete) {
				neighbor = new int [GlobalMaps.taxonIdentifier.taxonCount()];
				
				for (int i = gtAllBS.nextClearBit(0); i < n ; i = gtAllBS.nextClearBit(i+1)) {
					
					for (int j = 0; ; j++){
						if ( i > this.distMatrix[i][j] || gtAllBS.get(this.distMatrix[i][j])) {
							neighbor[i] = this.distMatrix[i][j];
							//System.err.println("Mapping " + i + " to "+ this.distMatrix[i][j]);
							break;
						}
						if (j > GlobalMaps.taxonIdentifier.taxonCount()) {
							throw new RuntimeException("Bug: this should not be reached");
						}
					}
				}
			}
			for (TNode node : tr.postTraverse()) {				
				// System.err.println("Node is:" + node);
				if (node.isLeaf()) {
					String nodeName = node.getName(); //GlobalMaps.TaxonNameMap.getSpeciesName(node.getName());

					STITreeCluster cluster = new STITreeCluster();
					Integer taxonID = GlobalMaps.taxonIdentifier.taxonId(nodeName);
					cluster.addLeaf(taxonID);

					STITreeCluster remaining = cluster.complementaryCluster();
					remaining.getBitSet().and(gtAll.getBitSet());

					addBipartitionToX(gtAllBS, cluster,
							remaining, neighbor);

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

					//((STINode)node).setData(new GeneTreeBitset(node.isRoot()? -2: -1));
					stack.add(cluster);


					STITreeCluster remaining = cluster.complementaryCluster();
					remaining.getBitSet().and(gtAll.getBitSet());


					if (addBipartitionToX(gtAllBS, cluster, remaining, neighbor)) {
						//System.err.println("Extra bipartition added: " + spm.getSTClusterForGeneCluster(cluster) +" | "+spm.getSTClusterForGeneCluster(remaining));
					}


					if (addTripartition) {
						if (remaining.getClusterSize() != 0) {
							childbslist.add(remaining);
						}
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

	}

	private boolean addBipartitionToX(BitSet gtAllBS,
			STITreeCluster c1, STITreeCluster c2, int[] neighbor) {

		boolean added = false;

		// TODO: should this be treated differently?
		if (c1.getClusterSize() == 1) {
			added |= addToClusters(c1, c1.getClusterSize());
		}

		STITreeCluster c1c = new STITreeCluster (c1);
		STITreeCluster c2c = new STITreeCluster (c2);
		BitSet b1c = c1c.getBitSet();
		BitSet b2c = c2c.getBitSet();
		
		if (neighbor != null && gtAllBS != null) {
			for (int i = gtAllBS.nextClearBit(0); i < n ; i = gtAllBS.nextClearBit(i+1)) {
				if (b1c.get(neighbor[i])) {
					b1c.set(i);
				} else if (b2c.get(neighbor[i])) {
					b2c.set(i);
				} else {
					throw new RuntimeException("neighbor in neither side");
				}
				//System.err.println("mapping " + i + " to " +  neighbor[i]);
			}
		}		
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


	private boolean addCompletedBipartionToX(STITreeCluster c1, STITreeCluster c2) {
		int size = c1.getClusterSize();
		boolean added = addToClusters(c1, size);  
		size  = c2.getClusterSize();
		added |= addToClusters(c2, size);
		return added;
	}


	public void computeTreePartitions(AbstractInference<Tripartition> inference, boolean addExtra) {

		int k = inference.trees.size();
		n = GlobalMaps.taxonIdentifier.taxonCount();

		int haveMissing = 0;
		for (Tree tree : inference.trees) {
			if (tree.getLeafCount() != n) {
				haveMissing++;
			}
		}
		System.err.println( haveMissing + " trees have missing taxa");

		float[][] geneDists = null; 		
		if ( (haveMissing > k/20) || addExtra) {
			System.err.println("Calculating quartet distance matrix (for completion of X)");
			geneDists = calculateDistances(inference);
		}

		if (haveMissing > k/20 ) {
			this.distMatrix = sortByDistance(geneDists);
			System.err.println("Will attempt to complete bipartitions from X before adding using a distance matrix.");
		}

		if (addExtra) {
			this.distSTMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
			float[][] denum = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					int stI =  spm.getSpeciesIdForTaxon(i);
					int stJ =  spm.getSpeciesIdForTaxon(j);
					this.distSTMatrix[stI][stJ] += geneDists[i][j]; 
					this.distSTMatrix[stJ][stI] = this.distSTMatrix[stI][stJ];
					denum[stI][stJ] ++;
					denum[stJ][stI] ++;
				}
			}
			for (int i = 0; i < spm.getSpeciesCount(); i++) {
				for (int j = 0; j < spm.getSpeciesCount(); j++) {
					this.distSTMatrix[i][j] = denum[i][j] == 0 ? 0 : 
						this.distSTMatrix[i][j] / denum[i][j];
				}
				this.distSTMatrix[i][i] = 1;
				//System.err.println(Arrays.toString(this.distSTMatrix[i]));
			}
			System.err.println("Species tree distances calculated ...");
		}



		/*for (String s :spm.getSTTaxonIdentifier().getAllTaxonNames()) {
            System.err.print(String.format("%1$8s",s));
        }
        System.err.println();
        for (int i = 0; i < spm.getSpeciesCount(); i++) {
            for (int j = 0; j< spm.getSpeciesCount(); j++) {
                System.err.print(String.format("%1$8.3f",distSTMatrix[i][j]));
            }
            System.err.println();
        }*/

		Map<Tripartition, Integer> geneTreeTripartitonCount = new HashMap<Tripartition, Integer>(k * n);

		/*
		 * Calculate gene tree clusters and bipartitions
		 */
		STITreeCluster all = new STITreeCluster();
		all.getBitSet().set(0, n);
		addToClusters(all, GlobalMaps.taxonIdentifier.taxonCount());		
		traverseTrees(inference.trees, true, geneTreeTripartitonCount, this.distMatrix != null);

		this.setAlgorithm(geneTreeTripartitonCount.size(), k);
		
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

			for (Tree tr : inference.trees) {
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


		int s = 0;
		for (Integer c : geneTreeTripartitonCount.values()) {
			s += c;
		}
		System.err.println("Number of gene trees: " + k);
		System.err.println("Tripartitons in gene trees (count): "
				+ geneTreeTripartitonCount.size());
		System.err.println("Tripartitons in gene trees (sum): " + s);
		System.err.println("Number of Clusters: " + clusters.getClusterCount());

		inference.weightCalculator.initializeWeightContainer(
				geneTreeTripartitonCount.size() * 2);
		// System.err.println("sigma n is "+sigmaN);

	}

	private void setAlgorithm(int geneTreeTripartitonCountSize, int k){
		if (k <= 0 || geneTreeTripartitonCountSize <= 0) {
			throw new RuntimeException("gene tree tripartition size or k not set properly");
		}
		if (this.algorithm == -1) {
			this.algorithm = (n <= 32 || (geneTreeTripartitonCountSize < k*6)) ? 2 : 1;
		} else {
			throw new RuntimeException("Algorithm already set");
		}
	}

	private Integer[][] sortByDistance(float[][] geneDists) {
		Integer [][] ret = new Integer[n][n];
		for (int i = 0; i < n; i++) {
			Integer [] indices = new Integer[n];
			for (int j = 0; j < n; j++) {
					indices[j] = j;
			}
			final float[] js = geneDists[i];
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


	private void updateDistanceForTwoNodes(Integer treeall, List<Integer> left,
			List<Integer> right, float[][] matrix) {
		int c = treeall - left.size() - right.size();
		c = c*(c-1)/2;
		for (Integer l : left) {
			for (Integer r : right) {
				matrix[l][r] += c;
				matrix[r][l] = matrix[l][r];
			}
		}
	}

	private float[][] calculateDistances(AbstractInference<Tripartition> inference) {
		Deque<List<Integer>> stack = new ArrayDeque<List<Integer>>();
		float [][] matrix = new float[n][n];
		int [][] denom = new int [n][n];

		for (Tree tree : inference.trees) {
			Integer treeall = tree.getLeafCount();
			for (TNode node : tree.postTraverse()) {
				if (node.isLeaf()) {
					ArrayList<Integer> tmp = new ArrayList<Integer>();
					tmp.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
					stack.push(tmp);
				} else if (node.isRoot() && node.getChildCount() == 3){
					List<Integer> left = stack.pop();
					List<Integer> middle = stack.pop();
					List<Integer> right = stack.pop();
					updateDistanceForTwoNodes(treeall, left, right, matrix);
					updateDistanceForTwoNodes(treeall, left, middle, matrix);
					updateDistanceForTwoNodes(treeall, middle, right, matrix);
					left.addAll(right);
					left.addAll(middle);
					stack.push(left);
				} else {
					List<Integer> left = stack.pop();
					List<Integer> right = stack.pop();
					updateDistanceForTwoNodes(treeall, left, right, matrix);
					left.addAll(right);
					right.clear();
					stack.push(left);
				}
				if (node.isRoot()) {
					List<Integer> all = stack.pop();
					int c = all.size() - 2;
					for (Integer l : all) {
						for (Integer r : all) {
							denom[l][r] += c*(c-1)/2;
							denom[r][l] = denom[l][r];
						}
					}
				}
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (denom[i][j] == 0)
					matrix[i][j] = 0;
				else
					matrix[i][j] = matrix[i][j] / denom[i][j];
				matrix[j][i] = matrix[i][j];
			}
		}

		return matrix;
	}

	public void addExtraBipartitionsByInput(IClusterCollection extraClusters,
			List<Tree> trees, boolean extraTreeRooted) {

		traverseTrees(trees, false, null, this.distMatrix != null);
		int s = extraClusters.getClusterCount();
		/*
		 * for (Integer c: clusters2.keySet()){ s += clusters2.get(c).size(); }
		 */
		System.err
		.println("Number of Clusters after additions from extra trees: "
				+ s);
	}

	public void addExtraBipartitionByDistance() {
		ArrayList<Integer> inds = new ArrayList<Integer> (this.distSTMatrix.length);
		for (int i = 0; i < this.distSTMatrix.length; i++) {
			inds.add(i);
		}
		for (final float[] fs : this.distSTMatrix) {
			Collections.sort(inds, new Comparator<Integer>() {

				@Override
				public int compare(Integer i1, Integer i2) {
					return -(Float.compare(fs[i1],fs[i2]));
				}
			});
			
			BitSet stBS = new BitSet(spm.getSpeciesCount());
			float previous = fs[inds.get(1)];
			float lastStep = 0;
			for (int sp : inds) {
				stBS.set(sp);
				if (previous - fs[sp] < lastStep * 0) {
					continue;
				}
				STITreeCluster g = spm.getGeneClusterForSTCluster(stBS);
				this.addCompletedBipartionToX(g, g.complementaryCluster());
				lastStep = previous - fs[sp];
				previous = fs[sp];
			}
			//System.err.println(this.clusters.getClusterCount());
		}

		System.err.println("Number of Clusters after addition by distance: " + clusters.getClusterCount());
	}

	private void addTripartition(STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster remaining, TNode node,
			Map<Tripartition, Integer> geneTreeTripartitonCount) {

		Tripartition trip = new Tripartition(l_cluster, r_cluster, remaining);
		geneTreeTripartitonCount.put(trip,
				geneTreeTripartitonCount.containsKey(trip) ? 
						geneTreeTripartitonCount.get(trip) + 1 : 1);
	}


	@Override
	public void addExtraBipartitionByExtension() {
		this.addExtraBipartitionByDistance();       
	}



}
