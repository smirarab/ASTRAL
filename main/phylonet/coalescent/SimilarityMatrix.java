package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

/**
 * Implements a Distance method
 * @author smirarab
 *
 */
public class SimilarityMatrix {
	
	private float[][] similarityMatrix;
	private List<TreeSet<Integer>> orderedTaxonBySimilarity;
	private Integer n;
	
	public SimilarityMatrix(int n) {
		this.n = n;
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
	
	private List<TreeSet<Integer>> sortByDistance(float[][] refMatrix) {
		List<TreeSet<Integer>> ret = new ArrayList<TreeSet<Integer>>(n);
		List<Integer> range = Utils.getRange(n);
		for (int i = 0; i < n; i++) {
			final float[] js = refMatrix[i];
			TreeSet<Integer> indices = sortColumn(range, js);
			ret.add(indices);
		}
		return ret;
	}

	private TreeSet<Integer> sortColumn(List<Integer> range, final float[] js) {
		TreeSet<Integer> indices = new TreeSet<Integer>(new Comparator<Integer>() {

			@Override
			public int compare(Integer o1, Integer o2) {
				if (o1 == o2) {
					return 0;
				}
				int comp = Float.compare(js[o1], js[o2]) ;
				return  comp == 0 ? - o1.compareTo(o2) : - comp;
			}
		});
		indices.addAll(range);
		return indices;
	}
	
	private void assureOrderedTaxa () {
		if (this.orderedTaxonBySimilarity == null) {
			this.orderedTaxonBySimilarity = this.sortByDistance(this.similarityMatrix);
		}
	}

	/**
	 * Returns the id of the closest taxon that is either present in presentBS or
	 * has a smaller id than mssingId
	 * @param presentBS
	 * @param missingId
	 * @return
	 */
	int getClosestPresentTaxonId(BitSet presentBS, int missingId) {
		this.assureOrderedTaxa();
		int closestId = -1;
		for (Integer other: this.orderedTaxonBySimilarity.get(missingId)){
			if ( missingId > other // other is already added
					|| presentBS.get(other) // other was in original tree
					) {
				closestId = other;
				break;
			}
		}
		if (closestId == -1) {
			throw new RuntimeException("Bug: this should not be reached");
		}
		return closestId;
	}

	
	/*private void updateQuartetDistanceForPair (Integer treeall, BitSet left,
			BitSet right, float[][] matrix) {
		long c = treeall - left.cardinality() - right.cardinality();
		c = c*(c-1)/2;
		for (int l = left.nextSetBit(0); l >= 0; l=left.nextSetBit(l+1)) {
			for (int r = right.nextSetBit(0); r >= 0; r=right.nextSetBit(r+1)) {
				matrix[l][r] += c;
				matrix[r][l] = matrix[l][r];
			}
		}
	}*/
	
	private void updateQuartetDistanceTri (BitSet left, BitSet right, float[][] matrix, double d) {
		for (int l = left.nextSetBit(0); l >= 0; l=left.nextSetBit(l+1)) {
			for (int r = right.nextSetBit(0); r >= 0; r=right.nextSetBit(r+1)) {
				matrix[l][r] += d;
				matrix[r][l] = matrix[l][r];
			}
		}
	}
	
	/*
	void populateByQuartetDistance(List<STITreeCluster> treeAllClusters, SetWeightCalculator swc) {
		this.similarityMatrix = new float[n][n];
		long [][] denom = new long [n][n];
		
		for (Tripartition t : swc.finalTripartitions) {
			updateQuartetDistanceForPair(t.cluster1.getBitSet(), t.cluster2.getBitSet(), t.cluster3.getBitSet(), similarityMatrix);
			updateQuartetDistanceForPair(t.cluster2.getBitSet(), t.cluster1.getBitSet(), t.cluster3.getBitSet(), similarityMatrix);
			updateQuartetDistanceForPair(t.cluster3.getBitSet(), t.cluster1.getBitSet(), t.cluster2.getBitSet(), similarityMatrix);
			
		}
		int k = 0;
		for (STITreeCluster treeallCL :  treeAllClusters) {
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
				if (i == j) {
					similarityMatrix[i][j] = 1;
				}
				similarityMatrix[j][i] = similarityMatrix[i][j];
			}
			System.err.println(Arrays.toString(similarityMatrix[i]));
		}
	}
	*/
	void populateByQuartetDistance(List<STITreeCluster> treeAllClusters, List<Tree> geneTrees) {

		this.similarityMatrix = new float[n][n];
		long [][] denom = new long [n][n];

		int k = 0;
		for (Tree tree :  geneTrees) {
			
			for (TNode node : tree.postTraverse()) {
				if (node.isLeaf()) {
					BitSet tmp = new BitSet(n);
					tmp.set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
					((STINode)node).setData(tmp);
				} else {
					
					BitSet newbs = new BitSet(n);
					for (TNode cn: node.getChildren()) {
						BitSet c = (BitSet) ((STINode)cn).getData();
						newbs.or(c);
					}
					 
					((STINode)node).setData(newbs);
					
				}
			}
		}
			
		for (Tree tree :  geneTrees) {
			STITreeCluster treeallCL = treeAllClusters.get(k++);
			
			Integer treeall = treeallCL.getClusterSize();
			
			for (TNode node : tree.postTraverse()) {
				if (node.isLeaf()) { 
					continue;
				}
				BitSet cluster = (BitSet) ((STINode)node).getData();
				BitSet others = (BitSet) treeallCL.getBitSet().clone();
				others.andNot(cluster);
				ArrayList<BitSet> children = new ArrayList<BitSet>();
				long totalPairs = 0;
				long totalUnresolvedPairs = 0;
				for (TNode cn: node.getChildren()) {
					BitSet c = (BitSet) ((STINode)cn).getData();
					children.add(c);
					long cc = c.cardinality();
					totalPairs += cc*(cc-1);
					totalUnresolvedPairs += cc * (treeall - cc); 
				}
				if (others.cardinality() != 0) {
					children.add(others);
					long cc = others.cardinality();
					totalPairs += cc*(cc-1);
					totalUnresolvedPairs += cc * (treeall - cc);
				}
				totalPairs /= 2;
				totalUnresolvedPairs /= 2;
				
				
				for (int j = 0; j < children.size(); j++ ) {
					BitSet left = children.get(j);
					long lc = left.cardinality();
					long lcu = lc * (treeall - lc);
					long lcp = lc*(lc-1)/2;
					for (int i = j+1; i < children.size(); i++ ) {
						BitSet right = children.get(i);
						long rc = right.cardinality();
						long rcu = rc * (treeall - lc - rc);
						long rcp = rc*(rc-1)/2;
						double sim = (totalPairs - lcp - rcp) // the number of fully resolved quartets
								//+ (totalUnresolvedPairs - lcu - rcu) / 2.0 // we count partially resolved quartets
								; 
						updateQuartetDistanceTri( left, right, similarityMatrix, sim);
					}
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
				if (i == j) {
					similarityMatrix[i][j] = 1;
				}
				similarityMatrix[j][i] = similarityMatrix[i][j];
			}
			//System.err.println(Arrays.toString(similarityMatrix[i]));
		}
	}
	
	SimilarityMatrix convertToSpeciesDistance(SpeciesMapper spm) {
		float [][] STsimMatrix = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
		float[][] denum = new float[spm.getSpeciesCount()][spm.getSpeciesCount()];
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				int stI =  spm.getSpeciesIdForTaxon(i);
				int stJ =  spm.getSpeciesIdForTaxon(j);
				STsimMatrix[stI][stJ] += this.similarityMatrix[i][j]; 
				STsimMatrix[stJ][stI] = STsimMatrix[stI][stJ];
				denum[stI][stJ] ++;
				denum[stJ][stI] ++;
			}
		}
		for (int i = 0; i < spm.getSpeciesCount(); i++) {
			for (int j = 0; j < spm.getSpeciesCount(); j++) {
				STsimMatrix[i][j] = denum[i][j] == 0 ? 0 : 
					STsimMatrix[i][j] / denum[i][j];
			}
			STsimMatrix[i][i] = 1;
			//System.err.println(Arrays.toString(this.distSTMatrix[i]));
		}
		System.err.println("Species tree distances calculated ...");
		
		SimilarityMatrix ret = new SimilarityMatrix(STsimMatrix.length);
		ret.similarityMatrix = STsimMatrix;
		
		return ret;
	}
	
	SimilarityMatrix getInducedMatrix(HashMap<String, Integer> randomSample) {
		
		int sampleSize = randomSample.size();
		float[][] sampleSimMatrix = new float [sampleSize][sampleSize];
		
		for (Entry<String, Integer> row : randomSample.entrySet()) {
			int rowI = GlobalMaps.taxonIdentifier.taxonId(row.getKey());
			int i = row.getValue();
			for (Entry<String, Integer> col : randomSample.entrySet()) {
				int colJ = GlobalMaps.taxonIdentifier.taxonId(col.getKey());
				sampleSimMatrix[i][col.getValue()] = this.similarityMatrix[rowI][colJ];
			}
		}
		SimilarityMatrix ret = new SimilarityMatrix(sampleSize);
		ret.similarityMatrix = sampleSimMatrix;
		return ret;
	}

	//TODO: generate iterable, not list
	Iterable<BitSet> getQuadraticBitsets() {
		List<BitSet> newBitSets = new ArrayList<BitSet>();
		ArrayList<Integer> inds = new ArrayList<Integer> (n);
		for (int i = 0; i < n; i++) {
			inds.add(i);
		}
		for (final float[] fs : this.similarityMatrix) {
			Collections.sort(inds, new Comparator<Integer>() {

				@Override
				public int compare(Integer i1, Integer i2) {
					if (i1 == i2) {
						return 0;
					}
					int vc = Float.compare(fs[i1],fs[i2]);
					if (vc != 0) {
						return - vc;
					}
					return i1 > i2 ? 1 : -1;
				}
			});
			BitSet stBS = new BitSet(n);
			//float previous = fs[inds.get(1)];
			//float lastStep = 0;
			for (int sp : inds) {
				stBS.set(sp);
				/*if (previous - fs[sp] < 0) {
						continue;
					}*/
				newBitSets.add((BitSet) stBS.clone());
				//lastStep = previous - fs[sp];
				//previous = fs[sp];
			}
			//System.err.println(this.clusters.getClusterCount());
		}
		return newBitSets;
	}
	
	
	List<BitSet> resolveByUPGMA(List<BitSet> bsList, boolean original, boolean randomize) {
		
		List<BitSet> internalBSList;
		if (original) {
			internalBSList = new ArrayList<BitSet>(bsList);
		} else {
			internalBSList = new ArrayList<BitSet>();
		}
		
		int size = bsList .size();
		List<TreeSet<Integer>> indsBySim = new ArrayList<TreeSet<Integer>>(size);
		List<float[]> sims = new ArrayList<float[]>(size);
		List<Integer> range = Utils.getRange(size);
		List<Integer> weights = new ArrayList<Integer>(size);
		
		for (int i = 0; i < size; i++) {
			if (!original) {
				BitSet internalBS = new BitSet(size);
				internalBS.set(i);
				internalBSList.add(internalBS);
			}
			
			final float[] is = new float[size];// this.similarityMatrix[i].clone();
			BitSet bsI = bsList.get(i);
			weights.add(bsI.cardinality());
			sims.add(is);
			
			for (int j = 0; j < size; j++) {
				
				BitSet bsJ = bsList.get(j);
				int c = 0;
				if (i == j) {
					is[j] = 1;
					continue;
				}
				for (int k = bsI.nextSetBit(0); k >= 0; k = bsI.nextSetBit(k + 1)) {
					for (int l = bsJ.nextSetBit(0); l >= 0; l = bsJ.nextSetBit(l + 1)) {
						is[j] += this.similarityMatrix[k][l];
						c++;
					}
				}
				if (c == 0) {
					throw new RuntimeException("Error: "+bsI + " "+bsJ);
				}
				is[j] /= c;
			}
		
			range.remove(i);
			TreeSet<Integer> sortColumn = this.sortColumn(range, is);
			range.add(i,i);
			indsBySim.add(sortColumn);
		}
		
		return upgmaLoop(weights, internalBSList, indsBySim, sims, size, randomize);
	}
	
	List<BitSet> UPGMA() {
		
		List<BitSet> bsList = new ArrayList<BitSet>(n);
		List<TreeSet<Integer>> indsBySim = new ArrayList<TreeSet<Integer>>(n);
		List<float[]> sims = new ArrayList<float[]>(n);
		List<Integer> range = Utils.getRange(n);
		List<Integer> weights = Utils.getOnes(n);
		
		for (int i = 0; i< n; i++) {
			BitSet bs = new BitSet();
			bs.set(i);
			bsList.add(bs);
			final float[] is = this.similarityMatrix[i].clone();
			sims.add(is);
			range.remove(i);
			TreeSet<Integer> sortColumn = this.sortColumn(range, is);
			range.add(i,i);
			indsBySim.add(sortColumn);
		}
		
		return upgmaLoop(weights, bsList, indsBySim, sims, n, false);
	}

	private List<BitSet> upgmaLoop(List<Integer> weights, List<BitSet> bsList,
			List<TreeSet<Integer>> indsBySim, List<float[]> sims, int left, boolean randomize) {
		List<BitSet> ret = new ArrayList<BitSet>();
		while ( left > 2) {
			int closestI = -1;
			int closestJ = -1;
			float bestHit = -1;
			for (int i = 0; i < indsBySim.size(); i++) {
				if (indsBySim.get(i) == null)
					continue;
				int j = indsBySim.get(i).first();
				if (sims.get(i)[j] > bestHit || (randomize & sims.get(i)[i] == bestHit & GlobalMaps.random.nextBoolean())) {
					bestHit = sims.get(i)[j];
					closestI = i;
					closestJ = j;
				}
			}
			BitSet bs = (BitSet) bsList.get(closestI).clone();
			bs.or(bsList.get(closestJ));
			bsList.set(closestJ,null);
			bsList.set(closestI,bs);
			
			float[] jDist = sims.get(closestJ);
			float[] iDist = sims.get(closestI).clone();
			for (int k = 0; k < sims.size(); k++) {
				if (k == closestJ || sims.get(k) == null) {
					continue;
				}
				
				if ( k != closestI) {
					float newSimToI = (iDist[k] * weights.get(closestI) + jDist[k] * weights.get(closestJ))/( weights.get(closestI)+ weights.get(closestJ));
					
					indsBySim.get(k).remove(closestI);
					sims.get(k)[closestI] = newSimToI;
					indsBySim.get(k).add(closestI);
					
					indsBySim.get(closestI).remove(k);
					sims.get(closestI)[k] = newSimToI;
					indsBySim.get(closestI).add(k);
				}
			
				indsBySim.get(k).remove(closestJ);
				sims.get(k)[closestJ] = -1;
				//indsBySim.get(k).add(closestJ);
			}
			
			sims.set(closestJ,null);
			indsBySim.set(closestJ,null);
			weights.set(closestI, weights.get(closestI)+weights.get(closestJ));
			weights.set(closestJ,null);
			ret.add(bs);
			left--;
		}
		return ret;
	}
}
