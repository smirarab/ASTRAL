package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

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
	
	public SimilarityMatrix(float[][] from) {
		this.n = from.length;
		this.similarityMatrix = new float[from.length][from[0].length];
		for (int i = 0; i < from.length; i++) {
			float[] l = from[i];
			for (int j = 0; j < l.length; j++) {
				this.similarityMatrix[i][j] = l[j];
			}
		}
	}
	
	public int getSize() {
		return n;
	}
	
	public float get(int i, int j) {
		return this.similarityMatrix[i][j];
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

	
	
	private void updateQuartetDistanceTri(BitSet left,
			BitSet right, float[][] matrix, float d) {
		for (int l = left.nextSetBit(0); l >= 0; l=left.nextSetBit(l+1)) {
			for (int r = right.nextSetBit(0); r >= 0; r=right.nextSetBit(r+1)) {
				matrix[l][r] += d;
				matrix[r][l] = matrix[l][r];
			}
		}
	}
	
	void populateByQuartetDistance(List<STITreeCluster> treeAllClusters, List<Tree> geneTrees) {
		this.similarityMatrix = new float[n][n];
		long [][] denom = new long [n][n];
		//fillZero2D(this.similarityMatrix);
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				denom[i][j] = 0L;
			}
		}
		if(CommandLine.timerOn) {
			System.err.println("TIME TOOK FROM LAST NOTICE SimilarityMatrix 145-148: " + (double)(System.nanoTime()-CommandLine.timer)/1000000000);
			CommandLine.timer = System.nanoTime();
			System.out.println(n);
		}
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

		ArrayList<Future<float[][][]>> futures = new ArrayList<Future<float[][][]>>();
		for (int i = 0; i < geneTrees.size(); i+=geneTrees.size()/CommandLine.numThreads) {
			int start = i;
			int end = Math.min(start + geneTrees.size()/CommandLine.numThreads, geneTrees.size());
			futures.add((Future) CommandLine.eService.submit(new populateByQuartetDistanceLoop(start, end, treeAllClusters, geneTrees)));
		}
		for (Future future: futures) {
			float[][][] res = null;
			try {
				res = (float[][][])future.get();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < n; j++) {
					similarityMatrix[i][j] += res[0][i][j];
					denom[i][j] += res[1][i][j];
				}
			}
		}
		if(CommandLine.timerOn) {
			System.err.println("TIME TOOK FROM LAST NOTICE SimilarityMatrix 161-164: " + (double)(System.nanoTime()-CommandLine.timer)/1000000000);
			CommandLine.timer = System.nanoTime();
		}
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (denom[i][j] == 0)
					similarityMatrix[i][j] = 0F;
				else
					similarityMatrix[i][j] = similarityMatrix[i][j] / (denom[i][j]/2);
				if (i == j) {
					similarityMatrix[i][j] = 1F;
				}
				similarityMatrix[j][i] = similarityMatrix[i][j];
			}
		}
	}
	public class populateByQuartetDistanceLoop implements Callable<float[][][]>{
		List<STITreeCluster> treeAllClusters;
		List<Tree> geneTrees;
		int start;
		int end;
		public populateByQuartetDistanceLoop(int start, int end, List<STITreeCluster> treeAllClusters, List<Tree> geneTrees) {
			this.geneTrees = geneTrees;
			this.start = start;
			this.end = end;
			this.treeAllClusters = treeAllClusters;
		}
		public float[][][] call() {
			float[][] array = new float[n][n];
			float[][] denom = new float[n][n];
			for(int w = start; w < end; w++) {
				STITreeCluster treeallCL = treeAllClusters.get(w);
				Tree tree = geneTrees.get(w);
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
								float sim = (totalPairs - lcp - rcp) // the number of fully resolved quartets
										//+ (totalUnresolvedPairs - lcu - rcu) / 3.0 // we count partially resolved quartets
										; 
								updateQuartetDistanceTri( left, right, array, sim);
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
			return new float[][][]{array, denom};
		}
	}
	
	
	SimilarityMatrix getInducedMatrix(HashMap<String, Integer> randomSample, TaxonIdentifier id) {
		
		int sampleSize = randomSample.size();
		float[][] sampleSimMatrix = new float [sampleSize][sampleSize];
		
		for (Entry<String, Integer> row : randomSample.entrySet()) {
			int rowI = id.taxonId(row.getKey());
			int i = row.getValue();
			for (Entry<String, Integer> col : randomSample.entrySet()) {
				int colJ = id.taxonId(col.getKey());
				sampleSimMatrix[i][col.getValue()] = this.similarityMatrix[rowI][colJ];
			}
		}
		SimilarityMatrix ret = new SimilarityMatrix(sampleSize);
		ret.similarityMatrix = sampleSimMatrix;
		return ret;
	}

	SimilarityMatrix getInducedMatrix(List<Integer> sampleOrigIDs) {
		
		int sampleSize = sampleOrigIDs.size();
		SimilarityMatrix ret = new SimilarityMatrix(sampleSize);
		ret.similarityMatrix = new float [sampleSize][sampleSize];
		
		int i = 0;
		for (Integer rowI : sampleOrigIDs) {
			int j = 0;
			for (Integer colJ : sampleOrigIDs) {
				ret.similarityMatrix[i][j] = this.similarityMatrix[rowI][colJ];
				j++;
			}
			i++;
		}
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
			//Float previous = fs[inds.get(1)];
			//Float lastStep = 0;
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
	
	
	List<BitSet> resolveByUPGMA(List<BitSet> bsList, boolean original) {
		
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
					is[j] = 1F;
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
		
		return upgmaLoop(weights, internalBSList, indsBySim, sims, size,false);
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
			List<TreeSet<Integer>> indsBySim, List<float[]> sims, int left,boolean randomize) {
		List<BitSet> ret = new ArrayList<BitSet>();
		while ( left > 2) {
			int closestI = -1;
			int closestJ = -1;
			Float bestHit = -1F;
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
					Float newSimToI = (iDist[k] * weights.get(closestI) + jDist[k] * weights.get(closestJ))/( weights.get(closestI)+ weights.get(closestJ));
					
					indsBySim.get(k).remove(closestI);
					sims.get(k)[closestI] = newSimToI;
					indsBySim.get(k).add(closestI);
					
					indsBySim.get(closestI).remove(k);
					sims.get(closestI)[k] = newSimToI;
					indsBySim.get(closestI).add(k);
				}
			
				indsBySim.get(k).remove(closestJ);
				sims.get(k)[closestJ] = -1F;
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
	public void fillZero2D(Float[][] array) {
		for(int i = 0; i < array.length; i++) {
			for(int j = 0; j < array[0].length; j++) {
				array[i][j] = 0F;
			}
		}
	}
}
