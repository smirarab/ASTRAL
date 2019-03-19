package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.HashOnlyTreeCluster;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class Polytree {
	
	static long time = 0;
	
	static long F(int[] x, int[] y, int[] z){
		long a = x[0], b = x[1], c = x[2], d = y[0], e = y[1], f = y[2], g = z[0], h = z[1], i = z[2];
		return a * ( (a + e + i - 3l)  * e * i + (a + f + h - 3l)  * f * h )
			 + b * ((b + d + i - 3l)  * d * i + (b + f + g - 3l)  * f * g )
			 + c * ((c + d + h - 3l)  * d * h + (c + e + g - 3l)  * e * g );
	}
	
	static long U(int[] x, int[] y, int[] z){
		long a = x[0], d = x[1], b = y[0], e = y[1], c = z[0], f = z[1], s = a + b + c - 3;
		return ((s - c + f) * b * f + (s - b + e) * e * c) * a + ((s - a + d) * d - 2 * s * a) * b * c;
	}
	
	/**
	 * A node in the input gene trees. Used temporarily to build the polytree. 
	 * @author smirarab
	 *
	 */
	final class PTNode{
		PTNode parent;
		ArrayList<PTNode> children;
		PTCluster cluster;
		PTPartition partition;
		boolean isUsed = false;
		
		/**
		 * To be used for leaves
		 * @param n
		 */

		PTNode(TNode n){
			HashOnlyTreeCluster c = new HashOnlyTreeCluster(GlobalMaps.taxonIdentifier.taxonId(n.getName()));
			cluster = Polytree.this.clusters.get(c);
			children = new ArrayList<PTNode>();
		}
		PTNode(ArrayList<PTNode> ch, HashOnlyTreeCluster s){
			children = ch;
			HashOnlyTreeCluster c =  new HashOnlyTreeCluster();
			ArrayList<STITreeCluster> cs = new ArrayList<STITreeCluster>();
			for (PTNode child: children){
				child.parent = this;
				c = c.disjointClusterMerge(child.cluster.clusterRef);
				cs.add(child.cluster.clusterRef);
			}
			cluster = findCluster(c, this);
			if (c.equals(s) == false){ // If this is not the root
				HashOnlyTreeCluster xc = s.subclusterComplement(c);
				cs.add(findCluster(xc, null).clusterRef);
			}
			if (cs.size() >= 3){
				AbstractPartition p = AbstractPartition.createPartition(cs);
				partition = Polytree.this.partitions.get(p);
				if (partition != null){
					partition.cardinality++;
				}
				else partition = new PTPartition(p, this);
			}
		}
		PTCluster findCluster(HashOnlyTreeCluster c, PTNode n){
			PTCluster cluster = Polytree.this.clusters.get(c);
			if (cluster != null) {
				if (cluster.firstNode == null) cluster.firstNode = n;
				return cluster;
			}
			else return new PTCluster(c, n);
		}
		boolean isFirstResolutionOfCluster(){
			if (partition == null) return false; // leaf node
			return partition.firstNode == this;
		}
		boolean isFirstClusterAppearence(){
			return cluster.firstNode == this;
		}
		
		/**
		 * This function sets the isUsed flag based on whether
		 *   1. the cluster is being seen for the first time
		 *      or 
		 *   2. the partition is being seen for the first time
		 */
		void setClusterFlag(){
			for (PTNode child: children){
				child.setClusterFlag();
			}
			if (isFirstResolutionOfCluster()){
				isUsed = true;
				for (PTNode child: children){
					child.setClusterFlagsByClusterPrecedence();
				}
			}
		}
		void setClusterFlagsByClusterPrecedence(){
			if (this.isUsed) return;
			if (isFirstClusterAppearence()){
				this.isUsed = true;
				for (PTNode child: children){
					child.setClusterFlagsByClusterPrecedence();
				}
			}
			else {
				cluster.intersectionAlreadyComputed = true;
				PTNode n = cluster.firstNode;
				if (n == null || n.isUsed) return;
				n.isUsed = true;
				for (PTNode child: n.children){
					child.setClusterFlagsByClusterPrecedence();
				}
			}
		}
		void buildInstructionQueue(){
			/*
			 * -1 - new tree
			 * 
			 * If the cluster is used, the five LS bits are flags, with meanings: 
			 * 1 - true = compute based on stack; 
			 *     false = get from list stack: value >> 5 - number of children
			 * 2 - true = store on stack
			 * 4 - true = store on list
			 * 8 - true = compute partition
			 * 16 - true = partition seen multiple times 
			 * and the remaining bits are simply the number of children. 
			 * When 16 is set, the following value is simply the cardinality fo the partition.
			 * 
			 * If the cluster is not used, bit 1 is not set and 
			 *   the remaining bits give the position in the list where the results for this can be found. 
			 * list:
			 * value >> 1 - the position on the list to fetch from
			 */
			for (PTNode child: children){
				child.buildInstructionQueue();
			}
			if (isUsed){
				int v = (children.size() << 5) | 1;
				if (addToStack()) v = v | 2;
				if (addToList()) {
					v = v | 4;
					cluster.listPos = Polytree.this.listSize++;
				}
				if (isFirstResolutionOfCluster()){
					v = v | 8;
					if (partition.cardinality > 1){
						v = v | 16;
						Polytree.this.queueBuilder.add(v);
						Polytree.this.queueBuilder.add(partition.cardinality); // The cardinality is saved on the queue
					}
					else Polytree.this.queueBuilder.add(v);
				}
				else Polytree.this.queueBuilder.add(v);
			}
			else {
				if (parent != null && parent.isUsed) Polytree.this.queueBuilder.add(cluster.listPos << 1);
			}
		}
		
		/**
		 * The intersection for this cluster has to be saved on stack
		 * @return
		 */
		boolean addToStack(){
			return isUsed && parent != null && parent.isUsed;
		}
		
		/**
		 * The intersection for this cluster has to be saved on the end of the list 
		 * @return
		 */
		boolean addToList(){
			return isUsed && isFirstClusterAppearence() && cluster.intersectionAlreadyComputed;
		}
	}
	
	final class PTCluster{
		HashOnlyTreeCluster clusterRef;
		PTNode firstNode; // The first gene tree node that matched this cluster
		boolean intersectionAlreadyComputed = false; 
		int listPos = -1;
		
		PTCluster(HashOnlyTreeCluster c, PTNode n){
			clusterRef = c;
			firstNode = n;
			Polytree.this.clusters.put(c, this);
		}
		PTCluster(HashOnlyTreeCluster c){
			clusterRef = c;
			firstNode = null;
			intersectionAlreadyComputed = true;
			listPos = Polytree.this.listSize++;
			Polytree.this.clusters.put(c, this);
		}
	}
	
	final class PTPartition{
		PTNode firstNode;
		int cardinality = 1; // cardinality of the partition in the gene tree set
		
		PTPartition(AbstractPartition p, PTNode n){
			firstNode = n;
			Polytree.this.partitions.put(p, this);
		}
	}
	
	public static final class PTNative{
		private static final int batchSize = 32;
		

		static native void cppInit(int n, int listSize, int[] q, long[][] c);
		static native void cppBatchCompute(long[] result, long[][] a, long[][] b, long[][] c);

	}
	
	WQDataCollection dataCollection;
	HashMap<HashOnlyTreeCluster, PTCluster> clusters = new HashMap<HashOnlyTreeCluster, PTCluster>();	
	HashMap<AbstractPartition, PTPartition> partitions = new HashMap<AbstractPartition, PTPartition>();
	ArrayList<PTNode> nodeRoots = new ArrayList<PTNode>();

	ArrayList<Integer> queueBuilder = new ArrayList<Integer>();
	int[] queue;
	int listSize = 0;
	long maxScore = 0;
	private boolean useNativeMethod;
	
	public Polytree(List<Tree> trees, WQDataCollection dataCollection){
		
		this.dataCollection = dataCollection;
		long t = System.currentTimeMillis();
		
		// Create singleton clusters and add to map
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			HashOnlyTreeCluster c = new HashOnlyTreeCluster(i);
			new PTCluster(c);
		}
		
		// Represent gene trees as PTNodes
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			nodeRoots.add(buildTree(tr.getRoot(), new HashOnlyTreeCluster(tit.next())));
		}
		
		// Set the isUsedFlag on gene tree nodes
		for (PTNode n: nodeRoots){
			n.setClusterFlag();
		}
		
		// For each node of each gene tree, decide how it
		//  should be treated when calculating weights.
		// Options (non-exclusive) are:
		//   - compute the intersection for the cluster or retrieve it from the list
		//   - compute the weight for the partition if it's the first resolution
		for (PTNode n: nodeRoots){
			queueBuilder.add(-1);
			n.buildInstructionQueue();
		}
		

		queue = mapToInt(queueBuilder);
		clusters = null;
		partitions = null;
		queueBuilder = null;

		STITreeCluster c = (new STITreeCluster(GlobalMaps.taxonIdentifier)).complementaryCluster();
		maxScore = computeUpperbound(c.getBitSet());
		System.err.println("Polytree max score: " + maxScore / 4);
		System.err.println("Polytree building time: " + (System.currentTimeMillis() - t) / 1000.0D + " seconds.");
		
		try {
			System.loadLibrary("Astral");
			System.err.println("Using native AVX batch computing method.");
			useNativeMethod = true;
		}
		catch (Throwable e) {
			useNativeMethod = false;
			//e.printStackTrace(); 
			System.err.println("Fail to load native library "+System.mapLibraryName("Astral")+"; use Java default computing method.");
		}
		
		if (useNativeMethod) {
			int m = trees.size();
			long b[][] = new long[m][];
			Iterator<STITreeCluster> ti = dataCollection.treeAllClusters.iterator();
			for (int i = 0; i < m; i++) {
				b[i] = ti.next().getBitSet().getArray();
			}
			PTNative.cppInit(GlobalMaps.taxonIdentifier.taxonCount(), listSize, queue, b);
		}
	}
	
	private int[] mapToInt(List<Integer> list) {
		int[] ret = new int [list.size()]; 
		for (int i = 0; i < ret.length; i++)
			ret[i] = list.get(i);
		return ret;
	}
	

	private PTNode buildTree(TNode node, HashOnlyTreeCluster s){
		if (node.isLeaf()) return new PTNode(node);
		else {
			ArrayList<PTNode> cs = new ArrayList<PTNode>();
			for (TNode ch: node.getChildren()){
				cs.add(buildTree(ch, s));
			}
			return new PTNode(cs, s);
		}
	}

	public Long[] WQWeightByTraversal(Tripartition[] trips){
		long t = System.nanoTime();
		Long[] ret = new Long[trips.length];
		if (!useNativeMethod) {
			int i = 0;
			for (Tripartition trip: trips) {
				ret[i++] = this.WQWeightByTraversal(trip);
			}
		} else {
			for (int i = 0; i < trips.length; i += PTNative.batchSize) {
				int size = (trips.length - i > PTNative.batchSize) ? PTNative.batchSize : trips.length - i;
				long[] result = new long[size];
				long[][] a = new long[size][];
				long[][] b = new long[size][];
				long[][] c = new long[size][];
				for (int j = 0; j < size; j++) {
					a[j] = trips[(i + j)].cluster1.getBitSet().getArray();
					b[j] = trips[(i + j)].cluster2.getBitSet().getArray();
					c[j] = trips[(i + j)].cluster3.getBitSet().getArray();
				}
				PTNative.cppBatchCompute(result, a, b, c);
				for (int j = 0; j < size; j++) {
					ret[i+j] = result[j];
				}
			}
		}
		Polytree.time += System.nanoTime() - t;
		//System.err.println((System.nanoTime() - t));
		return ret;
		
	}
	public Long WQWeightByTraversal(Tripartition trip){

		if (trip == null)
		{
			System.err.println("why here?");
		}
		if (trip.cluster1 == trip.cluster2) return this.computeUpperbound(trip.cluster1.getBitSet());
		//long t = System.nanoTime();
		BitSet[] b = new BitSet[]{trip.cluster1.getBitSet(), trip.cluster2.getBitSet(), trip.cluster3.getBitSet()};
		return this.WQWeightByTraversal(b);

	}
	
	public Long WQWeightByTraversal(BitSet[] b){
		int[][] stack, list;
		stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		list = new int[listSize][3];
		int[] treeTotal = new int[3];
		long[] sx = new long[3], sxy = new long[3];
		long weight = 0;
		int stackEnd = 0, listEnd = GlobalMaps.taxonIdentifier.taxonCount();
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (int i = 0, i_end = GlobalMaps.taxonIdentifier.taxonCount(); i < i_end; i++){
			list[i][0] = b[0].get(i) ? 1 : 0;
			list[i][1] = b[1].get(i) ? 1 : 0;
			list[i][2] = b[2].get(i) ? 1 : 0;
		}
		for (int i = 0, i_end = queue.length; i < i_end; i++){
			int cmd = queue[i];
			if (cmd == -1) {
				BitSet all = tit.next().getBitSet();
				treeTotal[0] = b[0].intersectionSize(all);
				treeTotal[1] = b[1].intersectionSize(all);
				treeTotal[2] = b[2].intersectionSize(all);
				continue;
			}
			if ((cmd & 1) != 0){
				int numChildren = cmd >> 5;
				long tempWeight = 0;
				int[] p = stack[stackEnd], q;
				p[0] = treeTotal[0];
				p[1] = treeTotal[1];
				p[2] = treeTotal[2];
				for (int j = stackEnd - numChildren; j < stackEnd; j++){
					q = stack[j];
					p[0] -= q[0];
					p[1] -= q[1];
					p[2] -= q[2];
				}
				if ((cmd & 8) != 0){
					if (numChildren == 2){
						tempWeight = F(stack[stackEnd - 2], stack[stackEnd - 1], stack[stackEnd]);
					}
					else{
						sx[0] = 0;
						sx[1] = 0;
						sx[2] = 0;
						sxy[0] = 0;
						sxy[1] = 0;
						sxy[2] = 0;
						for (int j = stackEnd - numChildren; j <= stackEnd; j++){
							q = stack[j];
							sx[0] += q[0];
							sx[1] += q[1];
							sx[2] += q[2];
							sxy[0] += q[1] * q[2];
							sxy[1] += q[2] * q[0];
							sxy[2] += q[0] * q[1];
						}
						for (int j = stackEnd - numChildren; j <= stackEnd; j++){
							q = stack[j];
							tempWeight += ((sx[1] - q[1]) * (sx[2] - q[2]) - sxy[0] + q[1] * q[2]) * q[0] * (q[0] - 1L)
								+ ((sx[2] - q[2]) * (sx[0] - q[0]) - sxy[1] + q[2] * q[0]) * q[1] * (q[1] - 1L)
								+ ((sx[0] - q[0]) * (sx[1] - q[1]) - sxy[2] + q[0] * q[1]) * q[2] * (q[2] - 1L);
						}
					}
					if ((cmd & 16) != 0) weight += tempWeight * queue[++i];
					else weight += tempWeight;
				}
				stackEnd -= numChildren;
				if ((cmd & 2) != 0){
					q = stack[stackEnd++];
					q[0] = treeTotal[0] - p[0];
					q[1] = treeTotal[1] - p[1];
					q[2] = treeTotal[2] - p[2];
				}
				if ((cmd & 4) != 0){
					q = list[listEnd++];
					q[0] = treeTotal[0] - p[0];
					q[1] = treeTotal[1] - p[1];
					q[2] = treeTotal[2] - p[2];
				}
			}
			else {
				int[] p = list[cmd >> 1], q = stack[stackEnd++];
				q[0] = p[0];
				q[1] = p[1];
				q[2] = p[2];
			}
		}
		//System.err.println(weight);
		return weight;
	}
	
	public Long computeUpperbound(BitSet b){
		int[][] stack, list;
		stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		list = new int[listSize][3];
		int[] treeTotal = new int[3];
		long[] sx = new long[3], sxy = new long[3];
		long weight = 0;
		int stackEnd = 0, listEnd = GlobalMaps.taxonIdentifier.taxonCount();
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (int i = 0, i_end = GlobalMaps.taxonIdentifier.taxonCount(); i < i_end; i++){
			list[i][0] = b.get(i) ? 1 : 0;
			list[i][1] = 1;
		}
		for (int i = 0, i_end = queue.length; i < i_end; i++){
			int cmd = queue[i];
			if (cmd == -1) {
				BitSet all = tit.next().getBitSet();
				treeTotal[0] = b.intersectionSize(all);
				treeTotal[1] = all.cardinality();
				continue;
			}
			if ((cmd & 1) != 0){
				int numChildren = cmd >> 5;
				long tempWeight = 0;
				int[] p = stack[stackEnd], q;
				p[0] = treeTotal[0];
				p[1] = treeTotal[1];
				for (int j = stackEnd - numChildren; j < stackEnd; j++){
					q = stack[j];
					p[0] -= q[0];
					p[1] -= q[1];
				}
				if ((cmd & 8) != 0){
					if (numChildren == 2){
						tempWeight = U(stack[stackEnd - 2], stack[stackEnd - 1], stack[stackEnd]);
					}
					else{
						sx[0] = 0;
						sx[1] = 0;
						sxy[0] = 0;
						sxy[1] = 0;
						for (int j = stackEnd - numChildren; j <= stackEnd; j++){
							q = stack[j];
							sx[0] += q[0];
							sx[1] += q[1];
							sxy[0] += q[0] * q[0];
							sxy[1] += q[0] * q[1];
						}
						for (int j = stackEnd - numChildren; j <= stackEnd; j++){
							q = stack[j];
							tempWeight += ((sx[0] - q[0]) * (sx[1] - q[1]) - sxy[1] + q[0] * q[1]) * q[0] * (q[0] - 1L)
								+ ((sx[0] - q[0]) * (sx[0] - q[0]) - sxy[0] + q[0] * q[0]) * (q[1] * (q[1] - 1L) / 2 - q[0] * (q[0] - 1L));
						}
					}
					if ((cmd & 16) != 0) weight += tempWeight * queue[++i];
					else weight += tempWeight;
				}
				stackEnd -= numChildren;
				if ((cmd & 2) != 0){
					q = stack[stackEnd++];
					q[0] = treeTotal[0] - p[0];
					q[1] = treeTotal[1] - p[1];
				}
				if ((cmd & 4) != 0){
					q = list[listEnd++];
					q[0] = treeTotal[0] - p[0];
					q[1] = treeTotal[1] - p[1];
				}
			}
			else {
				int[] p = list[cmd >> 1], q = stack[stackEnd++];
				q[0] = p[0];
				q[1] = p[1];
			}
		}
		return weight;
	}
}
