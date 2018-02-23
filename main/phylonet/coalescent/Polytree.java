package phylonet.coalescent;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.coalescent.WQWeightCalculator.CondensedTraversalWeightCalculator;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
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
	
	final class PTNode{
		PTNode parent;
		ArrayList<PTNode> children;
		PTCluster cluster;
		PTPartition partition;
		boolean called = false;
		
		PTNode(TNode n){
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(GlobalMaps.taxonIdentifier.taxonId(n.getName()));
			cluster = Polytree.this.clusters.get(c);
			children = new ArrayList<PTNode>();
		}
		PTNode(ArrayList<PTNode> ch, STITreeCluster s){
			children = ch;
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			ArrayList<STITreeCluster> cs = new ArrayList<STITreeCluster>();
			for (PTNode child: children){
				child.parent = this;
				c.getBitSet().xor(child.cluster.clusterRef.getBitSet());
				cs.add(child.cluster.clusterRef);
			}
			cluster = findCluster(c, this);
			if (c.equals(s) == false){
				STITreeCluster xc = new STITreeCluster(GlobalMaps.taxonIdentifier);
				xc.getBitSet().xor(c.getBitSet());
				xc.getBitSet().xor(s.getBitSet());
				cs.add(findCluster(xc, null).clusterRef);
			}
			if (cs.size() >= 3){
				AbstractPartition p = AbstractPartition.createPartition(cs);
				if (Polytree.this.partitions.containsKey(p)){
					partition = Polytree.this.partitions.get(p);
					partition.cnt++;
				}
				else partition = new PTPartition(p, this);
			}
		}
		PTCluster findCluster(STITreeCluster c, PTNode n){
			if (Polytree.this.clusters.containsKey(c)) {
				PTCluster cluster = Polytree.this.clusters.get(c);
				if (cluster.firstNode == null) cluster.firstNode = n;
				return cluster;
			}
			else return new PTCluster(c, n);
		}
		boolean isPartitionNode(){
			if (partition == null) return false;
			return partition.firstNode == this;
		}
		boolean isClusterNode(){
			return cluster.firstNode == this;
		}
		void addAllPartitions(){
			for (PTNode child: children){
				child.addAllPartitions();
			}
			if (isPartitionNode()){
				called = true;
				for (PTNode child: children){
					child.addAllClusters();
				}
			}
		}
		void addAllClusters(){
			if (called) return;
			if (isClusterNode()){
				called = true;
				for (PTNode child: children){
					child.addAllClusters();
				}
			}
			else {
				cluster.listUsed = true;
				PTNode n = cluster.firstNode;
				if (n == null || n.called) return;
				n.called = true;
				for (PTNode child: n.children){
					child.addAllClusters();
				}
			}
		}
		void buildQueue(){
			/*
			 * 1 - true = compute based on stack; false = get from list
			 * stack:
			 * value >> 5 - number of children
			 * 2 - true = store on stack
			 * 4 - true = store on list
			 * 8 - true = compute partition
			 * 16 - true = partition cnt>1 
			 * list:
			 * value >> 1 - the position on the list to fetch from
			 */
			for (PTNode child: children){
				child.buildQueue();
			}
			if (called){
				int v = (children.size() << 5) | 1;
				if (addToStack()) v = v | 2;
				if (addToList()) {
					v = v | 4;
					cluster.listPos = Polytree.this.listSize++;
				}
				if (isPartitionNode()){
					v = v | 8;
					if (partition.cnt > 1){
						v = v | 16;
						Polytree.this.queueBuilder.add(v);
						Polytree.this.queueBuilder.add(partition.cnt);
					}
					else Polytree.this.queueBuilder.add(v);
				}
				else Polytree.this.queueBuilder.add(v);
			}
			else {
				if (parent != null && parent.called) Polytree.this.queueBuilder.add(cluster.listPos << 1);
			}
		}
		boolean addToStack(){
			return called && parent != null && parent.called;
		}
		boolean addToList(){
			return called && isClusterNode() && cluster.listUsed;
		}
	}
	
	final class PTCluster{
		STITreeCluster clusterRef;
		PTNode firstNode;
		boolean listUsed = false;
		int listPos = -1;
		
		PTCluster(STITreeCluster c, PTNode n){
			clusterRef = c;
			firstNode = n;
			Polytree.this.clusters.put(c, this);
		}
		PTCluster(STITreeCluster c){
			clusterRef = c;
			firstNode = null;
			listUsed = true;
			listPos = Polytree.this.listSize++;
			Polytree.this.clusters.put(c, this);
		}
	}
	
	final class PTPartition{
		PTNode firstNode;
		int cnt = 1;
		
		PTPartition(AbstractPartition p, PTNode n){
			firstNode = n;
			Polytree.this.partitions.put(p, this);
		}
	}
	
	public static final class PTNative{
		private static boolean useNativeMethod = false;
		private static final int batchSize = 32;
		private static Polytree pt = null;
		
		static {
			try {
				System.loadLibrary("Astral");
				System.err.println("Using native AVX batch computing method.");
				useNativeMethod = true;
			}
			catch (Throwable e) {
				useNativeMethod = false;
				System.err.println("Fail to load native library; use Java default computing method.");
			}
		}
		private static native void cppInit(int n, int listSize, int[] q, long[][] c);
		private static native long cppCompute(long[] a, long[] b, long[] c);
		private static native void cppBatchCompute(long[] result, long[][] a, long[][] b, long[][] c);
		public static void compute(ArrayList<VertexPair> todolist) {
			System.err.println("number of jobs: " + todolist.size());
			long t = System.nanoTime();
			if (!useNativeMethod) {
				for (VertexPair p: todolist) {
					BitSet[] b = {
						p.cluster1.getCluster().getBitSet(),
						p.cluster2.getCluster().getBitSet(),
						p.both.getCluster().complementaryCluster().getBitSet()
					};
					p.weight = pt.WQWeightByTraversal(b);
				}
				Polytree.time += System.nanoTime() - t;
				return;
			}
			for (int i = 0; i < todolist.size(); i += batchSize) {
				int size = (todolist.size() - i > batchSize) ? batchSize : todolist.size() - i;
				long[] result = new long[size];
				long[][] a = new long[size][];
				long[][] b = new long[size][];
				long[][] c = new long[size][];
				for (int j = 0; j < size; j++) {
					a[j] = todolist.get(i + j).cluster1.getCluster().getBitSet().getArray();
					b[j] = todolist.get(i + j).cluster2.getCluster().getBitSet().getArray();
					c[j] = todolist.get(i + j).both.getCluster().complementaryCluster().getBitSet().getArray();
				}
				cppBatchCompute(result, a, b, c);
				for (int j = 0; j < size; j++) {
					todolist.get(i + j).weight = result[j];
				}
			}
			Polytree.time += System.nanoTime() - t;
		}
	}
	
	WQDataCollection dataCollection;
	HashMap<STITreeCluster, PTCluster> clusters = new HashMap<STITreeCluster, PTCluster>();	
	HashMap<AbstractPartition, PTPartition> partitions = new HashMap<AbstractPartition, PTPartition>();
	ArrayList<PTNode> nodeRoots = new ArrayList<PTNode>();
	ArrayList<Integer> queueBuilder = new ArrayList<Integer>();
	int[][] stack, list;
	int[] queue;
	int listSize = 0;
	long[] sx = new long[3], sxy = new long[3];
	int[] treeTotal = new int[3];
	long maxScore = 0;
	
	public Polytree(List<Tree> trees, WQDataCollection dataCollection){
		PTNative.pt = this;
		this.dataCollection = dataCollection;
		long t = System.currentTimeMillis();
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(i);
			new PTCluster(c);
		}
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			nodeRoots.add(buildTree(tr.getRoot(), tit.next()));
		}
		for (PTNode n: nodeRoots){
			n.addAllPartitions();
		}
		for (PTNode n: nodeRoots){
			queueBuilder.add(-1);
			n.buildQueue();
		}
		
		stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		list = new int[listSize][3];
		queue = mapToInt(queueBuilder);
		clusters = null;
		partitions = null;
		queueBuilder = null;

		STITreeCluster c = (new STITreeCluster(GlobalMaps.taxonIdentifier)).complementaryCluster();
		maxScore = WQWeightByTraversal(new Tripartition(c, c, c, false), null);
		System.err.println("Polytree max score: " + maxScore / 4);
		System.err.println("Polytree building time: " + (System.currentTimeMillis() - t) / 1000.0D + " seconds.");
		
		if (PTNative.useNativeMethod) {
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
	
	private PTNode buildTree(TNode node, STITreeCluster s){
		if (node.isLeaf()) return new PTNode(node);
		else {
			ArrayList<PTNode> cs = new ArrayList<PTNode>();
			for (TNode ch: node.getChildren()){
				cs.add(buildTree(ch, s));
			}
			return new PTNode(cs, s);
		}
	}

	public Long WQWeightByTraversal(Tripartition trip, CondensedTraversalWeightCalculator algorithm){
		if (trip.cluster1 == trip.cluster2) return computeUpperbound(trip.cluster1.getBitSet());
		long t = System.nanoTime();
		BitSet[] b = new BitSet[]{trip.cluster1.getBitSet(), trip.cluster2.getBitSet(), trip.cluster3.getBitSet()};
		return WQWeightByTraversal(b);
	}
	
	public Long WQWeightByTraversal(BitSet[] b){
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
		return weight;
	}
	
	public Long computeUpperbound(BitSet b){
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
