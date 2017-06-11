package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import phylonet.coalescent.WQWeightCalculator.CondensedTraversalWeightCalculator;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class Polytree {
	static long time = 0;
	
	static long F(int[] x, int[] y, int[] z){
		long a = x[0], b = x[1], c = x[2], d = y[0], e = y[1], f = y[2], g = z[0], h = z[1], i = z[2];
		return a * ( (a + e + i - 3)  * e * i + (a + f + h - 3)  * f * h )
			 + b * ((b + d + i - 3)  * d * i + (b + f + g - 3)  * f * g )
			 + c * ((c + d + h - 3)  * d * h + (c + e + g - 3)  * e * g );
	}
	
	static long U(int[] x, int[] y, int[] z){
		long a = x[0], d = x[1], b = y[0], e = y[1], c = z[0], f = z[1], s = a + b + c - 3;
		return ((s - c + f) * b * f + (s - b + e) * e * c) * a + ((s - a + d) * d - 2 * s * a) * b * c;
	}
	
	final class PTNode{
		PTNode parent;
		PTCluster cluster, complementCluster;
		PTPartition partition;
		
		PTNode(PTNode p){
			parent = p;
		}
		PTNode(PTNode p, TNode n, STITreeCluster s){
			parent = p;
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			STITreeCluster xc = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(GlobalMaps.taxonIdentifier.taxonId(n.getName()));
			xc.getBitSet().xor(s.getBitSet());
			xc.getBitSet().xor(c.getBitSet());
			cluster = Polytree.this.clusters.get(c);
			if (Polytree.this.clusters.containsKey(xc)) complementCluster = Polytree.this.clusters.get(xc);
			else complementCluster = new PTCluster(xc, this);
			partition = null;
		}
		boolean needToComputeCluster(){
			if (cluster == null) return false;
			return (cluster.firstNode == this);
		}
		boolean needToComputeComplementCluster(){
			if (complementCluster == null) return false;
			return ((complementCluster.firstNode == this) && (complementCluster.used == true));
		}
		boolean needToComputePartition(){
			if (partition == null) return false;
			return (partition.firstNode == this);
		}
		boolean needToCompute(){
			return needToComputeCluster() || needToComputeComplementCluster() || needToComputePartition();
		}
		boolean needChildInfo(){
			return needToComputeCluster() || needToComputePartition();
		}
		boolean parentNeedChildInfo(){
			if (parent == null) return false;
			return parent.needChildInfo();
		}
		void addChidren(ArrayList<PTNode> children, STITreeCluster s){
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			STITreeCluster xc = new STITreeCluster(GlobalMaps.taxonIdentifier);
			ArrayList<STITreeCluster> cs = new ArrayList<STITreeCluster>();
			for (PTNode ch: children){
				STITreeCluster chc = ch.cluster.clusterRef;
				c.getBitSet().xor(chc.getBitSet());
				cs.add(chc);
			}
			if (Polytree.this.clusterRef.containsKey(c)) c = Polytree.this.clusterRef.get(c);
			xc.getBitSet().xor(s.getBitSet());
			xc.getBitSet().xor(c.getBitSet());
			if (Polytree.this.clusters.containsKey(c)) cluster = Polytree.this.clusters.get(c);
			else cluster = new PTCluster(c, this);
			if (xc.getBitSet().cardinality() != 0){
				if (Polytree.this.clusterRef.containsKey(xc)) xc = Polytree.this.clusterRef.get(xc);
				cs.add(xc);
				if (Polytree.this.clusters.containsKey(xc)) complementCluster = Polytree.this.clusters.get(xc);
				else complementCluster = new PTCluster(xc, this);
			}
			if (cs.size() >= 3){
				AbstractPartition p = AbstractPartition.createPartition(cs);
				if (Polytree.this.partitions.containsKey(p)){
					partition = Polytree.this.partitions.get(p);
					partition.cnt++;
				}
				else partition = new PTPartition(p, this);
			}
			else partition = null;
			if (needChildInfo()){
				for (PTNode ch: children){
					if (ch.needChildInfo() == false) ch.cluster.used = true;
				}
			}
		}
		void addComplementDependency(){
			if (needChildInfo() == false && complementCluster != null && complementCluster.used) cluster.used = true;
		}
	}
	
	final class PTCluster{
		STITreeCluster clusterRef;
		PTNode firstNode;
		boolean used = false;
		int listPos = -1;
		
		PTCluster(STITreeCluster c, PTNode n){
			clusterRef = c;
			firstNode = n;
			Polytree.this.clusters.put(c, this);
			Polytree.this.clusterRef.put(c, c);
			Polytree.this.clusterID.add(this);
		}
	}
	
	final class PTPartition{
		AbstractPartition partitionRef;
		PTNode firstNode;
		int cnt = 1;
		
		PTPartition(AbstractPartition p, PTNode n){
			partitionRef = p;
			firstNode = n;
		}
	}
	
	WQDataCollection dataCollection;
	HashMap<STITreeCluster, STITreeCluster> clusterRef = new HashMap<STITreeCluster, STITreeCluster>();
	HashMap<STITreeCluster, PTCluster> clusters = new HashMap<STITreeCluster, PTCluster>();	
	HashMap<AbstractPartition, PTPartition> partitions = new HashMap<AbstractPartition, PTPartition>();
	ArrayList<PTNode> nodes = new ArrayList<PTNode>();
	ArrayList<PTCluster> clusterID = new ArrayList<PTCluster>();
	ArrayList<Integer> queueBuilder = new ArrayList<Integer>();
	int[][] stack, list;
	int[] queue;
	int nodeIDCounter = 0;
	int[] sx = new int[3], sxy = new int[3], treeTotal = new int[3];
	long maxScore = 0;
	
	int debugClusterPos = 0;
	
	public Polytree(List<Tree> trees, WQDataCollection dataCollection){
		this.dataCollection = dataCollection;
		long t = System.currentTimeMillis();
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(i);
			(new PTCluster(c, null)).used = true;
			
			debugClusterPos++;
		}
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			traversal1(tr.getRoot(), tit.next(), null);
		}
		int listSize = 0;
		for (PTNode n: nodes){
			n.addComplementDependency();
		}
		for (PTCluster c: clusterID){
			if (c.used) c.listPos = listSize++;
		}
		tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			queueBuilder.add(-1);
			traversal2(tr.getRoot());
		}
		
		stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		list = new int[listSize][3];
		queue = mapToInt(queueBuilder);
		
		nodes = null;
		clusterRef = null;
		clusters = null;
		clusterID = null;
		partitions = null;
		queueBuilder = null;
		
		STITreeCluster c = (new STITreeCluster(GlobalMaps.taxonIdentifier)).complementaryCluster();
		maxScore = WQWeightByTraversal(new Tripartition(c, c, c, false), null);
		System.err.println("Polytree max score: " + maxScore / 4);
		System.err.println("Polytree building time: " + (System.currentTimeMillis() - t) / 1000.0D + " seconds.");
	}
	
	private int[] mapToInt(List<Integer> list) {
		int[] ret = new int [list.size()]; 
		for (int i = 0; i < ret.length; i++)
			ret[i] = list.get(i);
		return ret;
	}
	
	private PTNode traversal1(TNode node, STITreeCluster s, PTNode parent){
		PTNode n = null;
		if (node.isLeaf()){
			n = new PTNode(parent, node, s);
		}
		else {
			n = new PTNode(parent);
			ArrayList<PTNode> cs = new ArrayList<PTNode>();
			for (TNode ch: node.getChildren()){
				cs.add(traversal1(ch, s, n));
			}
			n.addChidren(cs, s);
		}
		nodes.add(n);
		return n;
	}

	private void traversal2(TNode node){
		for (TNode child: node.getChildren()){
		    traversal2(child);
		}
		PTNode n = nodes.get(nodeIDCounter);
		/*
		 * 1 - true = compute based on stack; false = get from list
		 * stack:
		 * 2 - true = store on stack
		 * 3 - true = store on list
		 * 4 - true = store complement on list
		 * 5 - true = compute partition
		 * 6 - true = partition cnt>1
		 * list:
		 * 2 - true = store on stack
		 * 3 - true = store complement on list 
		 */
		if (n.needToCompute() || n.parentNeedChildInfo()){
			int cmd;
			if (n.needChildInfo()){
				cmd = ((node.getChildCount() << 6) | 1);
				if (n.parentNeedChildInfo()) cmd = cmd | 2;
				if (n.needToComputeCluster() && n.cluster.used) cmd = cmd | 4;
				if (n.needToComputeComplementCluster() && n.complementCluster.used) cmd = cmd | 8;
				if (n.needToComputePartition()){
					cmd = cmd | 16;
					if (n.partition.cnt > 1){
						cmd = cmd | 32;
						queueBuilder.add(cmd);
						queueBuilder.add(n.partition.cnt);
					}
					else queueBuilder.add(cmd);
				}
				else queueBuilder.add(cmd);
			}
			else {
				cmd = n.cluster.listPos << 3;
				if (n.parentNeedChildInfo()) cmd = cmd | 2;
				if (n.needToComputeComplementCluster()) cmd = cmd | 4;
				queueBuilder.add(cmd);
			}
		}
		nodeIDCounter++;
	}

	public Long WQWeightByTraversal(Tripartition trip, CondensedTraversalWeightCalculator algorithm){
		if (trip.cluster1 == trip.cluster2) return computeUpperbound(trip.cluster1.getBitSet());
		long t = System.nanoTime();
		long weight = 0;
		int stackEnd = 0, listEnd = GlobalMaps.taxonIdentifier.taxonCount();
		BitSet[] b = new BitSet[]{trip.cluster1.getBitSet(), trip.cluster2.getBitSet(), trip.cluster3.getBitSet()};
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
				int numChildren = cmd >> 6;
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
				if ((cmd & 16) != 0){
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
					if ((cmd & 32) != 0) weight += tempWeight * queue[++i];
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
				if ((cmd & 8) != 0){
					q = list[listEnd++];
					q[0] = p[0];
					q[1] = p[1];
					q[2] = p[2];
				}
			}
			else {
				int[] p = list[cmd >> 3], q;
				if ((cmd & 2) != 0){
					q = stack[stackEnd++];
					q[0] = p[0];
					q[1] = p[1];
					q[2] = p[2];
				}
				if ((cmd & 4) != 0){
					q = list[listEnd++];
					q[0] = treeTotal[0] - p[0];
					q[1] = treeTotal[1] - p[1];
					q[2] = treeTotal[2] - p[2];
				}
			}
		}
		time += System.nanoTime() - t;
		return weight;
	}
	
	public Long computeUpperbound(BitSet b){
		long t = System.nanoTime();
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
				int numChildren = cmd >> 6;
				long tempWeight = 0;
				int[] p = stack[stackEnd], q;
				p[0] = treeTotal[0];
				p[1] = treeTotal[1];
				for (int j = stackEnd - numChildren; j < stackEnd; j++){
					q = stack[j];
					p[0] -= q[0];
					p[1] -= q[1];
				}
				if ((cmd & 16) != 0){
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
					if ((cmd & 32) != 0) weight += tempWeight * queue[++i];
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
				if ((cmd & 8) != 0){
					q = list[listEnd++];
					q[0] = p[0];
					q[1] = p[1];
				}
			}
			else {
				int[] p = list[cmd >> 3], q;
				if ((cmd & 2) != 0){
					q = stack[stackEnd++];
					q[0] = p[0];
					q[1] = p[1];
				}
				if ((cmd & 4) != 0){
					q = list[listEnd++];
					q[0] = treeTotal[0] - p[0];
					q[1] = treeTotal[1] - p[1];
				}
			}
		}
		time += System.nanoTime() - t;
		return weight;
	}
}
