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
	
	WQDataCollection dataCollection;
	HashMap<STITreeCluster, STITreeCluster> clusterRef = new HashMap<STITreeCluster, STITreeCluster>();
	HashMap<STITreeCluster, Integer> clusterID = new HashMap<STITreeCluster, Integer>();
	HashMap<AbstractPartition, Integer> partitionID = new HashMap<AbstractPartition, Integer>();
	ArrayList<Integer> nodeParent = new ArrayList<Integer>();
	ArrayList<Integer> nodeCluster = new ArrayList<Integer>();
	ArrayList<Integer> nodePartition = new ArrayList<Integer>();
	ArrayList<Integer> clusterNode = new ArrayList<Integer>();
	ArrayList<Integer> clusterUsage = new ArrayList<Integer>();
	ArrayList<Integer> clusterPos = new ArrayList<Integer>();
	ArrayList<Integer> partitionNode = new ArrayList<Integer>();
	ArrayList<Integer> partitionUsage = new ArrayList<Integer>();
	ArrayList<Integer> queueBuilder = new ArrayList<Integer>();
	int[][] stack, list;
	int[] queue;
	int nodeIDCounter = 0;
	int[] sx = new int[3], sxy = new int[3], treeTotal = new int[3];
	long maxScore = 0;
	
	public Polytree(List<Tree> trees, WQDataCollection dataCollection){
		this.dataCollection = dataCollection;
		long t = System.currentTimeMillis();
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(i);
			clusterID.put(c, i);
			clusterNode.add(-1);
			clusterUsage.add(1);
		}
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			traversal1(tr.getRoot(), tit.next());
		}
		int listSize = 0;
		for (int usage: clusterUsage){
			clusterPos.add((usage > 0) ? listSize++ : -1);
		}
		
		tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			queueBuilder.add(-1);
			traversal2(tr.getRoot());
		}
		
		stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		list = new int[listSize][3];
		queue = mapToInt(queueBuilder);
		
		clusterID = null;
		partitionID = null;
		nodeParent = null;
		nodeCluster = null;
		nodePartition = null;
		clusterNode = null;
		clusterUsage = null;
		clusterPos = null;
		partitionNode = null;
		partitionUsage = null;
		queueBuilder = null;
		
		STITreeCluster c = (new STITreeCluster()).complementaryCluster();
		maxScore = WQWeightByTraversal(new Tripartition(c, c, c, false), null) / 6;
		System.err.println("Polytree max score: " + maxScore / 4);
		System.err.println("Polytree building time: " + (System.currentTimeMillis() - t) / 1000.0D + " seconds.");
	}
	
	private int[] mapToInt(List<Integer> list) {
		int[] ret = new int [list.size()]; 
		for (int i = 0; i < ret.length; i++)
			ret[i] = list.get(i);
		return ret;
	}
	
	private STITreeCluster traversal1(TNode node, STITreeCluster s){
		STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
		if (node.isLeaf()){
			c.getBitSet().set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
			nodeParent.add(-1);
			nodeCluster.add(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
			nodePartition.add(-1);
		}
		else {
			ArrayList<STITreeCluster> cs = new ArrayList<STITreeCluster>();
			ArrayList<Integer> childrenID = new ArrayList<Integer>();
			for (TNode child: node.getChildren()){
		    	STITreeCluster ch = traversal1(child, s);
		    	cs.add(ch);
		    	c.getBitSet().or(ch.getBitSet());
		    	childrenID.add(nodeParent.size() - 1);
		    }
			
		    nodeParent.add(-1);
		    if (clusterRef.containsKey(c)){
		    	c = clusterRef.get(c);
		    	nodeCluster.add(clusterID.get(c));
		    }
		    else{
		    	clusterRef.put(c, c);
		    	clusterID.put(c, clusterID.size());
		    	nodeCluster.add(clusterID.size() - 1);
		    	clusterNode.add(nodeParent.size() - 1);
				clusterUsage.add(0);
		    }
		    if (c.equals(s) == false){
		    	STITreeCluster xc = new STITreeCluster(GlobalMaps.taxonIdentifier);
			    xc.getBitSet().xor(c.getBitSet());
			    xc.getBitSet().xor(s.getBitSet());
			    cs.add(xc);	
		    }
		    if (cs.size() >= 3){
		    	AbstractPartition p = AbstractPartition.createPartition(cs);
				if (partitionID.containsKey(p)){
					Integer pID = partitionID.get(p);
					nodePartition.add(pID);
					partitionUsage.set(pID, partitionUsage.get(pID) + 1);
				}
				else{
					partitionID.put(p, partitionID.size());
					nodePartition.add(partitionID.size() - 1);
			    	partitionNode.add(nodeParent.size() - 1);
			    	partitionUsage.add(1);
				}
		    }
		    else{
		    	nodePartition.add(-1);
		    }
		    
		    for (int childID: childrenID){
		    	nodeParent.set(childID, nodeParent.size() - 1);
		    	if ((nodePartition.get(nodeParent.size() - 1) != -1 && (partitionNode.get(nodePartition.get(nodeParent.size() - 1)) == nodeParent.size() - 1))
		    		&& (nodePartition.get(childID) == -1 || partitionNode.get(nodePartition.get(childID)) != childID)){
		    		Integer clusterID = nodeCluster.get(childID);
		    		clusterUsage.set(clusterID, clusterUsage.get(clusterID) + 1);
		    	}
		    }
		    
		}
		return c;
	}

	private void traversal2(TNode node){
		ArrayList<Integer> childrenID = new ArrayList<Integer>();
		for (TNode child: node.getChildren()){
		    traversal2(child);
		    childrenID.add(nodeIDCounter - 1);
		}
		boolean pos0 = (nodePartition.get(nodeIDCounter) != -1 && partitionNode.get(nodePartition.get(nodeIDCounter)) == nodeIDCounter);
		boolean pos1 = (nodeParent.get(nodeIDCounter) != -1 && nodePartition.get(nodeParent.get(nodeIDCounter)) != -1
			&& partitionNode.get(nodePartition.get(nodeParent.get(nodeIDCounter))).intValue() == nodeParent.get(nodeIDCounter).intValue());
		if (pos0 == true || pos1 == true){
			if (pos0 == true){
				int cmd = (node.getChildCount() << 4) | 1;
				if (pos1 == true) cmd = cmd | 2;
				if (nodeCluster.get(nodeIDCounter) != -1 && clusterNode.get(nodeCluster.get(nodeIDCounter)) == nodeIDCounter
					&& clusterUsage.get(nodeCluster.get(nodeIDCounter)) > 0){
					cmd = cmd | 4;
				}
				if (nodePartition.get(nodeIDCounter) != -1 && partitionUsage.get(nodePartition.get(nodeIDCounter)) > 1){
					cmd = cmd | 8;
					queueBuilder.add(cmd);
					queueBuilder.add(partitionUsage.get(nodePartition.get(nodeIDCounter)));
				}
				else {
					queueBuilder.add(cmd);
				}
			}
			else {
				queueBuilder.add(clusterPos.get(nodeCluster.get(nodeIDCounter)) << 1);
			}
		}
		nodeIDCounter++;
	}

	public Long WQWeightByTraversal(Tripartition trip, CondensedTraversalWeightCalculator algorithm){
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
				int numChildren = cmd >> 4;
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
				if ((cmd & 8) != 0) weight += tempWeight * queue[++i];
				else weight += tempWeight;
			}
			else {
				int pos = cmd >> 1;
				stack[stackEnd][0] = list[pos][0];
				stack[stackEnd][1] = list[pos][1];
				stack[stackEnd][2] = list[pos][2];
				stackEnd++;
			}
		}
		time += System.nanoTime() - t;
		return weight;
	}
}
