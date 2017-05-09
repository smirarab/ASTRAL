package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import phylonet.coalescent.WQWeightCalculator.CondensedTraversalWeightCalculator;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class Polytree {
	static long time = 0, time1 = 0, time2 = 0;
	
	static long F(int[] x, int[] y, int[] z){
		long a = x[0], b = x[1], c = x[2], d = y[0], e = y[1], f = y[2], g = z[0], h = z[1], i = z[2];
		return a * ( (a + e + i - 3)  * e * i + (a + f + h - 3)  * f * h )
			 + b * ((b + d + i - 3)  * d * i + (b + f + g - 3)  * f * g )
			 + c * ((c + d + h - 3)  * d * h + (c + e + g - 3)  * e * g );
	}
	
	HashMap<STITreeCluster, Integer> clusterID = new HashMap<STITreeCluster, Integer>();
	ArrayList<Integer> aDependerID = new ArrayList<Integer>();
	ArrayList<Integer> aDependeeID = new ArrayList<Integer>();
	ArrayList<Integer> aDependingFactor = new ArrayList<Integer>();
	HashMap<AbstractPartition, Integer> partitionCount = new HashMap<AbstractPartition, Integer>();
	ArrayList<Integer> aPartitionMultiplicity = new ArrayList<Integer>();
	ArrayList<Integer> aPartitionNumClusters = new ArrayList<Integer>();
	ArrayList<Integer> aPartitionClusterID = new ArrayList<Integer>();
	
	int[] dependerID, dependeeID, dependingFactor;
	int[] partitionMultiplicity, partitionNumClusters, partitionClusterID;
	private boolean randomResolveMultiInd; //TODO: Arbitrarily resolve a polytomy when this is true
											//       but all the children of the polytomy map to the same species. 
	
	int[][] overlap;

	public Polytree(List<Tree> trees, int polytomyMaxSize, boolean randomResolve){	
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(i);
			clusterID.put(c, i);
		}
		for (Tree tr: trees){
			STITreeCluster s = addSubtreeClusters(tr.getRoot(), null);
			addSubtreeClusters(tr.getRoot(), s);
			addSubtreePartitions(tr.getRoot(), s);
		}
		for (Entry<AbstractPartition, Integer> entry : partitionCount.entrySet()){
			STITreeCluster[] cs = entry.getKey().getClusters();
			if (polytomyMaxSize != 0 && cs.length > polytomyMaxSize) continue;
			aPartitionMultiplicity.add(entry.getValue());
			aPartitionNumClusters.add(cs.length);
			for (STITreeCluster c: cs){
				aPartitionClusterID.add(clusterID.get(c));
			}
		}
		dependerID = mapToInt(aDependerID);
		dependeeID = mapToInt(aDependeeID);	
		dependingFactor = mapToInt(aDependingFactor);
		aDependerID = null;
		aDependeeID = null;
		aDependingFactor = null;
		partitionMultiplicity = mapToInt(aPartitionMultiplicity);
		partitionNumClusters = mapToInt(aPartitionNumClusters);
		partitionClusterID = mapToInt(aPartitionClusterID);
		aPartitionMultiplicity = null;
		aPartitionNumClusters = null;
		aPartitionClusterID = null;
		overlap = new int[clusterID.size()][3];
	}
	
	private int[] mapToInt(List<Integer> list) {
		int[] ret = new int [list.size()]; 
		for (int i = 0; i < ret.length; i++)
			ret[i] = list.get(i);
		return ret;
	}
	
	private STITreeCluster addSubtreeClusters(TNode node, STITreeCluster s){
		STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
		ArrayList<STITreeCluster> children = new ArrayList<STITreeCluster>();
		if (node.isLeaf()){
			c.getBitSet().set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
		}
		else {
		    for (TNode child: node.getChildren()){
		    	STITreeCluster ch = addSubtreeClusters(child, s);
		    	children.add(ch);
		    	c.getBitSet().or(ch.getBitSet());
		    }
		}
		if (s == null){
			if (clusterID.containsKey(c) == false){
				clusterID.put(c, clusterID.size());
				for (STITreeCluster ch: children){
					aDependerID.add(clusterID.get(c));
					aDependeeID.add(clusterID.get(ch));
					aDependingFactor.add(1);
			    }
			}
		}
		else {
			if (node.isRoot() == false){
				STITreeCluster xc = new STITreeCluster(GlobalMaps.taxonIdentifier);
				xc.getBitSet().xor(c.getBitSet());
				xc.getBitSet().xor(s.getBitSet());
				
				if (clusterID.containsKey(xc) == false){
					clusterID.put(xc, clusterID.size());
					aDependerID.add(clusterID.get(xc));
					aDependeeID.add(clusterID.get(s));
					aDependingFactor.add(1);
					aDependerID.add(clusterID.get(xc));
					aDependeeID.add(clusterID.get(c));
					aDependingFactor.add(-1);
				}
			}
		}
		return c;
	}
	
	private STITreeCluster addSubtreePartitions(TNode node, STITreeCluster s){
		STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier);
		ArrayList<STITreeCluster> clusters = new ArrayList<STITreeCluster>();
		if (node.isLeaf()){
			c.getBitSet().set(GlobalMaps.taxonIdentifier.taxonId(node.getName()));
			return c;
		}
		for (TNode child: node.getChildren()){
	    	STITreeCluster ch = addSubtreePartitions(child, s);
	    	clusters.add(ch);
	    	c.getBitSet().or(ch.getBitSet());
	    }
	    if (node.isRoot()){
	    	if (node.getChildCount() < 3) return c;
	    }
	    else{
	    	STITreeCluster xc = new STITreeCluster(GlobalMaps.taxonIdentifier);
			xc.getBitSet().xor(c.getBitSet());
			xc.getBitSet().xor(s.getBitSet());
		    clusters.add(xc);
	    }
	    AbstractPartition p = AbstractPartition.createPartition(clusters);
	    if (partitionCount.containsKey(p) == false) partitionCount.put(p, 1);
	    else partitionCount.put(p, partitionCount.get(p) + 1);
	    return c;
	}
	
	public Long WQWeightByTraversal(Tripartition trip, CondensedTraversalWeightCalculator algorithm){
		long t = System.nanoTime();
		long weight = 0;
		long[] sx = new long[3], sxy = new long[3];
		int[] q;
		int[] from, to;
		BitSet[] b = new BitSet[]{trip.cluster1.getBitSet(), trip.cluster2.getBitSet(), trip.cluster3.getBitSet()};
		for (int i = 0, i_end = GlobalMaps.taxonIdentifier.taxonCount(); i < i_end; i++){
			overlap[i][0] = b[0].get(i) ? 1 : 0;
			overlap[i][1] = b[1].get(i) ? 1 : 0;
			overlap[i][2] = b[2].get(i) ? 1 : 0;
		}
		for (int i = GlobalMaps.taxonIdentifier.taxonCount(), i_end = overlap.length; i < i_end; i++){
			overlap[i][0] = 0;
			overlap[i][1] = 0;
			overlap[i][2] = 0;
		}
		time1 += System.nanoTime() - t;
		for (int i = 0, i_end = dependerID.length; i < i_end; i++){
			from = overlap[dependeeID[i]];
			to = overlap[dependerID[i]];
			if (dependingFactor[i] == 1){
				from[0] += to[0];
				from[1] += to[1];
				from[2] += to[2];
			}
			else {
				from[0] -= to[0];
				from[1] -= to[1];
				from[2] -= to[2];
			}
		}
		time2 += System.nanoTime() - t;
		for (int i = 0, j = 0, i_end = partitionMultiplicity.length; i < i_end; i++){
			if (partitionNumClusters[i] == 3){
				weight += F(overlap[partitionClusterID[j]], overlap[partitionClusterID[j + 1]], overlap[partitionClusterID[j + 2]])
						* partitionMultiplicity[i];
			}
			else{
				long tempWeight = 0;
				sx[0] = 0;
				sx[1] = 0;
				sx[2] = 0;
				sxy[0] = 0;
				sxy[1] = 0;
				sxy[2] = 0;
				for (int p = j, p_end = j + partitionNumClusters[i]; p < p_end; p++){
					q = overlap[partitionClusterID[p]];
					sx[0] += q[0];
					sx[1] += q[1];
					sx[2] += q[2];
					sxy[0] += q[1] * q[2];
					sxy[1] += q[2] * q[0];
					sxy[2] += q[0] * q[1];
				}
				for (int p = j, p_end = j + partitionNumClusters[i]; p < p_end; p++){
					q = overlap[partitionClusterID[p]];
					tempWeight += ((sx[1] - q[1]) * (sx[2] - q[2]) - sxy[0] + q[1] * q[2]) * q[0] * (q[0] - 1L)
						+ ((sx[2] - q[2]) * (sx[0] - q[0]) - sxy[1] + q[2] * q[0]) * q[1] * (q[1] - 1L)
						+ ((sx[0] - q[0]) * (sx[1] - q[1]) - sxy[2] + q[0] * q[1]) * q[2] * (q[2] - 1L);
				}
				weight += tempWeight * partitionMultiplicity[i];
			}
			j += partitionNumClusters[i];
		}
		time += System.nanoTime() - t;
		return weight;
	}
}
