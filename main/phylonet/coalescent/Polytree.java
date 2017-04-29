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

public class Polytree {
	
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
	private int plytomyMaxSize; //TODO: Change code to take this into account. 
	private boolean randomResolveMultiInd; //TODO: Arbitrarily resolve a polytomy when this is true
											//       but all the children of the polytomy map to the same species. 
	
	public Polytree(List<Tree> trees, int plytomyMaxSize, boolean randomResolve){
		
		this.plytomyMaxSize = plytomyMaxSize; 
		this.randomResolveMultiInd  = randomResolve;
		
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
			aPartitionMultiplicity.add(entry.getValue());
			aPartitionNumClusters.add(cs.length);
			for (STITreeCluster c: cs){
				aPartitionClusterID.add(clusterID.get(c));
			}
		}
		mapToInt(dependerID, aDependerID);
		mapToInt(dependeeID, aDependeeID);
		
		mapToInt(dependingFactor, aDependingFactor);
		aDependerID = null;
		aDependeeID = null;
		aDependingFactor = null;
		mapToInt(partitionMultiplicity , aPartitionMultiplicity);
		mapToInt(partitionNumClusters , aPartitionNumClusters);
		mapToInt(partitionClusterID , aPartitionClusterID);
		aPartitionMultiplicity = null;
		aPartitionNumClusters = null;
		aPartitionClusterID = null;
	}
	
	private void mapToInt(int [] ret, List<Integer> list) {
		ret = new int [list.size()]; 
		for (int i = 0; i < ret.length; i++ )
			ret[i] = list.get(i);
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
		int[][] overlap = new int[clusterID.size()][3];
		long weight = 0;
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			overlap[i][0] = trip.cluster1.getBitSet().get(i) ? 1 : 0;
			overlap[i][1] = trip.cluster2.getBitSet().get(i) ? 1 : 0;
			overlap[i][2] = trip.cluster3.getBitSet().get(i) ? 1 : 0;
		}
		for (int i = 0; i < dependerID.length; i++){
			overlap[dependerID[i]][0] += overlap[dependeeID[i]][0] * dependingFactor[i];
			overlap[dependerID[i]][1] += overlap[dependeeID[i]][1] * dependingFactor[i];
			overlap[dependerID[i]][2] += overlap[dependeeID[i]][2] * dependingFactor[i];
		}
		for (int i = 0, j = 0; i < partitionMultiplicity.length; i++){
			if (partitionNumClusters[i] == 3){
				weight += algorithm.F(overlap[partitionClusterID[j]], overlap[partitionClusterID[j + 1]], overlap[partitionClusterID[j + 2]])
						* partitionMultiplicity[i];
			}
			else{
				for (int p = j; p < j + partitionNumClusters[i]; p++){
					for (int q = p + 1; q < j + partitionNumClusters[i]; q++){
						for (int r = p + 1; r < j + partitionNumClusters[i]; r++){
							weight += algorithm.F(overlap[partitionClusterID[p]], overlap[partitionClusterID[q]], overlap[partitionClusterID[r]])
									* partitionMultiplicity[i];
						}
					}
				}
			}
			j += partitionNumClusters[i];
		}
		
		/*
			// The following case is relevant only for polytomies.

			int [] nzc = {0,0,0};
			int [] newSides = {0,0,0};
			for (int side = 0; side < 3; side++) {
				for (int i = top - 1; i >= top + gtb; i--) {
					if (stack[i][side] > 0) {
						newSides[side] += stack[i][side];
						overlap[nzc[side]][side] = stack[i][side]; 
						overlapind[nzc[side]++][side] = i;
					}
				}
				stack[top][side] = allsides[side] - newSides[side];

				if (stack[top][side] > 0) {
					overlap[nzc[side]][side] = stack[top][side]; 
					overlapind[nzc[side]++][side] = top;
				}
				stack[top + gtb][side] = newSides[side];
			}

			for (int i = nzc[0] - 1; i >= 0; i--) {
				for (int j = nzc[1] - 1; j >= 0; j--) {
					if (overlapind[i][0] != overlapind[j][1])
						for (int k = nzc[2] - 1; k >= 0; k--) {
							if ((overlapind[i][0] != overlapind[k][2]) &&
								(overlapind[j][1] != overlapind[k][2]))
								weight += F(overlap[i][0], overlap[j][1], overlap[k][2]);
						}
				}
			}
		*/
		return weight;
	}
}
