package phylonet.coalescent;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import phylonet.coalescent.MGDInference_DP.Vertex;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

public class DuplicationWeightCounter {
	

	HashMap<STBipartition, Integer> weights = new HashMap<STBipartition, Integer>();
	
	String [] taxa;
	String [] gtTaxa;
	String [] stTaxa;

	private List<Set<STBipartition>> geneTreeSTBBySize;
	private List<Set<STBipartition>> X;
	private HashMap<STITreeCluster, Set<STBipartition>> clusterToSTBs = 
		new HashMap<STITreeCluster, Set<STBipartition>>();

	private Map<STBipartition,Integer> STBCountInGeneTrees;	
	
	private STITreeCluster makeCluster(String rep) {
		String [] reptaxa = rep.split("");
		STITreeCluster cluster = new STITreeCluster(taxa);
		for (int i = 0; i < reptaxa.length; i++) {
			String string = reptaxa[i];
			if (string != "") 
				cluster.addLeaf(string);
		}
		return cluster;
	}
	
	void computeTreeSTBipartitions(List<Tree> trees, String[] leaves, 
			Map<String, String> taxonMap, Map<Integer, Set<Vertex>> clusters) {

		geneTreeSTBBySize = new ArrayList<Set<STBipartition>>(leaves.length);
		X = new ArrayList<Set<STBipartition>>(leaves.length);
		
		for (int i = 0; i <= leaves.length; i++) {
			geneTreeSTBBySize.add(new HashSet<STBipartition>());
			X.add(new HashSet<STBipartition>());
		}
		STBCountInGeneTrees = new HashMap<STBipartition, Integer>();
		
		for (Tree tr : trees) {
			Map<TNode,STITreeCluster> map = new HashMap<TNode, STITreeCluster>();
			for (Iterator<TNode> nodeIt = tr.postTraverse()
					.iterator(); nodeIt.hasNext();) {
				TNode node = nodeIt.next();
		        BitSet bs = new BitSet(leaves.length);
	            if(node.isLeaf())
	            {
	                int i = 0;
	                String nodeName = node.getName();
	                if (taxonMap != null) {
	                	nodeName = taxonMap.get(nodeName);
	                }
	                for(i = 0; i < leaves.length; i++)
	                    if(nodeName.equals(leaves[i]))
	                        break;

	                bs.set(i);
	                STITreeCluster tb = new STITreeCluster(leaves);
	                tb.setCluster((BitSet) bs.clone());	                
	                map.put(node, tb);
	            } else {
	                int childCount = node.getChildCount();
	                BitSet childbslist[] = new BitSet[childCount];
	                int index = 0;
	                for(Iterator iterator3 = node.getChildren().iterator(); iterator3.hasNext();)
	                {
	                    TNode child = (TNode)iterator3.next();
	                    BitSet childCluster = map.get(child).getCluster();
	                    bs.or(childCluster);
	                    childbslist[index++] = childCluster;
	                }
	                
	                if (index > 2) {
	                	throw new RuntimeException("None bifurcating tree: "+
	                			tr+ "\n" +
	                			node);
	                }
	                
	                STITreeCluster cluster = new STITreeCluster(leaves);
	                cluster.setCluster((BitSet) bs.clone());
	                
	                
	                map.put(node, cluster);

	                STITreeCluster l_cluster = new STITreeCluster(leaves);
	                l_cluster.setCluster((BitSet) childbslist[0].clone());
	                
	                STITreeCluster r_cluster = new STITreeCluster(leaves);
	                r_cluster.setCluster((BitSet) childbslist[1].clone());
	                
	                if (l_cluster.isDisjoint(r_cluster)) {
	                	
	                	STBipartition stb = new STBipartition(r_cluster, l_cluster);
	                
	                	geneTreeSTBBySize.get(cluster.getClusterSize()).add(stb);
	                	//X.get(cluster.getClusterSize()).add(stb);
	                	addSTBToX(clusters, l_cluster, r_cluster);
	                	
	                	STBCountInGeneTrees.put(stb, 
	                			STBCountInGeneTrees.containsKey(stb)? 
	                					STBCountInGeneTrees.get(stb)+1: 1);
	                } else {
	                	BitSet and  = (BitSet) l_cluster.getCluster().clone();
	                	and.and(r_cluster.getCluster());
	                	
	                	BitSet l_Minus_r  = (BitSet) and.clone();
	                	l_Minus_r.xor(l_cluster.getCluster());
	                	STITreeCluster lmr = new STITreeCluster(stTaxa);
	                	lmr.setCluster(l_Minus_r);

	                	BitSet r_Minus_l  = (BitSet) and.clone();
	                	r_Minus_l.xor(r_cluster.getCluster());
	                	STITreeCluster rml = new STITreeCluster(stTaxa);
	                	rml.setCluster(r_Minus_l);	                	
	                	
	                	if (rml.getClusterSize() > 0) {
	                		addSTBToX(clusters,l_cluster, rml);
	                	}
	                	if (lmr.getClusterSize() > 0) {
	                		addSTBToX(clusters,r_cluster, lmr);
	                	}
	                }
	            }	           
			}

		}
		System.out.println("STBs in gene trees: " + STBCountInGeneTrees.keySet().size());
	}
	
	void addExtraBipartitions(Map<Integer, Set<Vertex>> clusters, String[] stTaxa) {
		// Set X = set of STBs observed in gene trees
		/*List<Vertex> allclusters = new ArrayList<Vertex>();
		for (Set<Vertex> list : clusters.values()) {
			allclusters.addAll(list);
		}
		for (Vertex vertex1 : allclusters) {
			for (Vertex vertex2 : allclusters) {
				if (! vertex1._cluster.isDisjoint(vertex2._cluster)) {
					continue;
				}
				addSTBToX(clusters, vertex1._cluster, vertex2._cluster);
			}
			
		}
*/	}

	private void addSTBToX(Map<Integer, Set<Vertex>> clusters, STITreeCluster left,
			STITreeCluster right) {
		STBipartition stb = new STBipartition(left, right);
		int size = left.getClusterSize()+ right.getClusterSize();
		if (!X.get(size).contains(stb)){
			X.get(size).add(stb);
			STBCountInGeneTrees.put(stb, 0);
			STITreeCluster c = new STITreeCluster(left);
			c = c.merge(right);
			Vertex nv = new Vertex();
			nv._cluster = c;
			if (!clusters.containsKey(size)) {
				clusters.put(size, new HashSet<Vertex>());
			}
			if (! clusters.get(size).contains(nv)){						
				nv._min_cost = -1;
				clusters.get(size).add(nv);
			}
			Set<STBipartition> stbs = clusterToSTBs.get(c);
			stbs = (stbs== null)? new HashSet<STBipartition>() : stbs;
			stbs.add(stb);
			clusterToSTBs.put(c, stbs);
		}
	}

	void calculateWeights(String[] leaves) {
		weights = new HashMap<STBipartition, Integer>();
		weights.putAll(STBCountInGeneTrees);						
		
		for (int i = leaves.length; i > 1; i--) {
			Set<STBipartition> biggerSTBs = X.get(i);
			for (STBipartition biggerSTB : biggerSTBs) {
				for (int j = i - 1; j > 1; j--) {
					Set<STBipartition> smallerSTBs = geneTreeSTBBySize.get(j);
					for (STBipartition smallerSTB : smallerSTBs) {
						if (smallerSTB.isDominatedBy(biggerSTB)) {
							weights.put(biggerSTB, weights.get(biggerSTB)+weights.get(smallerSTB));
						}
					}
				}
			}
		}
	}

	/*public void calculateWeights(String[] leaves, List<Tree> trees) {
	
		weights = new HashMap<STBipartition, Integer>();
		weights.putAll(STBCountInGeneTrees);						
		
		for (int i =0 ; i < leaves.length; i++) {
			Set<STBipartition> gtSTBs = geneTreeSTBBySize.get(i);
			for (STBipartition gtSTB : gtSTBs) {
				STITreeCluster c = new STITreeCluster(gtSTB.cluster1);
				c.merge(gtSTB.cluster2);
				Set<STBipartition> stSTBs = clusterToSTBs.get(c);
			}
		}
		
		for (int i = leaves.length; i > 1; i--) {
			Set<STBipartition> biggerSTBs = X.get(i);
			for (STBipartition biggerSTB : biggerSTBs) {
				for (int j = i - 1; j > 1; j--) {
					Set<STBipartition> smallerSTBs = geneTreeSTBBySize.get(j);
					for (STBipartition smallerSTB : smallerSTBs) {
						if (smallerSTB.isDominatedBy(biggerSTB)) {
							weights.put(biggerSTB, weights.get(biggerSTB)+weights.get(smallerSTB));
						}
					}
				}
			}
		}
	}
*/
	public DuplicationWeightCounter(String[] stTaxa){
		taxa = stTaxa;		
	}

	public DuplicationWeightCounter(String[] gtTaxa, String[] stTaxa){
		this.gtTaxa = gtTaxa;		
		this.stTaxa = stTaxa;
	}
	
	public boolean BiPartitionExists(STITreeCluster cluster1, STITreeCluster cluster2) {
		STBipartition bi = new STBipartition(cluster1,cluster2);
		return weights.containsKey(bi);
	}
	public int getBiPartitionDPWeight(STITreeCluster cluster1, STITreeCluster cluster2) {
		
		STBipartition bi = new STBipartition(cluster1,cluster2);
		//boolean a = bi.equals(bc);
		if (!weights.containsKey(bi)){
			System.out.println(bi);	
			throw new RuntimeException("Weight not found");
		}
		return weights.get(bi);
	}

	public Set<STBipartition> getClusterBiPartitions(STITreeCluster cluster) {
						
		return clusterToSTBs.get(cluster);
	}
	static class STBipartition {
		
		STITreeCluster cluster1;
		STITreeCluster cluster2;
		
		public STBipartition(STITreeCluster c1, STITreeCluster c2) {
			cluster1 = c1;
			cluster2 = c2;
		}
		@Override
		public boolean equals(Object obj) {
			STBipartition stb2 = (STBipartition) obj; 
			return ((stb2.cluster1.equals(this.cluster1) && stb2.cluster2.equals(this.cluster2)) ||
					stb2.cluster2.equals(this.cluster1) && stb2.cluster1.equals(this.cluster2));					
		}
		@Override
		public int hashCode() {
			return cluster1.hashCode() + cluster2.hashCode();
		}
		@Override
		public String toString() {		
			return cluster1.toString()+"|"+cluster2.toString();
		}
		
		public boolean isDominatedBy(STBipartition dominant) {
			return (dominant.cluster1.containsCluster(this.cluster1) && dominant.cluster2.containsCluster(this.cluster2)) ||
					(dominant.cluster2.containsCluster(this.cluster1) && dominant.cluster1.containsCluster(this.cluster2));
		}
		
	}
	
}