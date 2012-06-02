package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import javax.management.RuntimeErrorException;



import phylonet.coalescent.MGDInference_DP.Vertex;
import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class DuplicationWeightCounter {
	

	HashMap<STBipartition, Integer> weights;
	
	String [] gtTaxa;
	String [] stTaxa;

	private List<Set<STBipartition>> geneTreeSTBBySize;
	
	//private List<Set<STBipartition>> X;
	
	private HashMap<STITreeCluster, Set<STBipartition>> X;
	
	//private HashMap<TNode, STBipartition> gtNodeToSTBs;

	private Map<STBipartition,Integer> STBCountInGeneTrees;
	
	//private Map<Tree,List<STBipartition>> treeSTBs;

	//private List<Set<STBipartition>> goodSTBs;

	private boolean rooted;
	
	private HashMap<STBipartition, Set<STBipartition>> alreadyProcessed = new HashMap<STBipartition, Set<STBipartition>>();
	
	private void addToClusters (Map<Integer, Set<Vertex>> clusters, STITreeCluster c, int size) {
		Vertex nv = new Vertex();
		nv._cluster = c;
		if (!clusters.get(size).contains(nv)){						
			nv._min_cost = -1;
			nv._el_num = -1; 			
			clusters.get(size).add(nv);
		}
	}
	

	void addExtraBipartitionsByInput(Map<Integer, Set<Vertex>> clusters,
			List<Tree> trees,Map<String, String> taxonMap , boolean extraTreeRooted) {

		int sigmaN = 0;
		int k = trees.size();
		String[] leaves = stTaxa;
		int n = leaves.length;
		
		STITreeCluster all = clusters.get(n).iterator().next()._cluster;
				
		for (Tree tr : trees) {			
			String[] treeLeaves = tr.getLeaves();			
			STITreeCluster treeAll = new STITreeCluster(leaves);
			for (int i = 0; i < treeLeaves.length; i++) {
				treeAll.addLeaf(treeLeaves[i]);
			}		
			sigmaN += tr.getLeafCount() - 1;
			Map<TNode,STITreeCluster> map = new HashMap<TNode, STITreeCluster>(n);
			for (Iterator<TNode> nodeIt = tr.postTraverse()
					.iterator(); nodeIt.hasNext();) {
				TNode node = nodeIt.next();		        
	            if(node.isLeaf())
	            {
	                String nodeName = node.getName();
	                if (taxonMap != null) {
	                	nodeName = taxonMap.get(nodeName);
	                }

	                STITreeCluster tb = new STITreeCluster(leaves);
	                tb.addLeaf(nodeName);	               

        			map.put(node, tb);
	                
	                if (!extraTreeRooted) {
		                addSTBToX(clusters, tb, treeComplementary(treeAll, tb), null,node);
	                }
	            } else {
	                int childCount = node.getChildCount();
	                STITreeCluster childbslist[] = new STITreeCluster[childCount];
			        BitSet bs = new BitSet(leaves.length);
	                int index = 0;
	                for(Iterator iterator3 = node.getChildren().iterator(); iterator3.hasNext();)
	                {
	                    TNode child = (TNode)iterator3.next();
	                    childbslist[index++] = map.get(child);
	                    bs.or(map.get(child).getCluster());
	                }
	                	                		                
	                STITreeCluster cluster = new STITreeCluster(leaves);
	                cluster.setCluster((BitSet) bs.clone());	
	                
	                int size = cluster.getClusterSize();
	                
	                map.put(node, cluster);	                
	                
	                
	                if (extraTreeRooted) {

	                	if (index > 2) {	                	
		                	throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);
		                }

		                STITreeCluster l_cluster = childbslist[0];
		                
		                STITreeCluster r_cluster = childbslist[1];
		                		                
		                addSTBToX(clusters, l_cluster, r_cluster, cluster,node);
	                } else {		                	
	                	if (childCount == 2) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster r_cluster = childbslist[1];
			                
			                STITreeCluster allMinuslAndr_cluster = treeComplementary(treeAll, cluster);
			                		                
			                STITreeCluster lAndr_cluster = cluster;
			                
			                // add Vertex STBs
			                addSTBToX(clusters, l_cluster, r_cluster, cluster,node);
			                if (allMinuslAndr_cluster.getClusterSize() != 0) {
			                	addSTBToX(clusters, r_cluster, allMinuslAndr_cluster, null,node);
			                	addSTBToX(clusters, l_cluster, allMinuslAndr_cluster, null,node);
				                addSTBToX(clusters, lAndr_cluster, allMinuslAndr_cluster,  all,node);
			                }			                

	                	} else if (childCount == 3 && node.isRoot()) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster m_cluster =  childbslist[1];
			                
			                STITreeCluster r_cluster = childbslist[2];

			                addSTBToX(clusters, l_cluster, r_cluster, null,node);
			                addSTBToX(clusters, r_cluster, m_cluster, null,node);
			                addSTBToX(clusters, l_cluster, m_cluster, null,node);
	                	} else {
	                		throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);	
	                	}	                	
	                }
	            }	           
			}

		}
		int s = 0;
		for (Integer c: clusters.keySet()){
			s += clusters.get(c).size();
		}
		System.out.println("Number of Clusters After additions from extra Trees: " +s);
		System.out.println("Number of STBs After additions from extra Trees: " +STBCountInGeneTrees.size());
	}
	
	private STITreeCluster treeComplementary(STITreeCluster treeAll, STITreeCluster c){
		STITreeCluster res = c.complementaryCluster();
		res.getCluster().and(treeAll.getCluster());
		return res;
	}
	
	int computeTreeSTBipartitions(List<Tree> trees, 
			Map<String, String> taxonMap, Map<Integer, Set<Vertex>> clusters) {

		int sigmaN = 0;
		int k = trees.size();
		String[] leaves = stTaxa;
		int n = leaves.length;
		geneTreeSTBBySize = new ArrayList<Set<STBipartition>>(leaves.length);
		//X = new ArrayList<Set<STBipartition>>(leaves.length);
		//goodSTBs = new ArrayList<Set<STBipartition>>(leaves.length);
		
		for (int i = 0; i <= leaves.length; i++) {
			geneTreeSTBBySize.add(new HashSet<STBipartition>());
			//X.add(new HashSet<STBipartition>());
			//goodSTBs.add(new HashSet<STBipartition>());
			clusters.put(i, new HashSet<Vertex>());
		}
		STBCountInGeneTrees = new HashMap<STBipartition, Integer>(k*n);
		
		X = new HashMap<STITreeCluster, Set<STBipartition>>(k*n);
		
		//gtNodeToSTBs = new HashMap<TNode, STBipartition>(k*n);
		
		STITreeCluster all = new STITreeCluster(stTaxa);
		String as[];
		int j = (as = stTaxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i];
			all.addLeaf(t);
		}
		
		addToClusters(clusters, all, leaves.length);
				
		for (Tree tr : trees) {			
			String[] treeLeaves = tr.getLeaves();			
			STITreeCluster treeAll = new STITreeCluster(leaves);
			for (int i = 0; i < treeLeaves.length; i++) {
				treeAll.addLeaf(treeLeaves[i]);
			}			
			sigmaN += tr.getLeafCount() - 1;
			Map<TNode,STITreeCluster> map = new HashMap<TNode, STITreeCluster>(n);
			for (Iterator<TNode> nodeIt = tr.postTraverse()
					.iterator(); nodeIt.hasNext();) {
				TNode node = nodeIt.next();		        
	            if(node.isLeaf())
	            {
	                String nodeName = node.getName();
	                if (taxonMap != null) {
	                	nodeName = taxonMap.get(nodeName);
	                }

	                STITreeCluster tb = new STITreeCluster(leaves);
	                tb.addLeaf(nodeName);	                	                	        			
	        			
        			addToClusters(clusters, tb, 1);

        			map.put(node, tb);
	                
	                if (!rooted) {
		                addSTB(clusters, tb, treeComplementary(treeAll, tb), null, node);
	                }
	            } else {
	                int childCount = node.getChildCount();
	                STITreeCluster childbslist[] = new STITreeCluster[childCount];
			        BitSet bs = new BitSet(leaves.length);
	                int index = 0;
	                for(Iterator iterator3 = node.getChildren().iterator(); iterator3.hasNext();)
	                {
	                    TNode child = (TNode)iterator3.next();
	                    childbslist[index++] = map.get(child);
	                    bs.or(map.get(child).getCluster());
	                }
	                	                		                
	                STITreeCluster cluster = new STITreeCluster(leaves);
	                cluster.setCluster((BitSet) bs.clone());	
	                
	                int size = cluster.getClusterSize();
	                
        			addToClusters(clusters, cluster, size);
	                map.put(node, cluster);	                
	                
	                
	                if (rooted) {

	                	if (index > 2) {	                	
		                	throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);
		                }

		                STITreeCluster l_cluster = childbslist[0];
		                
		                STITreeCluster r_cluster = childbslist[1];
		                		                
		                addSTB(clusters, l_cluster, r_cluster, cluster, node);
	                } else {		                	
	                	if (childCount == 2) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster r_cluster = childbslist[1];
			                
			                STITreeCluster allMinuslAndr_cluster = treeComplementary(treeAll, cluster);
			                		                
			                STITreeCluster lAndr_cluster = cluster;
			                
			                if (allMinuslAndr_cluster.getClusterSize() != 0) {
				                // add Vertex STBs			                
			                	addSTB(clusters, l_cluster, r_cluster, cluster, node);
			                	addSTB(clusters, r_cluster, allMinuslAndr_cluster, null, node);
			                	addSTB(clusters, l_cluster, allMinuslAndr_cluster, null, node);
			                
			                	// Add the Edge STB
			                	addSTB(clusters, lAndr_cluster, allMinuslAndr_cluster,  all, node);
			                }

	                	} else if (childCount == 3 && node.isRoot()) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster m_cluster =  childbslist[1];
			                
			                STITreeCluster r_cluster = childbslist[2];

			                addSTB(clusters, l_cluster, r_cluster, null, node);
			                addSTB(clusters, r_cluster, m_cluster, null, node);
			                addSTB(clusters, l_cluster, m_cluster, null, node);
	                	} else {
	                		throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);
	                	}
	                }
	            }	           
			}

		}
		System.out.println("STBs in gene trees: " + STBCountInGeneTrees.keySet().size());
		
		int s = 0;
		for (Integer c: clusters.keySet()){
			s += clusters.get(c).size();
		}
		System.out.println("Number of gene tree Clusters: " +s);
		
		weights = new HashMap<STBipartition, Integer>(STBCountInGeneTrees.size());
		
		//System.out.println(average(hash.values()));
		
		return sigmaN;
	}
	
	private double average(Collection<Set<STBipartition>> values) {
		double sum = 0.0;
		for (Set<STBipartition> set : values) {
			sum += set.size();
		}
		return sum/values.size();
	}
	//private HashMap<Integer, Set<STBipartition>> hash = new HashMap<Integer, Set<STBipartition>>();
	private void addSTB(Map<Integer, Set<Vertex>> clusters,
			STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster cluster, TNode node) {
		
		//System.err.println("trying: " + l_cluster + " | " + r_cluster);
		if (cluster == null) {
			cluster = new STITreeCluster(l_cluster);
			cluster.getCluster().or(r_cluster.getCluster());
		}
		int size = cluster.getClusterSize();
		if (l_cluster.isDisjoint(r_cluster)) {
			
			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);		
			geneTreeSTBBySize.get(size).add(stb);
			//gtNodeToSTBs.put(node,stb);
			((STINode)node).setData(stb);
			addSTBToX(clusters, stb);	
			//System.out.println(stb + " hashes to " + stb.hashCode());
			//if (! hash.containsKey(stb.hashCode()))
			//	hash.put(stb.hashCode(), new HashSet<STBipartition>());
			//hash.get(stb.hashCode()).add(stb);			
			STBCountInGeneTrees.put(stb, 
					STBCountInGeneTrees.containsKey(stb)? 
							STBCountInGeneTrees.get(stb)+1: 1);
		} else {
			// This case could happen for multiple-copy
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
			
			if (!rml.getCluster().isEmpty()) {
				addToClusters(clusters, rml, rml.getClusterSize());
				addSTBToX(clusters, new STBipartition(l_cluster, rml, cluster));
			}
			if (!lmr.getCluster().isEmpty()) {
				addToClusters(clusters, lmr, lmr.getClusterSize());
				addSTBToX(clusters, new STBipartition(lmr, r_cluster, cluster));
			}
		}
	}
	
	private void addSTBToX(Map<Integer, Set<Vertex>> clusters,
			STITreeCluster l_cluster, STITreeCluster r_cluster,
			STITreeCluster cluster, TNode node) {
		
		if (cluster == null) {
			cluster = new STITreeCluster(l_cluster);
			cluster.getCluster().or(r_cluster.getCluster());
		}
		if (l_cluster.isDisjoint(r_cluster)) {
			
			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);		
			((STINode)node).setData(stb);
			addSTBToX(clusters, stb);	
		} else {
			// This case could happen for multiple-copy
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
			
			if (!rml.getCluster().isEmpty()) {
				addToClusters(clusters, rml, rml.getClusterSize());
				addSTBToX(clusters, new STBipartition(l_cluster, rml, cluster));
			}
			if (!lmr.getCluster().isEmpty()) {
				addToClusters(clusters, lmr, lmr.getClusterSize());
				addSTBToX(clusters, new STBipartition(lmr, r_cluster, cluster));
			}
		}
	
		
	}


	void addExtraBipartitionsByHeuristics(Map<Integer, Set<Vertex>> clusters) {
		HashSet<STBipartition> bipToAddToX = new HashSet<STBipartition>();
		//goodSTBs = X;
		//if (true) return;
/*		int added = 0;
		for (int i=1; i<goodSTBs.size(); i++)
		{			
			Set<STBipartition> curr_set = goodSTBs.get(i);
			for (STBipartition stb1:curr_set)
			{
				//if (Math.random() < 0.70) continue;
				for (int j=i; j<goodSTBs.size(); j++)
				{
					Set<STBipartition> other_set = goodSTBs.get(j);
					//if (Math.random() < 0.70) continue;
					for (STBipartition stb2:other_set) {
						//System.out.println(stb1 +" **AND** " + stb2);
						if (stb1.cluster1.getClusterSize() < 3 ||
								stb1.cluster2.getClusterSize() < 3 ||
								stb2.cluster1.getClusterSize() < 3 ||
								stb2.cluster2.getClusterSize() < 3) {							
							if (tryToAdd(stb1,stb2,bipToAddToX) != null) added++;						
						}
					}
					System.err.println(bipToAddToX.size() + " " + i);
				}
			}
		}
				
		for (STBipartition stb: bipToAddToX)
		{
			//System.err.println( "Adding: " + stb);
			addSTBToX(clusters, stb);
		}
		System.out.println("\n\nAdded " + added+ " bipartitions:\n");

*/	/*	int s = 0;
		for (Integer c: clusters.keySet()){
			s += clusters.get(c).size();
		}
		System.out.println("Number of Clusters After Addition: " +s);
*/
	}

	
	private STBipartition tryAddingExtraSTB_AndreRule(STBipartition stb1, STBipartition stb2, Set<STBipartition> bipToAddToX)
	{
		if (stb1.equals(stb2)) return null;
		if ( stb1.isDominatedBy(stb2) || stb2.isDominatedBy(stb1) )
			return null;

		if ( stb1.c.isDisjoint(stb2.c) ) return null;

		if ( stb1.cluster1.isDisjoint(stb2.cluster2) && 
				stb1.cluster2.isDisjoint(stb2.cluster1))
		{			
			STITreeCluster cl1 = new STITreeCluster(stb1.cluster1);
			cl1 = cl1.merge(stb2.cluster1);
			STITreeCluster cl2 = new STITreeCluster(stb1.cluster2);
			cl2 = cl2.merge(stb2.cluster2);
			STITreeCluster cl = new STITreeCluster(stb1.c);
			cl = cl.merge(stb2.c);
			STBipartition r = new STBipartition(cl1,cl2,cl);
			bipToAddToX.add(r);
			return r;
		}
		else if ( stb1.cluster1.isDisjoint(stb2.cluster1) && 
				stb1.cluster2.isDisjoint(stb2.cluster2) )
		{
			STITreeCluster cl1 = new STITreeCluster(stb1.cluster1);
			cl1 = cl1.merge(stb2.cluster2);
			STITreeCluster cl2 = new STITreeCluster(stb1.cluster2);
			cl2 = cl2.merge(stb2.cluster1);
			STITreeCluster cl = new STITreeCluster(stb1.c);
			cl = cl.merge(stb2.c);
			STBipartition r = new STBipartition(cl1,cl2,cl);
			bipToAddToX.add(r);
			return r;
		}
		return null;
	}

	private void addSTBToX(Map<Integer, Set<Vertex>> clusters, STBipartition stb) {
		//System.err.println("Adding "+stb);
		int size = stb.c.getClusterSize();
		if (!X.containsKey(stb.c) || !X.get(stb.c).contains(stb)){
			//X.get(size).add(stb);
			STBCountInGeneTrees.put(stb, 0);
			STITreeCluster c = stb.c;
			addToClusters(clusters, c, size);
			Set<STBipartition> stbs = X.get(c);
			stbs = (stbs== null)? new HashSet<STBipartition>() : stbs;
			stbs.add(stb);
			X.put(c, stbs);
		}
	}

	void preCalculateWeights(List<Tree> trees,List<Tree> extraTrees, Map<String, String> taxonMap) {		
		/*weights.putAll(STBCountInGeneTrees);						
		
		for (int i = leaves.length; i > 1; i--) {
			Set<STBipartition> biggerSTBs = X.get(i);			
			for (STBipartition biggerSTB : biggerSTBs) {
				int weight = weights.get(biggerSTB);
				for (int j = i - 1; j > 1; j--) {
					Set<STBipartition> smallerSTBs = geneTreeSTBBySize.get(j);
					for (STBipartition smallerSTB : smallerSTBs) {
						if (smallerSTB.isDominatedBy(biggerSTB)) {
							weight+= STBCountInGeneTrees.get(smallerSTB);
						}
					}
				}
				weights.put(biggerSTB, weight);
			}
		}*/
		if (rooted && taxonMap == null &&
				stTaxa.length > trees.size()) {
			calculateWeightsByLCA(trees,trees);
			if (extraTrees != null) {
				calculateWeightsByLCA(extraTrees,trees);
			}
		}
		
	}

	void calculateWeightsByLCA (List<Tree> stTrees,List<Tree> gtTrees)
	{		
		
		for (Tree stTree:stTrees){
			SchieberVishkinLCA lcaLookup = new SchieberVishkinLCA(stTree);			
			for (Tree gtTree: gtTrees) {
				Stack<TNode> stack = new Stack<TNode>();
				for (TNode gtNode : gtTree.postTraverse()) {
					if (gtNode.isLeaf()) {
						stack.push(stTree.getNode(gtNode.getName()));
					} else {
						TNode rightLCA = stack.pop();
						TNode leftLCA = stack.pop();
						// If gene trees are incomplete, we can have this case
						if (rightLCA == null || leftLCA == null) {
							stack.push(null);
							continue;
						}
						TNode lca = lcaLookup.getLCA(leftLCA, rightLCA);
						stack.push(lca);
						if (lca != leftLCA && lca != rightLCA) {
							// LCA in stTree dominates gtNode in gene tree gtTree
							STBipartition stSTB = (STBipartition) ((STINode)lca).getData();
							STBipartition gtSTB = (STBipartition) ((STINode)gtNode).getData();
							Set<STBipartition> alreadyProcessedSTBs = alreadyProcessed.get(gtSTB);
							
							if (alreadyProcessedSTBs == null) {
								alreadyProcessedSTBs = new HashSet<STBipartition>(gtTrees.size()/4);
								alreadyProcessed.put(gtSTB, alreadyProcessedSTBs);
							}
							
							if (alreadyProcessedSTBs.contains(stSTB)) {
								continue;
							}
							
							weights.put(stSTB, 
									(weights.containsKey(stSTB)?
											weights.get(stSTB):0) 
											+ STBCountInGeneTrees.get(gtSTB));							
							alreadyProcessedSTBs.add(stSTB);
						}
					}									
				}
			}
		}
	}
	
	int calculateMissingWeight(STBipartition biggerSTB) {
		//System.out.println("why am I here? "+biggerSTB);
		Integer count = STBCountInGeneTrees.get(biggerSTB);
		int weight = count != null? count : 0;
		for (int j = biggerSTB.c.getClusterSize() - 1; j > 1; j--) {			
			Set<STBipartition> smallerSTBs = geneTreeSTBBySize.get(j);
			for (STBipartition smallerSTB : smallerSTBs) {
				if (smallerSTB.isDominatedBy(biggerSTB)) {
					weight+= STBCountInGeneTrees.get(smallerSTB);
				}
			}			
		}		
		if (!rooted && biggerSTB.c.getClusterSize() == stTaxa.length) {
			weight -= 2 * count;
			if (biggerSTB.cluster1.getClusterSize() != 1 &&
					biggerSTB.cluster2.getClusterSize() != 1) {
				weight -= 2 * count;				
			}
		}
		return weight;
	}
	
	public DuplicationWeightCounter(String[] stTaxa, boolean rooted){
		this.stTaxa = stTaxa;
		this.rooted = rooted;
	}

	public DuplicationWeightCounter(String[] gtTaxa, String[] stTaxa, boolean rooted2){
		this.gtTaxa = gtTaxa;		
		this.stTaxa = stTaxa;
		this.rooted = rooted2;
	}
	
	public int getBiPartitionDPWeight(STITreeCluster cluster1, STITreeCluster cluster2, STITreeCluster cluster) {		
		STBipartition bi = new STBipartition(cluster1,cluster2, cluster);
		if (!weights.containsKey(bi)){
			weights.put(bi,calculateMissingWeight(bi));
		}
		return weights.get(bi);
	}

/*	public void addGoodSTB (STBipartition good, int size) {		
		goodSTBs.get(size).add(good);
	}
*/	
	public Set<STBipartition> getClusterBiPartitions(STITreeCluster cluster) {
						
		return X.get(cluster);
	}
	
	static public int cnt = 0;
	static class STBipartition {
		
		STITreeCluster cluster1;
		STITreeCluster cluster2;		
		private STITreeCluster c;
		private int _hash = 0;
		
		public STBipartition(STITreeCluster c1, STITreeCluster c2, STITreeCluster cluster) {
			if (c1.getCluster().nextSetBit(0) > c2.getCluster().nextSetBit(0)) {
				cluster1 = c1;
				cluster2 = c2;
			} else {
				cluster1 = c2;
				cluster2 = c1;				
			}
			c = cluster;				
		}
		@Override
		public boolean equals(Object obj) {
			STBipartition stb2 = (STBipartition) obj; 
			
			return this == obj ||
					((stb2.cluster1.equals(this.cluster1) && stb2.cluster2.equals(this.cluster2)));					
		}
		@Override
		public int hashCode() {
			if (_hash == 0) {
				_hash = cluster1.hashCode() * cluster2.hashCode();
			}
			return _hash;
		}
		@Override
		public String toString() {		
			return cluster1.toString()+"|"+cluster2.toString();
		}
		public boolean isDominatedBy(STBipartition dominant) {
			if (! dominant.c.containsCluster(this.c)) {
				return false;
			}
			//cnt++;
			return (dominant.cluster1.containsCluster(this.cluster1) && dominant.cluster2.containsCluster(this.cluster2)) ||
					(dominant.cluster2.containsCluster(this.cluster1) && dominant.cluster1.containsCluster(this.cluster2));
		}
		
	}
	
}