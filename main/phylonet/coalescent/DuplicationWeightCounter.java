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
import java.util.concurrent.RecursiveTask;



import phylonet.coalescent.MGDInference_DP.Vertex;
import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;
import phylonet.coalescent.MGDInference_DP.TaxonNameMap;

public class DuplicationWeightCounter {
	

	HashMap<STBipartition, Integer> weights;
	
	String [] gtTaxa;
	String [] stTaxa;

	private List<Set<STBipartition>> geneTreeSTBBySize;
	
	//private List<Set<STBipartition>> X;
	
	private HashMap<STITreeCluster, Set<STBipartition>> clusterToSTBs;
	
	//private HashMap<TNode, STBipartition> gtNodeToSTBs;

	private Map<STBipartition,Integer> geneTreeSTBCount;
	/**
	 *  Each time a STB is at the root of a gene trees, it will
	 *  be over-counted 4 times if it is of the form a|bcd..
	 *  or 2 times if it is of the form a|b. The following 
	 *  keeps track of gene tree root clusters.
	 */
	private Map<STBipartition,Integer> geneTreeRootSTBs;
	
	//private Map<Tree,List<STBipartition>> treeSTBs;

	//private List<Set<STBipartition>> goodSTBs;

	private boolean rooted;
	
	private HashMap<STBipartition, Set<STBipartition>> alreadyProcessed = new HashMap<STBipartition, Set<STBipartition>>();

	private TaxonNameMap taxonNameMap;

	private Map<Integer, Set<Vertex>> clusters;
	
	private boolean addToClusters (STITreeCluster c, int size) {
		Vertex nv = new Vertex();
		boolean ret = false;
		nv._cluster = c;
		if (!clusters.get(size).contains(nv)){	
			nv._el_num = -1; 			
			clusters.get(size).add(nv);
			ret = true;
			//System.err.println("Clusters: "+clusters);
		}
		
/*		Vertex nv_rev = new Vertex();
		nv_rev._cluster = c.complementaryCluster();
		size = stTaxa.length - size;
		if (!clusters.get(size).contains(nv_rev)){	
			nv_rev._min_cost = -1;
			nv_rev._el_num = -1; 			
			clusters.get(size).add(nv_rev);
			ret = true;
			//System.err.println("Clusters: "+clusters);
		}*/
		return ret;
	}
	

	void addExtraBipartitionsByInput(Map<Integer, Set<Vertex>> clusters,
			List<Tree> trees, boolean extraTreeRooted) {

		int sigmaN = 0;
		int k = trees.size();
		String[] leaves = stTaxa;
		int n = leaves.length;
		
		STITreeCluster all = clusters.get(n).iterator().next()._cluster;
		
		for (Tree tr : trees) {			
			String[] treeLeaves = tr.getLeaves();					
			STITreeCluster treeAll = new STITreeCluster(treeLeaves);				
			
			for (int i = 0; i < treeLeaves.length; i++) {
				String l = treeLeaves[i];
/*				if (taxonMap != null) {
					l = taxonMap.getTaxonName(l);
				}*/
				treeAll.addLeaf(l);
			}	
			sigmaN += tr.getLeafCount() - 1;
			Map<TNode,STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
			Map<TNode,List<String>> map2 = new HashMap<TNode, List<String>>(n);
			for (Iterator<TNode> nodeIt = tr.postTraverse()
					.iterator(); nodeIt.hasNext();) {
				TNode node = nodeIt.next();		        
	            if(node.isLeaf())
	            {
	                String treeName = node.getName();
	                String nodeName = getSpeciesName(treeName);
	                //System.err.println("Adding: "+nodeName);

	                STITreeCluster tb = new STITreeCluster(leaves);
	                tb.addLeaf(nodeName);	               

        			nodeToSTCluster.put(node, tb);
	                
	                if (!extraTreeRooted) {
	                	//TODO: fix the following (first null)
		                tryAddingSTBToX( tb, treeComplementary(null, leaves), null,node);
	                }
	            } else {
	                int childCount = node.getChildCount();
	                STITreeCluster childbslist[] = new STITreeCluster[childCount];
			        BitSet bs = new BitSet(leaves.length);
	                int index = 0;
	                for(Iterator iterator3 = node.getChildren().iterator(); iterator3.hasNext();)
	                {
	                    TNode child = (TNode)iterator3.next();
	                    childbslist[index++] = nodeToSTCluster.get(child);
	                    bs.or(nodeToSTCluster.get(child).getCluster());
	                }
	                	                		                
	                STITreeCluster cluster = new STITreeCluster(leaves);
	                cluster.setCluster((BitSet) bs.clone());	
	                
	                int size = cluster.getClusterSize();
	                
	                nodeToSTCluster.put(node, cluster);	                
	                
	                
	                if (extraTreeRooted) {

	                	if (index > 2) {	                	
		                	throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);
		                }

		                STITreeCluster l_cluster = childbslist[0];
		                
		                STITreeCluster r_cluster = childbslist[1];
		                		                
		                tryAddingSTBToX(l_cluster, r_cluster, cluster,node);
	                } else {		                	
	                	if (childCount == 2) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster r_cluster = childbslist[1];
			                // Fix the following (first null)
			                STITreeCluster allMinuslAndr_cluster = treeComplementary(null, leaves);
			                		                
			                STITreeCluster lAndr_cluster = cluster;
			                
			                // add Vertex STBs
			                tryAddingSTBToX( l_cluster, r_cluster, cluster,node);
			                if (allMinuslAndr_cluster.getClusterSize() != 0) {
			                	tryAddingSTBToX( r_cluster, allMinuslAndr_cluster, null,node);
			                	tryAddingSTBToX( l_cluster, allMinuslAndr_cluster, null,node);
				                tryAddingSTBToX( lAndr_cluster, allMinuslAndr_cluster,  all,node);
			                }			                

	                	} else if (childCount == 3 && node.isRoot()) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster m_cluster =  childbslist[1];
			                
			                STITreeCluster r_cluster = childbslist[2];

			                tryAddingSTBToX( l_cluster, r_cluster, null,node);
			                tryAddingSTBToX( r_cluster, m_cluster, null,node);
			                tryAddingSTBToX( l_cluster, m_cluster, null,node);
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
		System.err.println("Number of Clusters After additions from extra Trees: " +s);
		System.err.println("Number of STBs After additions from extra Trees: " +geneTreeSTBCount);
	}
	
	/*private STITreeCluster treeComplementary(List<String> treeNames, Cluster c , TaxonNameMap taxonMap){
		HashSet<String> set = new HashSet<String> ();
		set.add(cluster);
		return treeComplementary(treeNames, set, taxonMap);
	}*/
	private STITreeCluster treeComplementary(STITreeCluster gtCluster, String[] leaves){
		//System.err.print("Tree complementary of "+gtCluster);
		STITreeCluster newGTCluster = gtCluster.complementaryCluster();		
		//System.err.println(" is: "+newGTCluster.getCluster());
		STITreeCluster newSTCluster = new STITreeCluster(leaves);
		for (String s : newGTCluster.getClusterLeaves()) {
			newSTCluster.addLeaf(getSpeciesName(s));
		}
		//System.err.println("Tree complementary of "+gtCluster+" is: "+newSTCluster);
		return newSTCluster;
	}


	private String getSpeciesName(String geneName) {
		String stName = geneName;
		if (taxonNameMap != null) {
			stName = taxonNameMap.getTaxonName(geneName);
		}
		return stName;
	}
	List<STITreeCluster> treeAlls = new ArrayList< STITreeCluster>();
	
	int computeTreeSTBipartitions(List<Tree> trees, boolean duploss) {

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
		geneTreeSTBCount = new HashMap<STBipartition, Integer>(k*n);
		geneTreeRootSTBs = new HashMap<STBipartition, Integer>(k*n);
		
		clusterToSTBs = new HashMap<STITreeCluster, Set<STBipartition>>(k*n);
		
		//gtNodeToSTBs = new HashMap<TNode, STBipartition>(k*n);
		
		STITreeCluster all = new STITreeCluster(stTaxa);
		String as[];
		int j = (as = stTaxa).length;
		for (int i = 0; i < j; i++) {
			String t = as[i];
			all.addLeaf(t);
		}
		
		addToClusters(all, leaves.length);
				
		for (int t = 0; t < trees.size(); t++ ) {			
			Tree tr = trees.get(t);
			
			STITreeCluster allInducedByGT = new STITreeCluster(stTaxa);
						
			String[] gtLeaves = tr.getLeaves();		
			//STITreeCluster gtAllCluster = new STITreeCluster(gtLeaves);
			for (int i = 0; i < gtLeaves.length; i++) {
				String l = gtLeaves[i];
				//gtAllCluster.addLeaf(l);
				allInducedByGT.addLeaf(getSpeciesName(l));
			}
			treeAlls.add( allInducedByGT);
			int allInducedByGTSize = allInducedByGT.getClusterSize();
			sigmaN += duploss ? 
					(tr.getLeafCount() + 2 * allInducedByGTSize - 3)
					: (tr.getLeafCount() - 1);
			
			//System.err.println((tr.getLeafCount() - 1) + " "+ sigmaN);
			
			Map<TNode,STITreeCluster> nodeToSTCluster = new HashMap<TNode, STITreeCluster>(n);
			Map<TNode,STITreeCluster> nodeToGTCluster = new HashMap<TNode, STITreeCluster>(n);
			for (Iterator<TNode> nodeIt = tr.postTraverse()
					.iterator(); nodeIt.hasNext();) {
				TNode node = nodeIt.next();		
                //System.err.println("Node is:" + node);
	            if(node.isLeaf())
	            {
	                String nodeName = node.getName();

	                STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
	                gtCluster.addLeaf(nodeName);	                	                	        				                
	                nodeName = getSpeciesName(nodeName);	                
	                STITreeCluster cluster = new STITreeCluster(leaves);
	                cluster.addLeaf(nodeName);	                	                	        			

        			addToClusters(cluster, 1);

        			nodeToSTCluster.put(node, cluster);
        			nodeToGTCluster.put(node, gtCluster);
	                
	                if (!rooted) {
		                tryAddingGeneSTBs(cluster, treeComplementary(gtCluster,leaves), null, node, allInducedByGTSize);
	                }
	            } else {
	                int childCount = node.getChildCount();
	                STITreeCluster childbslist[] = new STITreeCluster[childCount];
			        BitSet bs = new BitSet(leaves.length);
			        BitSet gbs = new BitSet(leaves.length);
	                int index = 0;
	                for(Iterator iterator3 = node.getChildren().iterator(); iterator3.hasNext();)
	                {
	                    TNode child = (TNode)iterator3.next();
	                    childbslist[index++] = nodeToSTCluster.get(child);
	                    bs.or(nodeToSTCluster.get(child).getCluster());
	                    gbs.or(nodeToGTCluster.get(child).getCluster());
	                }
	                
	                STITreeCluster gtCluster = new STITreeCluster(gtLeaves);
	                gtCluster.setCluster(gbs);	
	                
	                STITreeCluster cluster = new STITreeCluster(leaves);
	                cluster.setCluster((BitSet) bs.clone());
	                	                
	                int size = cluster.getClusterSize();
	                
        			addToClusters(cluster, size);
	                nodeToSTCluster.put(node, cluster);	                
	                nodeToGTCluster.put(node, gtCluster);
	                

	                if (rooted) {

	                	if (index > 2) {	                	
		                	throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);
		                }

		                STITreeCluster l_cluster = childbslist[0];
		                
		                STITreeCluster r_cluster = childbslist[1];
		                		                
		                tryAddingGeneSTBs(l_cluster, r_cluster, cluster, node, allInducedByGTSize);
	                } else {		                	
	                	if (childCount == 2) {	                		
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster r_cluster = childbslist[1];			                			               
			                
			                STITreeCluster allMinuslAndr_cluster = treeComplementary(gtCluster,leaves);
			                		                
			                STITreeCluster lAndr_cluster = cluster;			                			               
			                
			                if (allMinuslAndr_cluster.getClusterSize() != 0) {
				                // add Vertex STBs			                
			                	tryAddingGeneSTBs(l_cluster, r_cluster, cluster, node, allInducedByGTSize);
			                	tryAddingGeneSTBs( r_cluster, allMinuslAndr_cluster, null, node, allInducedByGTSize);
			                	tryAddingGeneSTBs(l_cluster, allMinuslAndr_cluster, null, node, allInducedByGTSize);
			                
			                	// Add the Edge STB
			                	tryAddingGeneSTBs(lAndr_cluster, allMinuslAndr_cluster,  null, node, allInducedByGTSize);
			                }

	                	} else if (childCount == 3 && node.isRoot()) {
			                STITreeCluster l_cluster = childbslist[0];
			                
			                STITreeCluster m_cluster =  childbslist[1];
			                
			                STITreeCluster r_cluster = childbslist[2];

			                tryAddingGeneSTBs(l_cluster, r_cluster, null, node, allInducedByGTSize);
			                tryAddingGeneSTBs(r_cluster, m_cluster, null, node, allInducedByGTSize);
			                tryAddingGeneSTBs(l_cluster, m_cluster, null, node, allInducedByGTSize);
	                	} else {
	                		throw new RuntimeException("None bifurcating tree: "+
		                			tr+ "\n" + node);
	                	}
	                }
	            }	           
			}

		}
		
		//System.err.println("STBs in gene trees (count): " + geneTreeSTBCount);
		
		int s = 0;
		for (Integer c: geneTreeSTBCount.values()){
			s += c;
		}
		System.err.println("STBs in gene trees (count): " + geneTreeSTBCount.size());		
		System.err.println("STBs in gene trees (sum): " + s);
		
		s = 0;
		for (Integer c: clusters.keySet()){
			s += clusters.get(c).size();
		}
		System.err.println("Number of gene tree Clusters: " +s);
		
		weights = new HashMap<STBipartition, Integer>(geneTreeSTBCount.size());
		
		//System.err.println("sigma n is "+sigmaN);
		
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
	private void tryAddingGeneSTBs(
			STITreeCluster l_cluster,
			STITreeCluster r_cluster, STITreeCluster cluster, TNode node, 
			 int allInducedByGTSize) {
		//System.err.println("before adding: " + STBCountInGeneTrees);
		//System.err.println("Trying: " + l_cluster + "|" + r_cluster);
		if (cluster == null) {
			cluster = new STITreeCluster(l_cluster);
			cluster.getCluster().or(r_cluster.getCluster());
			//System.err.println("Cluster is: "+cluster);
		}
		int size = cluster.getClusterSize();
		if (l_cluster.isDisjoint(r_cluster)) {
			
			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);	
			geneTreeSTBBySize.get(size).add(stb);
			//gtNodeToSTBs.put(node,stb);
			((STINode)node).setData(stb);
			addSTBToX(stb);	
			//System.out.println(stb + " hashes to " + stb.hashCode());
			//if (! hash.containsKey(stb.hashCode()))
			//	hash.put(stb.hashCode(), new HashSet<STBipartition>());
			//hash.get(stb.hashCode()).add(stb);			
			geneTreeSTBCount.put(stb, 
					geneTreeSTBCount.containsKey(stb)? 
							geneTreeSTBCount.get(stb)+1: 1);
			
			if (size == allInducedByGTSize){
				if (! geneTreeRootSTBs.containsKey(stb)) {
					geneTreeRootSTBs.put(stb, 1);
				} else {
					geneTreeRootSTBs.put(stb, geneTreeRootSTBs.get(stb)+1);
				}
			}			
		} else {
			//System.err.println("Adding only to extra");
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
				addToClusters( rml, rml.getClusterSize());
				addSTBToX( new STBipartition(l_cluster, rml, cluster));
			}
			if (!lmr.getCluster().isEmpty()) {
				addToClusters(lmr, lmr.getClusterSize());
				addSTBToX(new STBipartition(lmr, r_cluster, cluster));
			}
		}
	}
	
	private void tryAddingSTBToX(
			STITreeCluster l_cluster, STITreeCluster r_cluster,
			STITreeCluster cluster, TNode node) {
		
		if (cluster == null) {
			cluster = new STITreeCluster(l_cluster);
			cluster.getCluster().or(r_cluster.getCluster());
		}
		if (l_cluster.isDisjoint(r_cluster)) {
			
			STBipartition stb = new STBipartition(l_cluster, r_cluster, cluster);		
			((STINode)node).setData(stb);
			addSTBToX( stb);	
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
				addToClusters( rml, rml.getClusterSize());
				addSTBToX( new STBipartition(l_cluster, rml, cluster));
			}
			if (!lmr.getCluster().isEmpty()) {
				addToClusters( lmr, lmr.getClusterSize());
				addSTBToX( new STBipartition(lmr, r_cluster, cluster));
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

	private void addSTBToX(STBipartition stb) {
		//System.err.println("Adding to X: "+stb+" "+stb.c +" "+clusterToSTBs.containsKey(stb.c));		
		if (clusterToSTBs.containsKey(stb.c) && clusterToSTBs.get(stb.c).contains(stb)){
			return;
		}
		int size = stb.c.getClusterSize();
		// TODO: following line is algorithmically harmless, 
		// but inefficient. is it necessary?
		geneTreeSTBCount.put(stb, 0);
		STITreeCluster c = stb.c;
		addToClusters( c, size);
		Set<STBipartition> stbs = clusterToSTBs.get(c);
		stbs = (stbs== null)? new HashSet<STBipartition>() : stbs;
		stbs.add(stb);
		clusterToSTBs.put(c, stbs);		
		//System.err.println("X updated: "+STBCountInGeneTrees);
		//System.err.println("X updated: "+clusterToSTBs);
	}

	void preCalculateWeights(List<Tree> trees,List<Tree> extraTrees) {		
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
		if (rooted && taxonNameMap == null &&
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
											+ geneTreeSTBCount.get(gtSTB));							
							alreadyProcessedSTBs.add(stSTB);
						}
					}									
				}
			}
		}
	}
	

	public DuplicationWeightCounter(String[] gtTaxa, String[] stTaxa, boolean rooted, TaxonNameMap taxonNameMap, Map<Integer, Set<Vertex>> clusters){
		this.gtTaxa = gtTaxa;		
		this.stTaxa = stTaxa;
		this.rooted = rooted;
		this.taxonNameMap = taxonNameMap;
		this.clusters = clusters;
	}
	
	public Integer getCalculatedBiPartitionDPWeight(STBipartition bi) {		
		if (!weights.containsKey(bi)){
			//weights.put(bi,calculateMissingWeight(bi));
			return null;
		}
		return weights.get(bi);
	}

/*	public void addGoodSTB (STBipartition good, int size) {		
		goodSTBs.get(size).add(good);
	}
*/	
	public Set<STBipartition> getClusterBiPartitions(STITreeCluster cluster) {
						
		return clusterToSTBs.get(cluster);
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
		
		public STBipartition getInducedSTB(STITreeCluster cluster) {
			STITreeCluster lf = new STITreeCluster(c.getTaxa());
			lf.setCluster((BitSet) this.cluster1.getCluster().clone());
			lf.getCluster().and(cluster.getCluster());
			
			STITreeCluster rf = new STITreeCluster(c.getTaxa());
			rf.setCluster((BitSet) this.cluster2.getCluster().clone());
			rf.getCluster().and(cluster.getCluster());

			STITreeCluster cf = new STITreeCluster(c.getTaxa());
			cf.setCluster((BitSet) this.c.getCluster().clone());
			cf.getCluster().and(cluster.getCluster());
			
			return new STBipartition(lf, rf, cf);
		}
		
	}
	class CalculateWeightTask extends RecursiveTask<Integer> {

		/**
		 * 
		 */
		private static final long serialVersionUID = -2614161117603289345L;
		private STBipartition biggerSTB;

		public CalculateWeightTask(STBipartition biggerSTB) {
			this.biggerSTB = biggerSTB;
		}
		
		int calculateMissingWeight() {
			//System.out.println("why am I here? "+biggerSTB);
			//System.err.print("Calculating weight for: " + biggerSTB);
			Integer count = geneTreeSTBCount.get(biggerSTB);
			count = count != null? count : 0;
			int weight = count;
			for (int j = biggerSTB.c.getClusterSize() - 1; j > 1; j--) {			
				Set<STBipartition> smallerSTBs = geneTreeSTBBySize.get(j);
				for (STBipartition smallerSTB : smallerSTBs) {
					if (smallerSTB.isDominatedBy(biggerSTB)) {
						weight+= geneTreeSTBCount.get(smallerSTB);
					}
				}			
			}
			//System.err.print(" ... " + weight);
			if (!rooted) {
				for (STBipartition rootSTB : geneTreeRootSTBs.keySet()) {
					int c = geneTreeRootSTBs.get(rootSTB);
					STBipartition inducedSTB = biggerSTB.getInducedSTB(rootSTB.c);
					if (inducedSTB.equals(rootSTB)){
						weight -= 2 * c;
						//System.err.print(" .. (" + rootSTB +" )" +c+" "+ weight);
						if (inducedSTB.cluster1.getClusterSize() != 1 &&
								inducedSTB.cluster2.getClusterSize() != 1) {
							weight -= 2 * c;			
							//System.err.print(" . " + weight);
						}						
					}
				}
			}
			//System.err.println();
			weights.put(biggerSTB,weight);
			//System.err.println("Weight of " + biggerSTB + " is " + weight);
			return weight;
		}
		
		@Override
		protected Integer compute() {
			return calculateMissingWeight();
		}
		
	}
	public void addAllPossibleSubClusters(STITreeCluster cluster, Map<Integer, HashSet<Vertex>> containedVertecies) {
		int size = cluster.getClusterSize();
		for (int i = cluster.getCluster().nextSetBit(0); i>=0 ;i = cluster.getCluster().nextSetBit(i+1)){
			STITreeCluster c = new STITreeCluster(cluster);
			c.getCluster().clear(i);

			Vertex nv = new Vertex();
			nv._cluster = c;
			if (!containedVertecies.get(size-1).contains(nv)){					
				nv._el_num = -1; 			
				containedVertecies.get(size-1).add(nv);
				//System.err.println("Clusters: "+clusters);
			}

			addAllPossibleSubClusters(c);
		}
	}

	public boolean addAllPossibleSubClusters(STITreeCluster cluster) {
		int size = cluster.getClusterSize();
		boolean ret = false;
		for (int i = cluster.getCluster().nextSetBit(0); i>=0 ;i = cluster.getCluster().nextSetBit(i+1)){
			STITreeCluster c = new STITreeCluster(cluster);
			c.getCluster().clear(i);
			ret |= addToClusters(c, size-1);
			ret |= addAllPossibleSubClusters(c);
		}
		return ret;
	}
	
	public Vertex getCompleteryVertx(Vertex x, STITreeCluster refCluster) {
		STITreeCluster c = x._cluster;		
		Vertex reverse = new Vertex();
		reverse._cluster = new STITreeCluster(refCluster);
		reverse._cluster.getCluster().xor(c.getCluster());
		int size = reverse._cluster.getClusterSize(); 
		if (!clusters.get(size).contains(reverse)){	
			reverse._el_num = -1; 			
			//clusters.get(size).add(reverse);
			return reverse;
			//System.err.println("Clusters: "+clusters);
		}
		return null;
	}

	public boolean addCompleteryVertx(Vertex x, STITreeCluster refCluster) {
		STITreeCluster c = x._cluster;		
		Vertex reverse = new Vertex();
		reverse._cluster = new STITreeCluster(refCluster);
		reverse._cluster.getCluster().xor(c.getCluster());
		int size = reverse._cluster.getClusterSize(); 
		if (!clusters.get(size).contains(reverse)){	
			reverse._el_num = -1; 			
			clusters.get(size).add(reverse);
			return true;
			//System.err.println("Clusters: "+clusters);
		}
		return false;
	}
}
