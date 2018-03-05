package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import phylonet.bits.BitVector;
import phylonet.lca.SchieberVishkinLCA;
import phylonet.tree.io.NewickReader;
import phylonet.tree.io.ParseException;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TMutableNode;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.util.Bipartitions;
import phylonet.util.BitSet;

public class Utils {

	public static ArrayList<TNode> getChildrenAsList(TNode node) {
		ArrayList<TNode> children = new ArrayList<TNode>();
		for (TNode child : node.getChildren()) {
			children.add(child);
		}
		return children;
	}

	/**
	 * For a given set of compatible clusters, it outputs the Tree object
	 * @param clusters
	 * @param identifier
	 * @param keepclusters a boolean indicating whether the data field of each node in the return 
	 *                     tree should be set to the corresponding input cluster
	 * @return
	 */
	public static Tree buildTreeFromClusters(Iterable<STITreeCluster> clusters, 
			TaxonIdentifier identifier, boolean keepclusters) {
		
        if ((clusters == null) || (!clusters.iterator().hasNext())) {
          throw new RuntimeException("Empty list of clusters. The function returns a null tree.");
        }
    
        //TaxonIdentifier spm = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier();
        STITree tree = new STITree<Double>();
        tree.getRoot().setData(identifier.newCluster().complementaryCluster());
    
        // Start from a star tree
        for (int i = 0; i < identifier.taxonCount(); i++) {
          tree.getRoot().createChild(identifier.getTaxonName(i)).setData(identifier.getClusterForNodeName(identifier.getTaxonName(i)));
        }
    
        /**
         *  Add each cluster to the start tree one by one
         */
        for (STITreeCluster tc : clusters) {
          // Let's start with easy cases
          if (tc.getClusterSize() <= 1 || tc.getClusterSize() == identifier.taxonCount())
            continue;
    
          /**  
           * Find the LCA of all nodes inside this cluster
           */
          Set<TNode> clusterLeaves = new HashSet<TNode>();
          for (String l : tc.getClusterLeaves()) {
            TNode node = tree.getNode(l);
            clusterLeaves.add(node);
          }
    
          SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(tree);
          TNode lca = lcaFinder.getLCA(clusterLeaves);
    
          // The set of clusters that will be moved from the LCA
          //    to become the children of this new node that we 
          //    will (potentially) create. 
          LinkedList<TNode> movedChildren = new LinkedList<TNode>();
          int remainingleaves = clusterLeaves.size();
          
          // Go through the children of the LCA
          for (TNode child : lca.getChildren()) {
            STITreeCluster childCluster = (STITreeCluster) ((STINode)child).getData();
            
            // If the child cluster is a substet of this cluster
            //    add it to the set of clusters that will be moved
            //    to become the children of this cluster
            if (tc.containsCluster(childCluster)) {
              movedChildren.add(child);
              remainingleaves -= childCluster.getClusterSize();
            }
    
          }
          
          // This should only happen if our set of clusters are not compatible. 
          // TODO: either update the documentation of the method or remove this check
          if (movedChildren.size() == 0 || remainingleaves != 0) {
              continue;
          }
    
          // Create a new child for this cluster and adopt
          //   the lca children that are a subset of this node
          //   as the new node's children
          STINode newChild = ((STINode)lca).createChild();
          newChild.setData(tc);
          while (!movedChildren.isEmpty()) {
            newChild.adoptChild((TMutableNode)movedChildren.get(0));
            movedChildren.remove(0);
          }
        }
    
        tree.setRooted(false);
        
        if (!keepclusters) {
        	for (TNode node: ((MutableTree)tree).postTraverse()) {
        		((STINode)node).setData(null);
        	}
        }
        return (Tree)tree;
      }

	
	/**
	 * Compute branch support given as set of trees
	 * @param support_tree
	 * @param trees
	 */
    public static final void computeEdgeSupports(STITree<Double> support_tree, Iterable<Tree> trees) {
    
        // generate leaf assignment
        Hashtable<String,Integer> leaf_assignment = new Hashtable<String,Integer>();
        for(TNode n : support_tree.getNodes()) {
            if(n.isLeaf()) {
                leaf_assignment.put(n.getName(), leaf_assignment.size());
            }
        }
    
        // generate all the bipartitions
        Hashtable<BitVector,TNode> support_partitions = new Hashtable<BitVector,TNode>();
        Bipartitions.computeBipartitions(support_tree, leaf_assignment, support_partitions);
    
        LinkedList<Hashtable<BitVector,TNode>> tree_partitions = new LinkedList<Hashtable<BitVector,TNode>>();
        for(Tree t : trees) {
            Hashtable<BitVector,TNode> th = new Hashtable<BitVector,TNode>();
            Bipartitions.computeBipartitions(t, leaf_assignment, th);
            tree_partitions.add(th);
        }
    
        // compute the ratios
        for(Map.Entry<BitVector,TNode> e : support_partitions.entrySet()) {
            BitVector bvcomp = new BitVector(e.getKey());
            bvcomp.not();
    
            int count = 0;
    
            for(Hashtable<BitVector,TNode> h : tree_partitions) {
                if(h.containsKey(e.getKey()) || h.containsKey(bvcomp)) {
                    count++;
                }
            }
            if (!e.getValue().isLeaf())
                ((STINode<Double>) e.getValue()).setData(((double) count) / tree_partitions.size() * 100);
        }
    
        return;
    }

    /**
     * 
     * @param trees
     * @param randomize
     * @param taxonIdentifier
     * @param keepclusters should we keep clusters as node objects
     * @return
     */
    public static final Tree greedyConsensus(Iterable<Tree> trees, boolean randomize,
    		TaxonIdentifier taxonIdentifier, boolean keepclusters) {
    	return greedyConsensus(trees,new double[]{0d}, randomize, 1, taxonIdentifier, keepclusters).iterator().next();
    }
    
    /**
     * Greedy consensus
     * @param trees
     * @param randomize
     * @param taxonIdentifier
     * @param threshold
     * @param geneTreeSkipProb 
     * @return
     */
    /*
    public static final Tree greedyConsensus(Iterable<Tree> trees, boolean randomize,
    		TaxonIdentifier taxonIdentifier, double threshold, double geneTreeKeepProb) {
    	return greedyConsensus(trees,new double[]{threshold}, randomize, 1, taxonIdentifier, geneTreeKeepProb).iterator().next();
	}
    */
    
    /***
     * Greedy consensus with a set of thresholds
     * @param trees
     * @param thresholds
     * @param randomzie
     * @param repeat
     * @param taxonIdentifier
     * @param keepclusters should we keep clusters as node objects
     * @return
     */
    public static final Collection<Tree> greedyConsensus(Iterable<Tree> trees, 
    		double[] thresholds, boolean randomzie, int repeat, 
    		TaxonIdentifier taxonIdentifier, boolean keepclusters) {
    	GlobalMaps.logTimeMessage("Utils 219-222: " + (double)(System.nanoTime()-GlobalMaps.timer)/1000000000);
			
    	List<Tree> outTrees = new ArrayList<Tree>();
        HashMap<STITreeCluster, Integer> count = new HashMap<STITreeCluster, Integer>();
        int treecount = 0;
        for (Tree tree : trees) {
            treecount++;
            List<STITreeCluster> geneClusters = Utils.getGeneClusters(tree, taxonIdentifier); //taxoncount changes
            for (STITreeCluster cluster: geneClusters) {
                if (count.containsKey(cluster)) {
                    count.put(cluster, count.get(cluster) + 1);
                    continue;
                }
            	
            	STITreeCluster comp = cluster.complementaryCluster();
            	if (count.containsKey(comp)) {
                    count.put(comp, count.get(comp) + 1);
                    continue;
                }
                count.put(cluster, 1);
            }
        }        
       
        GlobalMaps.logTimeMessage("Utils 240-243: " + (double)(System.nanoTime()-GlobalMaps.timer)/1000000000);
			
        ArrayList<Future<Tree>> futures = new ArrayList<Future<Tree>>();
        for (int gi = 0; gi < repeat; gi++) {
        	TreeSet<Entry<STITreeCluster,Integer>> countSorted = new 
        			TreeSet<Entry<STITreeCluster,Integer>>(new ClusterComparator(randomzie, taxonIdentifier.taxonCount()));
        
	        countSorted.addAll(count.entrySet());
	        
	        int ti = thresholds.length - 1;
	        double threshold = thresholds[ti];
	        List<STITreeCluster> clusters = new ArrayList<STITreeCluster>();   
	        for (Entry<STITreeCluster, Integer> entry : countSorted) {
	        	if (threshold > (entry.getValue()+.0d)/treecount) {	
	        		List<STITreeCluster> clusterCopy = new ArrayList<STITreeCluster>(clusters);
	        		futures.add(GlobalMaps.eService.submit(new greedyConsensusLoop(taxonIdentifier, keepclusters, clusterCopy)));
	        		ti--;
	        		if (ti < 0) {
	        			break;
	        		}
	        		threshold = thresholds[ti];
	        	}
	    		clusters.add(entry.getKey());
	        }
	        while (ti >= 0) {
        		List<STITreeCluster> clusterCopy = new ArrayList<STITreeCluster>(clusters);
        		futures.add(GlobalMaps.eService.submit(new greedyConsensusLoop(taxonIdentifier, keepclusters, clusterCopy)));
	    		ti--;
	        }
        }
        for(int i = 0; i < futures.size(); i++) {
        	try {
				outTrees.add(0, futures.get(i).get());
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }
        GlobalMaps.logTimeMessage("Utils 269-272: " + (double)(System.nanoTime()-GlobalMaps.timer)/1000000000);
			
        return outTrees;
    }
	
    public static class greedyConsensusLoop implements Callable<Tree>{
		TaxonIdentifier taxonIdentifier;
		boolean keepclusters;
		List<STITreeCluster> clusters;
		public greedyConsensusLoop(TaxonIdentifier taxonIdentifier, boolean keepclusters, List<STITreeCluster> clusters) {
			this.keepclusters = keepclusters;
			this.taxonIdentifier = taxonIdentifier;
			this.clusters = clusters;
		}
		public Tree call() {
			return Utils.buildTreeFromClusters(clusters, taxonIdentifier, keepclusters);
		}
	}
    /**
     * Gives you all clusters in the tree. These are equivalent
     *  of all bipartitions (if you know the set of all leaves in the tree). 
     *  
     * @param tree
     * @param taxonIdentifier
     * @return
     */
    public static List<STITreeCluster> getGeneClusters(Tree tree, 
    		TaxonIdentifier taxonIdentifier ){
        List<STITreeCluster> biClusters = new LinkedList<STITreeCluster>();
        Stack<BitSet> stack = new Stack<BitSet>();
        String[] leaves = taxonIdentifier.getAllTaxonNames();
        int ii=0;
        for (TNode node : tree.postTraverse()) {
            BitSet bs = new BitSet(leaves.length);
            if (node.isLeaf()) {
                // Find the index of this leaf.
//            	if(ii<2)
//            	System.err.println("name"+node.getName());
//            	ii++;
                int i = taxonIdentifier.taxonId(node.getName());                
                bs.set(i);              
                stack.add(bs);
            }
            else {
                int childCount = node.getChildCount();
                BitSet[] childbslist = new BitSet[childCount];
                int index = 0;
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = stack.pop();
                    bs.or(childCluster);
                    childbslist[index++] = childCluster;
                }             
                stack.add(bs);
            }
                          
            if(bs.cardinality()<leaves.length && bs.cardinality()>1){
                STITreeCluster tb = taxonIdentifier.newCluster();
                tb.setCluster((BitSet)bs.clone());
                //if(!biClusters.contains(tb)){
                biClusters.add(tb);
                //}
            }
            
        }
        
        return biClusters;              
    }
    
    public static List<Integer> getRange(int n) {
		List<Integer> range = new ArrayList<Integer>(n);
		for (int j = 0; j < n; j++) {
			range.add(j);
		}
		return range;
	}

	public static List<Integer> getOnes(int n) {
		List<Integer> range = new ArrayList<Integer>(n);
		for (int j = 0; j < n; j++) {
			range.add(1);
		}
		return range;
	}

    public static void main(String[] args) throws IOException{
        if ("--fixsupport".equals(args[0])) {
            String line;
            int l = 0;          
            BufferedReader treeBufferReader = new BufferedReader(new FileReader(args[1]));;
            List<Tree> trees = new ArrayList<Tree>();
            try {
                while ((line = treeBufferReader.readLine()) != null) {
                    l++;
                    if (line.length() > 0) {
                        line = line.replaceAll("\\)[^,);]*", ")");
                        NewickReader nr = new NewickReader(new StringReader(line));

                        Tree tr = nr.readTree();
                        trees.add(tr);
                        String[] leaves = tr.getLeaves();
                        for (int i = 0; i < leaves.length; i++) {
                            GlobalMaps.taxonIdentifier.taxonId(leaves[i]);
                        }
                    }
                }
                treeBufferReader.close();
            } catch (ParseException e) {
                treeBufferReader.close();
                throw new RuntimeException("Failed to Parse Tree number: " + l ,e);
            }
            int k = trees.size();
            System.err.println(k+" trees read from " + args[1]);
            STITree<Double> best = (STITree<Double>) trees.remove(k-1);
            STITree<Double> consensus = (STITree<Double>) trees.remove(k-2);
            
            for (TNode node : best.postTraverse()) {
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            for (TNode node : consensus.postTraverse()) {
                node.setParentDistance(TNode.NO_DISTANCE);
            }
            
            Utils.computeEdgeSupports(consensus, trees);
            Utils.computeEdgeSupports(best, trees);
            trees.add(consensus);
            trees.add(best);
            
            String outfile = args[1] + ".autofixed.tre";
            BufferedWriter outbuffer;
            if (outfile == null) {
                outbuffer = new BufferedWriter(new OutputStreamWriter(System.out));
            } else {
                outbuffer = new BufferedWriter(new FileWriter(outfile));
            }
            for (Tree tree : trees) {
                outbuffer.write(tree.toStringWD()+ " \n");
            }
            outbuffer.flush();
            System.err.println("File "+args[1]+" fixed and saved as " + outfile);
        } else if ("--alltrees".equals(args[0])) {
        	System.out.println(Utils.generateAllBinaryTrees(new String[]{"1","2","3","4","5"}));
        	System.out.println(Utils.generateAllBinaryTreeStrings(new String[]{"1","2","3","4","5"}));
        	System.out.println(Utils.generateAllBinaryTreeStrings(new String[]{"11","12","13","14"}));
        } 
        else {
            System.err.println("Command " + args[0]+ " not found.");
        }
    }
    
	public static String getLeftmostLeaf(TNode from){
		for (TNode node : from.postTraverse()) {
			if (node.isLeaf()) {
				return node.getName();
			}
		}
		throw new RuntimeException("not possible");	
	}
	
	//TODO: change to an iterable
	/**
	 * Given a tree with one sample from each side of a polytomy, 
	 *   this method returns a resolution of the original tree
	 * @param randomSample
	 * @param fullTree
	 * @return The resolution is returned as a set of bitsets
	 */
	public static List<BitSet> getBitsets(HashMap<String,Integer> randomSample, Tree fullTree) {
		
		ArrayList<BitSet> ret = new ArrayList<BitSet>();

		Stack<BitSet> stack = new Stack<BitSet>();
		for (TNode rgtn : fullTree.postTraverse()) {

			if (rgtn.isRoot() && rgtn.getChildCount() == 2) {
				continue;
			}
			BitSet bs = null;
			int legitchildcount = 0;
			if (rgtn.isLeaf()) {
				// Find the index of this leaf.
				if (randomSample.containsKey(rgtn.getName())) {
					bs = new BitSet(randomSample.size());
					int i =  randomSample.get(rgtn.getName());               
					bs.set(i); 
				}
			}
			else {
				int childCount = rgtn.getChildCount();
				for (int i = 0; i < childCount; i++) {
					BitSet pop = stack.pop();
					if (pop != null) {
						if (bs == null) {
							bs = new BitSet(randomSample.size());
						}
						bs.or(pop);
						legitchildcount++;
					}
				}
			}
			stack.push(bs);
			if (bs == null || legitchildcount < 2)
				continue;
			int bsc = bs.cardinality();
			if (bsc < 2 || bsc >= randomSample.size() - 1) {
				continue;
			}       
			ret.add(bs);
		}
		return ret;
	}
	
	public static void randomlyResolve(TNode node) {
			if (node.getChildCount() < 3) {
				return;
			}
			TNode first = node.getChildren().iterator().next();
			List<TNode> children = first.getSiblings();
			children.add(first);
			while (children.size() > 2) {
				TNode c1 = children.remove(GlobalMaps.random.nextInt(children.size()));
				TNode c2 = children.remove(GlobalMaps.random.nextInt(children.size()));
				TMutableNode mnode = (TMutableNode) node;
				TMutableNode newChild = mnode.createChild();
				newChild.adoptChild((TMutableNode) c1);
				newChild.adoptChild((TMutableNode) c2);
				children.add(newChild);
			}
		}

	
	
	public static class ClusterComparator implements Comparator<Entry<STITreeCluster,Integer>> {
		private BSComparator bsComparator;

		public ClusterComparator (boolean randomize, int size) {
			this.bsComparator = new BSComparator(randomize, size);
		}

		@Override
		public int compare(Entry<STITreeCluster, Integer> o1,
				Entry<STITreeCluster, Integer> o2) {
			return this.bsComparator.compare(
					o1.getKey().getBitSet(),o1.getValue(),
					o2.getKey().getBitSet(),o2.getValue());
					
		}
	}

	public static class BSComparator implements Comparator<Entry<BitSet,Integer>> {

		//private boolean random;
		List<Integer> inds;
		public BSComparator (boolean randomize, int size) {
			inds = new ArrayList<Integer>(); 
			for (int i = 0; i < size; i++) {
				inds.add(i);
			}
			if (randomize) {
				Collections.shuffle(inds, GlobalMaps.random);
			}
		}
		@Override
		public int compare(Entry<BitSet, Integer> o1,
				Entry<BitSet, Integer> o2) {
			return compare(o1.getKey(),o1.getValue(), o2.getKey(), o2.getValue());
		}
		private int compare(BitSet k1, Integer v1, BitSet k2, Integer v2) {
			int a = v2.compareTo(v1);
			if (a != 0) {
				return a;
			}
			if  (k1.equals(k2)) {
				return 0;
			}
			for (int ind : this.inds) {
				boolean j = k1.get(ind);
				boolean jj = k2.get(ind);
				if (j != jj) {
					return (j) ? 1 : -1;
				} 
			}
			throw new RuntimeException("hmm! this should never be reached");
		}
	}
	
	public static List<Tree> generateAllBinaryTrees(String[] leaves){
		List<Tree> alltrees = new ArrayList<Tree>();

		int index = 0;
		//build a basic 3 leaf unrooted tree
		STITree threeLeafTree = new STITree(false);
		STINode root = threeLeafTree.getRoot();
		root.createChild(leaves[index++]);
		STINode innode = root.createChild();
		innode.createChild(leaves[index++]);
		innode.createChild(leaves[index++]);
		alltrees.add(threeLeafTree);

		for(;index<leaves.length; index++){
			String leaf = leaves[index];
			List<Tree> temp = new ArrayList<Tree>();
			temp.addAll(alltrees);
			alltrees.clear();
			for(Tree preTree: temp){
				for(TNode n: preTree.postTraverse()){
					if(!n.isLeaf()){
						Iterator it = n.getChildren().iterator();
						STINode lchild = (STINode)(it.next());
						STITree newTree = new STITree(preTree);
						TNode peerChild = newTree.getNode(lchild.getID());
						TNode peerParent = peerChild.getParent();
						STINode newchild = ((STINode<Integer>)peerParent).createChild();
						newchild.adoptChild((TMutableNode)peerChild);
						newchild.createChild(leaf);
						alltrees.add(newTree);
						if(!n.isRoot()){
							STINode rchild = (STINode)(it.next());
							STITree newTree2 = new STITree(preTree);
							TNode peerChild2 = newTree2.getNode(rchild.getID());
							TNode peerParent2 = peerChild2.getParent();
							STINode newchild2 = ((STINode<Integer>)peerParent2).createChild();
							newchild2.adoptChild((TMutableNode)peerChild2);
							newchild2.createChild(leaf);
							alltrees.add(newTree2);
						}
					}
				}
			}
			temp.clear();
		}

		/*List<Tree> temp = new ArrayList<Tree>();
		temp.addAll(alltrees);
		alltrees.clear();
		for(Tree unrootedTree: temp){
			alltrees.addAll(((STITree)unrootedTree).getAllRootingTrees());
		}
		temp.clear();*/
		return alltrees;
	}
	

	public static List<String> generateAllBinaryTreeStrings(String[] leaves){
		List<StringBuffer> alltrees = new ArrayList<StringBuffer>();
		List<String> results = new ArrayList<String>();
		System.err.print("Computing all possible resolutions ");
		int index = 0;
		//build a basic 3 leaf unrooted tree
		StringBuffer threeLeafTree = new StringBuffer();
		threeLeafTree.append("("+leaves[index++]+",("+leaves[index++]+","+leaves[index++]+"));");
		alltrees.add(threeLeafTree);

		for(;index<leaves.length; index++){
			String leaf = leaves[index];
			String insert = ","+leaf+")";
			List<StringBuffer> temp = new ArrayList<StringBuffer>();
			temp.addAll(alltrees);
			alltrees.clear();
			Stack<Integer> stacks = new Stack<Integer>();
			Stack<Integer> stacke = new Stack<Integer>();
			for(StringBuffer preTree: temp){
				Matcher pattern = Pattern.compile(
						"([(])|([)][^,:;)]*)|([;])|(:)|([^,);(:]*)").matcher(preTree);
				while(!pattern.hitEnd()){
					pattern.find();
					String token = pattern.group();
					while ("".equals(token)) {
						pattern.find();
						token = pattern.group();
					}
					if (")".equals(token)) {
						Integer c1s = stacks.pop();
						Integer c1e = stacke.pop();
						Integer c2s = stacks.pop();
						Integer c2e = stacke.pop();
						stacks.push(c2s-1);
						stacke.push(c1e+1);
						StringBuffer newTree1 = new StringBuffer(preTree);
						newTree1.insert(c2e, insert);
						newTree1.insert(c2s, "(");
						if (index == leaves.length - 1) {
							results.add(newTree1.toString());
						} else {
							alltrees.add(newTree1);
						}
						StringBuffer newTree2 = new StringBuffer(preTree);
						newTree2.insert(c1e, insert);
						newTree2.insert(c1s, "(");
						if (index == leaves.length - 1) {
							results.add(newTree2.toString());
						} else {
							alltrees.add(newTree2);
						}
					} else if (";".equals(token)) {
						if (index == leaves.length - 1) {
							results.remove(results.size() - 1);
						} else {
							alltrees.remove(alltrees.size() - 1);
						}
					} else if (!"(".equals(token) && !",".equals(token)) {
						stacks.push(pattern.start());
						stacke.push(pattern.end());
					} 
				}
				/*for (int i = preTree.indexOf(")", 0); i > 0; i = preTree.indexOf(")", i+1)) {
					
				}*/
			}
			temp.clear();
			System.err.print(".");
			System.err.flush();
		}
		System.err.println();

		/*List<Tree> temp = new ArrayList<Tree>();
		temp.addAll(alltrees);
		alltrees.clear();
		for(Tree unrootedTree: temp){
			alltrees.addAll(((STITree)unrootedTree).getAllRootingTrees());
		}
		temp.clear();*/
		return results;
	}
	
	public static Integer increment(Integer i) {
		return i == null? 1 : (i+1);
	}
}
