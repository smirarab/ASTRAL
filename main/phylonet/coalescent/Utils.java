package phylonet.coalescent;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

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

    public static Tree buildTreeFromClusters(Iterable<STITreeCluster> clusters) {
        if ((clusters == null) || (!clusters.iterator().hasNext())) {
          System.err.println("Empty list of clusters. The function returns a null tree.");
          return null;
        }
    
        MutableTree tree = new STITree<Double>();
    
        //String[] taxa = ((STITreeCluster)clusters.get(0)).getTaxa();
        for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++) {
          tree.getRoot().createChild(GlobalMaps.taxonIdentifier.getTaxonName(i));
        }
    
        for (STITreeCluster tc : clusters) {
          if ((tc.getClusterSize() <= 1) || (tc.getClusterSize() == GlobalMaps.taxonIdentifier.taxonCount()))
          {
            continue;
          }
    
          Set<TNode> clusterLeaves = new HashSet<TNode>();
          TNode node;
          for (String l : tc.getClusterLeaves()) {
            node = tree.getNode(l);
            clusterLeaves.add(node);
          }
    
          SchieberVishkinLCA lcaFinder = new SchieberVishkinLCA(tree);
          TNode lca = lcaFinder.getLCA(clusterLeaves);
    
          LinkedList<TNode> movedChildren = new LinkedList<TNode>();
          int nodes = clusterLeaves.size();
          for (TNode child : lca.getChildren()) {
            BitSet childCluster = new BitSet(GlobalMaps.taxonIdentifier.taxonCount());
            for (TNode cl : child.getLeaves()) {
              int i = GlobalMaps.taxonIdentifier.taxonId(cl.getName());
              childCluster.set(i);
            }
            
    
            BitSet temp = (BitSet)childCluster.clone();
            temp.and(tc.getBitSet());
            if (temp.equals(childCluster)) {
              movedChildren.add(child);
              nodes -= temp.cardinality();
            }
    
          }
          
          if (movedChildren.size() == 0 || nodes != 0) {
              continue;
          }
    
          STINode newChild = ((STINode)lca).createChild();
    
          while (!movedChildren.isEmpty()) {
            newChild.adoptChild((TMutableNode)movedChildren.get(0));
            movedChildren.remove(0);
          }
        }
    
        ((STITree<Double>)tree).setRooted(false);
        return (Tree)tree;
      }

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

    public static final Tree greedyConsensus(Iterable<Tree> trees) {
    
        HashMap<STITreeCluster, Integer> count = new HashMap<STITreeCluster, Integer>();
        for (Tree tree : trees) {
            List<STITreeCluster> bipartitionClusters = Utils.getClusters(tree);
            for (STITreeCluster cluster: bipartitionClusters) {
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
        TreeSet<Entry<STITreeCluster,Integer>> countSorted = new 
            TreeSet<Entry<STITreeCluster,Integer>>(
                new Comparator<Entry<STITreeCluster,Integer>>() {
    
            @Override
            public int compare(Entry<STITreeCluster, Integer> o1,
                    Entry<STITreeCluster, Integer> o2) {
               int a = o2.getValue().compareTo(o1.getValue());
               if (a != 0) {
                   return a;
               }
               if  (o2.getKey().equals(o1.getKey())) {
                   return 0;
               }
               int i = 0;
               while (i >= 0) {
                   int j = o1.getKey().getBitSet().nextSetBit(i);
                   int jj = o2.getKey().getBitSet().nextSetBit(i);
                   if (j != jj) {
                       return (j > jj) ? -1 : 1;
                   } else {
                       if (j == -1) {
                           break;
                       }
                       i = j + 1;
                   }
               }
               throw new RuntimeException("hmm! this should never be reached");
            }
        });
        countSorted.addAll(count.entrySet());
        List<STITreeCluster> clusters = new ArrayList<STITreeCluster>();       
        for (Entry<STITreeCluster, Integer> entry : countSorted) {
            clusters.add(entry.getKey());
        }
        return Utils.buildTreeFromClusters(clusters);
    }

    
    public static List<STITreeCluster> getClusters(Tree tree){
        List<STITreeCluster> biClusters = new LinkedList<STITreeCluster>();
        Map<TNode, BitSet> map = new HashMap<TNode, BitSet>();
        String[] leaves = GlobalMaps.taxonIdentifier.getAllTaxonNames();
        for (TNode node : tree.postTraverse()) {
            BitSet bs = new BitSet(leaves.length);
            if (node.isLeaf()) {
                // Find the index of this leaf.
                int i = GlobalMaps.taxonIdentifier.taxonId(node.getName());                
                bs.set(i);              
                map.put(node, bs);
            }
            else {
                int childCount = node.getChildCount();
                BitSet[] childbslist = new BitSet[childCount];
                int index = 0;
                for (TNode child : node.getChildren()) {
                    BitSet childCluster = map.get(child);
                    bs.or(childCluster);
                    childbslist[index++] = childCluster;
                }             
                map.put(node, bs);
            }
                          
            if(bs.cardinality()<leaves.length && bs.cardinality()>1){
                STITreeCluster tb = new STITreeCluster();
                tb.setCluster((BitSet)bs.clone());
                if(!biClusters.contains(tb)){
                    biClusters.add(tb);
                }
            }
            
        }
        
        return biClusters;              
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
        } else {
            System.err.println("Command " + args[0]+ " not found.");
        }
    }
}
