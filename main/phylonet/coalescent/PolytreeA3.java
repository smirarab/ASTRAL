package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;

public class PolytreeA3 extends Polytree{	


	/**
	 * A node in the input gene trees. Used temporarily to build the polytree. 
	 * @author smirarab
	 *
	 */
	final class PTNode extends Polytree.PTNode{

		/**
		 * To be used for leaves
		 * @param n
		 */
		PTNode(TNode n){
			super();
			STITreeCluster c = Factory.instance.newCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(GlobalMaps.taxonIdentifier.taxonId(n.getName()));
			cluster = PolytreeA3.this.clusters.get(c);
			children = new ArrayList();
		}
		PTNode(ArrayList ch, STITreeCluster s){
			super();
			children = ch;
			STITreeCluster c = Factory.instance.newCluster(GlobalMaps.taxonIdentifier);
			ArrayList<STITreeCluster> cs = new ArrayList<STITreeCluster>();
			for (phylonet.coalescent.Polytree.PTNode child: children){
				child.parent = this;
				c.getBitSet().xor(child.cluster.clusterRef.getBitSet());
				cs.add(child.cluster.clusterRef);
			}
			cluster = findCluster(c, this);
			if (c.equals(s) == false){
				STITreeCluster xc = Factory.instance.newCluster(GlobalMaps.taxonIdentifier);
				xc.getBitSet().xor(c.getBitSet());
				xc.getBitSet().xor(s.getBitSet());
				cs.add(findCluster(xc, null).clusterRef);
			}
			if (cs.size() >= 3){
				AbstractPartition p = AbstractPartition.createPartition(cs);
				if (PolytreeA3.this.partitions.containsKey(p)){
					partition = PolytreeA3.this.partitions.get(p);
					partition.cardinality++;
				}
				else partition = new PTPartition(p, this);
			}
		}
		PTCluster findCluster(STITreeCluster c, PTNode n){
			if (PolytreeA3.this.clusters.containsKey(c)) {
				PTCluster cluster = PolytreeA3.this.clusters.get(c);
				if (cluster.firstNode == null) cluster.firstNode = n;
				return cluster;
			}
			else return new PTCluster(c, n);
		}

		void addAllPartitions(){
			for (phylonet.coalescent.Polytree.PTNode child: children){
				((PTNode) child).addAllPartitions();
			}
			if (isFirstResolutionOfCluster()){
				isUsed = true;
				for (phylonet.coalescent.Polytree.PTNode child: children){
					((PTNode) child).setClusterFlagsByClusterPrecedence();
				}
			}
		}

	}

	int listSize = 0;
	long[] sx = new long[3], sxy = new long[3];
	int[] treeTotal = new int[3];

	public PolytreeA3(List<Tree> trees, WQDataCollection dataCollection){
		super(trees, dataCollection);
	}

	@Override
	protected void initalizePolyTree(List<Tree> trees, WQDataCollection dataCollection) {
		long t = System.currentTimeMillis();

		// Create singleton clusters and add to map
		for (int i = 0; i < GlobalMaps.taxonIdentifier.taxonCount(); i++){
			STITreeCluster c = Factory.instance.newCluster(GlobalMaps.taxonIdentifier);
			c.getBitSet().set(i);
			new PTCluster(c);
		}

		// Represent gene trees as PTNodes
		Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
		for (Tree tr: trees){
			nodeRoots.add(buildTree(tr.getRoot(), tit.next()));
		}

		// Set the isUsedFlag on gene tree nodes
		for (Polytree.PTNode nn: nodeRoots){
			PTNode n = (PTNode) nn;
			n.addAllPartitions();
		}

		// For each node of each gene tree, decide how it
		//  should be treated when calculating weights.
		// Options (non-exclusive) are:
		//   - compute the intersection for the cluster or retrieve it from the list
		//   - compute the weight for the partition if it's the first resolution
		for (Polytree.PTNode nn: nodeRoots){
			PTNode n = (PTNode) nn;
			queueBuilder.add(-1);
			n.buildInstructionQueue();
		}

		queue = mapToInt(queueBuilder);
		clusters = null;
		partitions = null;
		queueBuilder = null;

		STITreeCluster c = (Factory.instance.newCluster(GlobalMaps.taxonIdentifier)).complementaryCluster();
		maxScore = WQWeightByTraversal(new Tripartition(c, c, c, false));
		Logging.log("Polytree max score: " + maxScore / 4);
		Logging.log("Polytree building time: " + (System.currentTimeMillis() - t) / 1000.0D + " seconds.");
	}

	@Override
	protected PTNode newPTNode(TNode node) {
		return new PTNode(node);
	}

	@Override
	protected PTNode newPTNode(ArrayList cs, STITreeCluster s) {
		return new PTNode(cs,s);
	}


	public static long getTime() {
		return time;
	}

	public static void setTime(long time) {
		PolytreeA3.time = time;
	}
}
