/**
 * 
 */
package phylonet.coalescent;

import java.util.ArrayList;

import phylonet.tree.model.sti.STITreeCluster;

/**
 * @author chaos
 *
 */
public abstract class AbstractPartition {
	public static AbstractPartition createPartition(ArrayList<STITreeCluster> clusters){
		return createPartition(clusters.toArray(new STITreeCluster[]{}));
	}
	
	public static AbstractPartition createPartition(STITreeCluster[] clusters){
		if (clusters.length == 3) return new Tripartition(clusters[0], clusters[1], clusters[2]);
		else return new Polytomy(clusters);
	}
	
	public abstract STITreeCluster[] getClusters();
	
}
