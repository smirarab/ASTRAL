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
	
	public int[] getCardinalities(){
		STITreeCluster[] clusters = getClusters();
		int[] cardialities = new int[clusters.length];
		for (int i = 0; i < cardialities.length; i++){
			cardialities[i] = clusters[i].cardinality();
		}
		return cardialities;
	}
	
	public long selfScore(){
		long sx = 0, sx2 = 0, score = 0;
		int[] cardialities = getCardinalities();
		for (int x: cardialities){
			sx += x;
			sx2 += x * x;
		}
		for (int x: cardialities){
			score += x * (x - 1) * ((sx - x) * (sx - x) - sx2 + x * x);
		}
		return score / 2;
	}
}
