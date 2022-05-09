package phylonet.dl;

import java.util.List;

import phylonet.coalescent.AbstractComputeMinCostTask;
import phylonet.coalescent.GlobalMaps;
import phylonet.coalescent.IClusterCollection;
import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;

public class DLComputeMinCostTask extends AbstractComputeMinCostTask<STBipartition>{

	DLInference inference;
	DLDataCollection dataCollection;
	DLWeightCalculator weightCalculator;
	
	public DLComputeMinCostTask(DLInference inference, Vertex v,
			IClusterCollection clusters) {
		super(inference, (VertexASTRAL3) v, null);
		this.inference = inference;
		dataCollection = (DLDataCollection)inference.dataCollection;
		weightCalculator = (DLWeightCalculator) inference.weightCalculator;
	}
	
	protected double adjustWeight(long clusterLevelCost, Vertex smallV,
			Vertex bigv, Long Wdom) {
		double c;
		if (inference.getOptimizeDuploss() == 3) {
			
			int OverlappingGeneCount = 0;
			int someSideMissingXLCount = 0;
			int bothSidesPresentGeneCount = 0;
			for (int k = 0; k < inference.trees.size(); k++) {
				STITreeCluster treeAll = dataCollection.treeAlls.get(k);
				Tree tree = inference.trees.get(k);
				boolean pDisJoint = smallV.getCluster().isDisjoint(treeAll);
				boolean qDisJoint = bigv.getCluster().isDisjoint(treeAll);
				if (pDisJoint || qDisJoint) {
					someSideMissingXLCount +=  GlobalMaps.taxonNameMap == null ?
						DeepCoalescencesCounter.getClusterCoalNum_rooted(tree, this.v.getCluster()):
						DeepCoalescencesCounter.getClusterCoalNum_rootedMap(tree, this.v.getCluster());
				}
				if (!pDisJoint && !qDisJoint) {
					bothSidesPresentGeneCount += 1;
				}
				if (!pDisJoint || !qDisJoint) {
					OverlappingGeneCount += 1;
				}							
			}
			
			c = - ( clusterLevelCost - 3 * Wdom - OverlappingGeneCount - someSideMissingXLCount 
					+ inference.getDLbdWeigth() * (someSideMissingXLCount + OverlappingGeneCount + bothSidesPresentGeneCount) );
		} else {
			c = Wdom;
		}
		return c;
	}
	
	@Override
	protected long scoreBaseCase(boolean rooted, List<Tree> trees) {
		long _el_num = -1;
		if (inference.getOptimizeDuploss() == 3) {
			if (GlobalMaps.taxonNameMap == null) {
				_el_num = DeepCoalescencesCounter.getClusterCoalNum(trees,
						v.getCluster(), rooted);
				// System.out.println(v + " XL is " + _el_num);
			} else {
				_el_num = DeepCoalescencesCounter.getClusterCoalNumMap(trees,
						v.getCluster(), rooted);
			}
		} else {
			_el_num = 0;
		}
		return - _el_num;
	}

	@Override
	protected AbstractComputeMinCostTask<STBipartition> newMinCostTask(
			 Vertex v, IClusterCollection clusters) {
		return new DLComputeMinCostTask(inference, v, clusters);
	}
	
	@Override
	protected long calculateClusterLevelCost() {
		if (inference.getOptimizeDuploss() == 3) {
			return weightCalculator.calculateDLstdClusterCost(
					this.v.getCluster(), inference.trees);
		}
		return 0;
	}


	@Override
	protected STBipartition STB2T(VertexPair vp) {
		return new STBipartition(vp.cluster1.getCluster(), vp.cluster2.getCluster(), vp.both.getCluster());
	}
	
	int calculateDLbdAdjustment(Vertex smallV, Vertex bigv) {
		int e = 0;
		for (int k = 0; k < inference.trees.size(); k++) {
			STITreeCluster treeAll = dataCollection.treeAlls.get(k);
			boolean pDisJoint = smallV.getCluster().isDisjoint(treeAll);
			boolean qDisJoint = bigv.getCluster().isDisjoint(treeAll);
			if (!pDisJoint && !qDisJoint) {
				e += 1;
			}
			// System.err.println("E for " + v.getCluster() + " is "+e +
			// " and k is  " + k);
			// System.out.println(bigv + "|" + smallV+ " Total is " + e +
			// " extra is "+extraTerms + " for tree" +k);
		}
		return e;
	}

	@Override
	protected Long defaultWeightForFullClusters() {
		return null;
	}

}
