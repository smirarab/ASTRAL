package phylonet.coalescent;

import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;

/**
 * This class is used in the inference class, and does two things:
 *  -- Keeps an instance of the top level IClusterCollection
 *  -- Has methods for building the set X
 * @author smirarab
 *
 * @param <T>
 */
public abstract class AbstractDataCollection <T> {

	protected IClusterCollection clusters;

	protected boolean addToClusters(STITreeCluster c, int size) {
		
		//String bug = "111_0_4, 111_0_0, 81_0_2, 102_0_0, 102_0_2, 102_0_3, 102_0_1, 102_0_4, 111_0_3, 111_0_2, 111_0_1, 53_0_2, 53_0_3, 53_0_4, 44_0_0, 53_0_1, 119_0_4, 119_0_0, 119_0_2, 119_0_1, 119_0_3, 44_0_3, 53_0_0, 44_0_4, 44_0_1, 44_0_2, 81_0_3, 81_0_0, 81_0_4, 81_0_1";
		//String bug = "111_0_4, 111_0_0, 116_0_1, 81_0_2, 102_0_0, 102_0_2, 102_0_3, 102_0_1, 102_0_4, 105_0_4, 105_0_1, 111_0_3, 111_0_2, 111_0_1, 116_0_2, 116_0_0, 116_0_4, 116_0_3, 53_0_2, 53_0_3, 53_0_4, 44_0_0, 53_0_1, 105_0_0, 105_0_2, 105_0_3, 119_0_4, 119_0_0, 119_0_2, 119_0_1, 119_0_3, 44_0_3, 53_0_0, 44_0_4, 44_0_1, 44_0_2, 81_0_3, 81_0_0, 81_0_4, 81_0_1";
		//String[] bugs = bug.split(", ");
		//Set<String> clusterSet = new HashSet<String>(Arrays.asList(c.getClusterLeaves()));
		//Set<String> bugSet = new HashSet<String>(Arrays.asList(bugs));	
		/*if(bugSet.equals(clusterSet)){
            System.err.println(c.toString());
            throw new RuntimeException();
		}*/
	
		if (size == 0) return false;
		Vertex nv = c.new Vertex(size);
		return clusters.addCluster(nv, size);
	}
	
	protected boolean removeCluster(STITreeCluster c, int size) {
		if (size == GlobalMaps.taxonIdentifier.taxonCount()
		|| size == 0) {
			return false;
		}	
		Vertex nv = c.new Vertex(size);
		return clusters.removeCluster(nv, size);
	}

	void addAllPossibleSubClusters(STITreeCluster cluster) {
	    int size = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier().taxonCount();
		STITreeCluster c = new STITreeCluster(GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier());
		BitSet bs = c.getBitSet();
		while (true) {
			int tsb = bs.nextClearBit(0);
			if (tsb >= size) {
				break;
			}
			bs.set(tsb);
			bs.clear(0, tsb);
			
			
			c = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getGeneClusterForSTCluster(bs);;
			//c.setCluster((BitSet) bs.clone());
			
			addToClusters(c, c.getClusterSize());
		}
		Logging.log("Number of Clusters After Adding All possible clusters: "
						+ clusters.getClusterCount());
	}

	/*	
	void addAllPossibleSubClusters2(STITreeCluster cluster) {
		int size = cluster.getClusterSize();
		for (int i = cluster.getBitSet().nextSetBit(0); i >= 0; i = cluster
				.getBitSet().nextSetBit(i + 1)) {
			STITreeCluster c = new STITreeCluster(cluster);
			c.getBitSet().clear(i);
	
			Vertex nv = c.new Vertex();
			this.clusters.addCluster(nv, size - 1);
	
			addAllPossibleSubClusters2(c);
		}
	}*/
	
	
	public abstract void addExtraBipartitionsByInput(
			List<Tree> trees, boolean extraTreeRooted);
	
	public abstract void removeExtraBipartitionsByInput(List<Tree> extraTrees,
			boolean extraTreeRooted);
	
	
	public abstract void formSetX(AbstractInference<T> inference);
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		return (AbstractDataCollection) super.clone();
	}

}
