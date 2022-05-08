package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeClusterMP;
import phylonet.util.BitSet;

/**
 * This class is used in the inference class, and does two things:
 *  -- Keeps an instance of the top level IClusterCollection
 *  -- Has methods for building the set X
 * @author smirarab
 *
 * @param <T>
 */
public abstract class AbstractDataCollectionMP<T> extends AbstractDataCollection<T>{


	protected void addAllPossibleSubClusters(STITreeCluster cluster) {
	    int size = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier().taxonCount();
		STITreeCluster c = new STITreeClusterMP((TaxonIdentifierMP)GlobalMaps.taxonNameMap.getSpeciesIdMapper().getSTTaxonIdentifier());
		BitSet bs = c.getBitSet();
		while (true) {
			int tsb = bs.nextClearBit(0);
			if (tsb >= size) {
				break;
			}
			bs.set(tsb);
			bs.clear(0, tsb);
			
			
			c = GlobalMaps.taxonNameMap.getSpeciesIdMapper().getGeneClusterForSTCluster(bs);
			addToClusters(c, c.getClusterSize());
		}
		Logging.log("Number of Clusters After Adding All possible clusters: "
						+ clusters.getClusterCount());
	}

	
	/*@Override
	protected Object clone() throws CloneNotSupportedException {
		return (AbstractDataCollectionMP) super.clone();
	}*/

}
