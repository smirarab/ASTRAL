package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public abstract class Factory {
	public static Factory instance;
	
	public abstract STITreeCluster newCluster(STITreeCluster c1);
	
	public abstract STITreeCluster newCluster(TaxonIdentifier taxonIdentifier);
	
	public abstract STITreeCluster newCluster(TaxonIdentifier taxid, BitSet clone);

	public abstract SpeciesMapper newSpeciesMapper();

	public abstract TaxonIdentifier newTaxonIdentifier();
	
	public abstract SimilarityMatrix newSimilarityMatrix(int n);

	public abstract SimilarityMatrix newSimilarityMatrix(float[][] from);

	public abstract AbstractDataCollection newCounter(IClusterCollection clusters, AbstractInference inference);

	public abstract GreedyConsensus greedyCons();

}
