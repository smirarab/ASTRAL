package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.util.BitSet;

public class FactoryAstral3 extends Factory {

	@Override
	public STITreeCluster newCluster(STITreeCluster c1) {
			return new STITreeCluster(c1);
	}

	@Override
	public STITreeCluster newCluster(TaxonIdentifier taxonIdentifier) {
		return new STITreeCluster(taxonIdentifier);
	}

	@Override
	public STITreeCluster newCluster(TaxonIdentifier taxid, BitSet clone) {
		return new STITreeCluster(taxid, clone);
	}
	
	@Override
	public SpeciesMapper newSpeciesMapper() {
		return new SpeciesMapper(GlobalMaps.taxonIdentifier.taxonCount());
	}

	@Override
	public TaxonIdentifier newTaxonIdentifier() {
		return new TaxonIdentifier();
	}

	@Override
	public SimilarityMatrix newSimilarityMatrix(int n) {
		return new SimilarityMatrix(n);
	}

	@Override
	public SimilarityMatrix newSimilarityMatrix(float[][] from) {
		return new SimilarityMatrix(from);
	}

	@Override
	public AbstractDataCollection newCounter(IClusterCollection clusters, AbstractInference inference) {
		return new WQDataCollection((WQClusterCollection)clusters, inference);
	}


}
