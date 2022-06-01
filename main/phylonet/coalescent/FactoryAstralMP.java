package phylonet.coalescent;

import java.util.List;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeClusterMP;
import phylonet.util.BitSet;

public class FactoryAstralMP extends Factory {

	@Override
	public STITreeCluster newCluster(STITreeCluster c1) {
		if (c1 instanceof HashOnlyTreeCluster)
			return new HashOnlyTreeCluster((STITreeClusterMP)c1);
		else 
			return new STITreeClusterMP((STITreeClusterMP) c1);
	}

	@Override
	public STITreeCluster newCluster(TaxonIdentifier taxonIdentifier) {
		return new STITreeClusterMP(taxonIdentifier);
	}

	@Override
	public STITreeCluster newCluster(TaxonIdentifier taxid, BitSet clone) {
		return new STITreeClusterMP(taxid, clone);
	}


	@Override
	public SpeciesMapper newSpeciesMapper() {
		return new SpeciesMapper(GlobalMaps.taxonIdentifier.taxonCount());
	}

	@Override
	public TaxonIdentifier newTaxonIdentifier() {
		return new TaxonIdentifierMP();
	}

	@Override
	public SimilarityMatrix newSimilarityMatrix(int n) {
		return new SimilarityMatrixMP(n);
	}

	@Override
	public SimilarityMatrix newSimilarityMatrix(float[][] from) {
		return new SimilarityMatrixMP(from);
	}

	@Override
	public WQDataCollectionMP newCounter(IClusterCollection clusters, AbstractInference inference, boolean constrained) {
		if (constrained)
			return new WQDataCollectionConstrainedMP((HashClusterCollection)clusters, inference);
		else 
			return new WQDataCollectionMP((HashClusterCollection)clusters, inference);
	}

	@Override
	public GreedyConsensus greedyCons() {
		return GreedyConsensusMP.instance;
	}

	@Override
	public LoggerInterface newLogger() {
		return new LoggerThreaded();
	}

	@Override
	public AbstractInference newInference(Options inOptions, List<Tree> trees, List<Tree> extraTrees,
			List<Tree> toRemoveExtraTrees) {
		return new WQInferenceConsumerMP(inOptions, trees, extraTrees, toRemoveExtraTrees);
	}
}
