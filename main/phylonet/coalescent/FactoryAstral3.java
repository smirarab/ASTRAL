package phylonet.coalescent;

import java.util.List;

import phylonet.tree.model.Tree;
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
	public AbstractDataCollection newCounter(IClusterCollection clusters, AbstractInference inference, boolean constrained) {
		if (constrained)
			return new WQDataCollectionConstrained((WQClusterCollection)clusters, inference);
		else
			return new WQDataCollection((WQClusterCollection)clusters, inference);
	}

	@Override
	public GreedyConsensus greedyCons() {
		return GreedyConsensus.instance;
	}

	@Override
	public LoggerInterface newLogger() {
		return new LoggerStdErr();
	}

	@Override
	public AbstractInference newInference(Options inOptions, List<Tree> trees, List<Tree> extraTrees, List<Tree> toRemoveExtraTrees) {
		return  new WQInference(inOptions, trees, extraTrees, toRemoveExtraTrees);
	}


}
