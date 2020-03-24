package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class WQInferenceProducer extends AbstractInference<Tripartition> {
	
	private BlockingQueue<Tripartition> queueReadyTripartitions;
	public static int weightCount = 0;
	
	public WQInferenceProducer(Options options, List<Tree> trees,
			List<Tree> extraTrees, List<Tree> toRemoveExtraTrees) {
		super(options, trees, extraTrees, toRemoveExtraTrees);
	}
	
	public WQInferenceProducer(AbstractInference in) {
		super(in.options, in.trees, in.extraTrees, in.toRemoveExtraTrees);
		this.dataCollection = in.dataCollection;
		this.setQueueClusterResolutions(in.getQueueClusterResolutions());
		this.queueReadyTripartitions = (new LinkedBlockingQueue<Tripartition>());
	}



	public List<Solution> inferSpeciesTree() {

		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		return super.inferSpeciesTree();
	}


	@Override
	public int countWeights() {
		return weightCount;
	}

	/**
	 * Sets up data structures before starting DP
	 */
	void setup() {
		this.weightCalculator = null;
		this.setupMisc();
	}

	public BlockingQueue getQueueReadyTripartitions() {
		return queueReadyTripartitions;
	}



	@Override
	Long getTotalCost(Vertex all) {
		return 0l;
	}

	@Override
	AbstractComputeMinCostTask newComputeMinCostTask(
			AbstractInference<Tripartition> inference, Vertex all) {
		return new WQComputeMinCostTaskProducer((WQInferenceProducer) inference, all);
	}

	IClusterCollection newClusterCollection() {
		return new HashClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}

	AbstractDataCollection newCounter(IClusterCollection clusters) {
		return new WQDataCollection((HashClusterCollection)clusters, this);
	}


	@Override
	public double scoreSpeciesTreeWithGTLabels(Tree scorest, boolean initialize) {
		throw new RuntimeException("Not Implemented");
	}

	@Override
	void initializeWeightCalculator() {
		throw new RuntimeException("Not Implemented");

	}
	
	@Override
	void setupMisc() {
	}


	AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		throw new RuntimeException("Not Implemented");
	}

}
