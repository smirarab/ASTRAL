package phylonet.coalescent;

import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeClusterMP.VertexMP;

public class WQInferenceProducer extends AbstractInference<Tripartition> {
	
	private BlockingQueue<Tripartition> queueReadyTripartitions;
	public static int weightCount = 0;
	
	
	public WQInferenceProducer(WQInferenceConsumerMP in) {
		super(in.options, in.trees, in.extraTrees, in.toRemoveExtraTrees);
		this.dataCollection =  (AbstractDataCollection<Tripartition>) in.dataCollection;
		GlobalQueues.instance.setQueueClusterResolutions(GlobalQueues.instance.getQueueClusterResolutions());
		this.queueReadyTripartitions = (new LinkedBlockingQueue<Tripartition>());
	}



	public List<Solution> inferSpeciesTree() {

		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		return super.inferSpeciesTree();
	}


	@Override
	List<Solution> processSolutions(Vertex all) {
		// No need to process any solution on the producer side
		// This is an important overwrite.
		return null;
	}

	@Override
	public int countWeights() {
		return weightCount;
	}

	/**
	 * Sets up data structures before starting DP
	 */
	public void setup() {
		this.weightCalculator = null;
		this.setupMisc();
	}

	public BlockingQueue getQueueReadyTripartitions() {
		return queueReadyTripartitions;
	}



	@Override
	public Long getTotalCost(Vertex all) {
		return 0l;
	}
	
	@Override
	public AbstractComputeMinCostTask<Tripartition> newComputeMinCostTask(
			AbstractInference<Tripartition> inference, Vertex all, IClusterCollection clusters) {
		return new WQComputeMinCostTaskProducer((WQInferenceProducer) inference, (VertexMP) all);
	}

	public IClusterCollection newClusterCollection() {
		return new HashClusterCollection(GlobalMaps.taxonIdentifier.taxonCount());
	}

	public AbstractDataCollection<Tripartition> newCounter(IClusterCollection clusters) {
		return new WQDataCollectionMP((HashClusterCollection)clusters, this);
	}


	@Override
	public double scoreSpeciesTreeWithGTLabels(Tree scorest, boolean initialize) {
		throw new RuntimeException("Not Implemented");
	}

	@Override
	protected void initializeWeightCalculator() {
		throw new RuntimeException("Not Implemented");

	}
	
	@Override
	public void setupMisc() {
	}


	public AbstractWeightCalculator<Tripartition> newWeightCalculator() {
		throw new RuntimeException("Not Implemented");
	}



}
