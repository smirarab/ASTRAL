package phylonet.coalescent;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractInferenceProducer<T> extends AbstractInference<T> {
	
	private BlockingQueue<T> queueReadyTripartitions;
	public int weightCount = 0;
	
	public AbstractInferenceProducer(Options options, List<Tree> trees,
			List<Tree> extraTrees, List<Tree> toRemoveExtraTrees) {
		super(options, trees, extraTrees, toRemoveExtraTrees);
	}
	
	public AbstractInferenceProducer(AbstractInference in) {
		super(in.options, in.trees, in.extraTrees, in.toRemoveExtraTrees);
		this.dataCollection = in.dataCollection;
		this.setQueueClusterResolutions(in.getQueueClusterResolutions());
		this.queueReadyTripartitions = (new LinkedBlockingQueue<T>());
	}


	public List<Solution> inferSpeciesTree() {

		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		return super.inferSpeciesTree();
	}

	@Override
	public Iterable<VertexPair> getClusterResolutions(Vertex v) {
		throw new RuntimeException("Shouldn't be called");
	}

	@Override
	public int countWeights() {
		return this.weightCount;
	}

	/**
	 * Sets up data structures before starting DP
	 */
	void setup() {
		this.weightCalculator = null;
		this.setupMisc();
	}

	public BlockingQueue<T> getQueueReadyTripartitions() {
		return queueReadyTripartitions;
	}

	public abstract double scoreGeneTree(Tree scorest, boolean initialize) ;

	abstract IClusterCollection newClusterCollection();
	
	abstract AbstractDataCollection<T> newCounter(IClusterCollection clusters);
	
}
