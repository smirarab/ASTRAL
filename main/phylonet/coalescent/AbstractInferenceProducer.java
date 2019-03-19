package phylonet.coalescent;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractInferenceProducer<T> extends AbstractInference<T> {
	
	private LinkedBlockingQueue<T> queueReadyTripartitions;
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
		IClusterCollection containedVertecies = this.dataCollection.clusters.getContainedClusters(v);
		int clusterSize = v.getCluster().getClusterSize();
		
		Iterable<VertexPair> clusterResolutions;
		
		if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
			clusterResolutions = new ArrayList<VertexPair>();
			Vertex v1 = null;
			int smallestSize = 1;
			while (v1 == null) {
				Set<Vertex> cs = containedVertecies.getSubClusters(smallestSize);
				if (cs.size() != 0) {
					for(Vertex csi : cs) {
						if(csi.getCluster().getBitSet().nextSetBit(0) == 0) {
							v1 = csi;
							break;
						}
					}				
				}
				else 
					smallestSize++;
			}
			for (Vertex v2: containedVertecies.getSubClusters(GlobalMaps.taxonIdentifier.taxonCount()-smallestSize))
			{
				if (v1.getCluster().isDisjoint(v2.getCluster())) {
					VertexPair vp = new VertexPair(v1, v2, v);
					((ArrayList<VertexPair>) clusterResolutions).add(vp);
					break;
				}
				//System.out.println(v2.toString());
			}
			
			try {
				this.getQueueClusterResolutions().put(clusterResolutions);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			
			clusterResolutions = containedVertecies.getClusterResolutions();
			try {
				this.getQueueClusterResolutions().put(clusterResolutions);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	
		}
		return clusterResolutions;
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

	public LinkedBlockingQueue<T> getQueueReadyTripartitions() {
		return queueReadyTripartitions;
	}

	public abstract double scoreGeneTree(Tree scorest, boolean initialize) ;

	abstract IClusterCollection newClusterCollection();
	
	abstract AbstractDataCollection<T> newCounter(IClusterCollection clusters);
	
}
