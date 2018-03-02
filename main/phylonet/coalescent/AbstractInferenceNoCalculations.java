package phylonet.coalescent;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.MutableTree;
import phylonet.tree.model.TNode;
import phylonet.tree.model.Tree;
import phylonet.tree.model.sti.STINode;
import phylonet.tree.model.sti.STITree;
import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.util.Collapse;

public abstract class AbstractInferenceNoCalculations<T> extends AbstractInference<T> {
	
	LinkedBlockingQueue<Tripartition> queue1;
	LinkedBlockingQueue<Iterable<VertexPair>> queue4;
	public AbstractInferenceNoCalculations(Options options, List<Tree> trees,
			List<Tree> extraTrees) {
		super(options, trees, extraTrees);
		initDF();
	}
	
	public AbstractInferenceNoCalculations(AbstractInference in) {
		super(in.options, in.trees, in.extraTrees);
		this.dataCollection = in.dataCollection;
		initDF();
	}

	private void initDF() {
		df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		DecimalFormatSymbols dfs = DecimalFormatSymbols.getInstance();
		dfs.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(dfs);
	}

	public abstract double scoreGeneTree(Tree scorest, boolean initialize) ;
	

	public List<Solution> inferSpeciesTree() {

		if (! this.options.isRunSearch() ) {
			System.exit(0);
		}
		return super.inferSpeciesTree();
	}

	abstract IClusterCollection newClusterCollection();
	
	abstract AbstractDataCollection<T> newCounter(IClusterCollection clusters);
	
	abstract AbstractWeightCalculatorNoCalculations<T> newWeightCalculator();

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
				this.queue4.put(clusterResolutions);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			
			clusterResolutions = containedVertecies.getClusterResolutions();
			try {
				this.queue4.put(clusterResolutions);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
		return clusterResolutions;
	}
	
	/**
	 * Sets up data structures before starting DP
	 */
	void setup() {
		weightCalculator = newWeightCalculator();
		this.initializeWeightCalculator();
		this.setupMisc();
	}
}
