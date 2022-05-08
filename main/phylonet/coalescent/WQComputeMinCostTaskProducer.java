package phylonet.coalescent;

import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;

import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.tree.model.sti.STITreeClusterMP;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeClusterMP.VertexMP;

public class WQComputeMinCostTaskProducer extends  AbstractComputeMinCostTaskMP<Tripartition>{


	WQInferenceProducer inference;
	WQDataCollectionMP wqDataCollection;
	//public static final VertexPair  POISON_PILL = new VertexPair(null,null,null);
	
	public WQComputeMinCostTaskProducer(WQInferenceProducer inference, VertexMP v) {
		super(inference, v);
		this.inference = inference;
		this.wqDataCollection = (WQDataCollectionMP)inference.dataCollection;
	}

	//final byte getDoneState = 4;	
	//final byte getProducerDoneState = 3;
	
	
	@Override
	long computeMinCost() throws CannotResolveException {
	
		if ( v.isProdDone() ) {
			return 0;
		}
		//
	
		final int clusterSize = v.getCluster().getClusterSize();
	
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {

			v.setProdDone(true);
	
			return 0;
		}
	
		final IClusterCollection containedVertecies = this.inference.dataCollection.clusters.getContainedClusters(v);
		
		final BlockingQueue<VertexPair> clusterResolutions  = new LinkedBlockingQueue<VertexPair>();
		
		// At the top, there is only one resolution. Find it. 
		if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
			try {
			Vertex v1 = null;
			int smallestSize = -1;
			// If this seem overtly complicated, remember multi-ind cases. 
			while (v1 == null) {
				smallestSize++;
				Set<Vertex> cs = containedVertecies.getSubClusters(smallestSize);
				for(Vertex csi : cs) {
					if(csi.getCluster().getBitSet().nextSetBit(0) == 0) {
						v1 = csi;
						break;
					}
				}
			}
			for (Vertex v2: containedVertecies.getSubClusters(GlobalMaps.taxonIdentifier.taxonCount()-smallestSize))
			{
				if (v1.getCluster().isDisjoint(v2.getCluster())) {
					VertexPair vp = new VertexPair(v1, v2, v);
					(clusterResolutions).put(vp);
					break;
				}
			}
			} catch (Exception e) {
				throw new RuntimeException(e);
			}

			
		} else {

			// Find all cluster resolutions; use the hashing technique. 
			// Add them to the clusterResolution and weight queue at the same time  (locked)
			Future []  futures = new Future[clusterSize / 2];
			final Object lock = new Object();
			for (int j = 1; j <= (clusterSize / 2); j++) {
				final int i = j;
				futures[j-1] = Threading.submit(new Callable<Integer>() {	
					
					public Integer call() {
						Set<Vertex> small = ((HashClusterCollection)containedVertecies).getSmalls(i, clusterSize);
						if (small == null)
							return 0;

						for (Vertex v1b : small) {
							VertexMP v1 = (VertexMP) v1b;
							VertexMP v2 = ((HashClusterCollection)containedVertecies).getCompVertex(v,v1);
							if (v2 == null) {
								continue;
							}
							if (v.getCluster().hash2 == (v1.getCluster().hash2 + v2.getCluster().hash2) ) {
									//&& v.clusterSize == v1.clusterSize + v2.clusterSize // redundant check 
								if (v1.clusterSize != v2.clusterSize || v1.getCluster().hash1 < v2.getCluster().hash1) { // To avoid a pair of the same size twice
									VertexPair bi = new VertexPair(v1, v2, v);
									try {
										// Locking needed to make sure weights and clusters are added in the same order 
										synchronized (lock) {
											clusterResolutions.put(bi);
											inference.getQueueReadyTripartitions().put(STB2T(bi));
										}
									} catch (InterruptedException e) {
										throw new RuntimeException(e);
									}
								}
							}
						}
						return 1; 
					}
					
				});
			}
			// Cluster resolution order matters. Wait for all of them before recursing. 
			for(int i = 0; i < futures.length; i++) {
				try {
					futures[i].get();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

		}

		try {
			
			GlobalQueues.instance.getQueueClusterResolutions().put(clusterResolutions);
			
		} catch (Exception c) {
				throw new RuntimeException(c);
		}

		for (VertexPair bi: clusterResolutions) {
			WQInferenceProducer.weightCount += 1;
			
			newMinCostTask((VertexMP) bi.cluster2).compute();
			newMinCostTask((VertexMP) bi.cluster1).compute();
		}
		
		v.setProdDone(true);

	
		return 0;
	}

	@Override
	public Long getWeight(Tripartition t) {
		throw new RuntimeException("should not call");
	}

	

	@Override
	protected WQComputeMinCostTaskProducer newMinCostTask(VertexMP v) {
		return new WQComputeMinCostTaskProducer((WQInferenceProducer) inference, v);
	}
	
	@Override
	protected Tripartition STB2T(VertexPair vp) {
		return new Tripartition(
				(STITreeClusterMP)vp.cluster1.getCluster(), 
				(STITreeClusterMP)vp.cluster2.getCluster(), 
				(STITreeClusterMP)vp.both.getCluster().complementaryCluster());
	}

	@Override
	protected Long defaultWeightForFullClusters() {
		return 0l;
	}
	
	@Override
	protected long calculateClusterLevelCost() {
		return 0l;
	}


}
