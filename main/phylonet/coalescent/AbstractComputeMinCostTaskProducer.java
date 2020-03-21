package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.coalescent.IClusterCollection.VertexPair;
import java.util.*;
import java.util.concurrent.*; 

public abstract class AbstractComputeMinCostTaskProducer<T> extends  AbstractComputeMinCostTask<T>{


	AbstractInferenceProducer<T> inference;
	//public static final VertexPair  POISON_PILL = new VertexPair(null,null,null);
	
	public AbstractComputeMinCostTaskProducer(AbstractInferenceProducer<T> inference, Vertex v) {
		super(inference, v);
		this.inference = inference;
	}

	byte getDoneState() {
		return 3;
	}
	
	byte getOtherDoneState() {
		return 1;
	}
	
	@Override
	double computeMinCost() throws CannotResolveException {
	
		if ( v._done == this.getDoneState() || v._done == 4) {
			return v._max_score;
		}
		//
	
		final int clusterSize = v.getCluster().getClusterSize();
	
		// SIA: base case for singelton clusters.
		if (clusterSize <= 1 || spm.isSingleSP(v.getCluster().getBitSet())) {
	
			v._max_score = scoreBaseCase(inference.isRooted(), inference.trees);
	
			v._min_lc = (v._min_rc = null);
			if(v._done == getOtherDoneState())
				v._done = 4;
			else
				v._done = getDoneState();
	
			return v._max_score;
		}
	
		//Iterable<VertexPair> clusterResolutions = this.inference.getClusterResolutions(this.v);
		final IClusterCollection containedVertecies = this.inference.dataCollection.clusters.getContainedClusters(v);
		
		final BlockingQueue<VertexPair> clusterResolutions  = new LinkedBlockingQueue<VertexPair>();
		
		if (clusterSize == GlobalMaps.taxonIdentifier.taxonCount()) {
			try {
			Vertex v1 = null;
			int smallestSize = -1;
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

			Vertex vert = v;
			Future<Integer> []  futures = new Future[clusterSize / 2];
			final Object lock = new Object();
			for (int j = 1; j <= (clusterSize / 2); j++) {
				final int i = j;
				futures[j-1] = Threading.submit(new Callable() {	
					public Integer call() {

						ArrayList<Set<Vertex>> clusters = ((HashClusterCollection)containedVertecies).clusters; 
						HashMap<Long, Vertex> h1ToVertexMap = ((HashClusterCollection)containedVertecies).h1ToVertexMap;

						if (clusters.get(i) == null || clusters.get(i).size() == 0 ||
								clusters.get(clusterSize - i) == null || clusters.get(clusterSize - i).size() == 0) {
							return 0;
								}
						Set<Vertex> small = (clusters.get(i).size() < clusters.get(clusterSize - i).size()) ? clusters.get(i) : clusters.get(clusterSize - i);

						for (Vertex v1 : small) {
							Vertex v2 = h1ToVertexMap.get(v.getCluster().hash1 - v1.getCluster().hash1);
							if (v2 != null) {
								if (v.getCluster().hash2 == (v1.getCluster().hash2 + v2.getCluster().hash2)
										//&& v.clusterSize == v1.clusterSize + v2.clusterSize // redundant check 
								   ) {
									if (v1.clusterSize != v2.clusterSize || v1.getCluster().hash1 < v2.getCluster().hash1) { // To avoid a pair of the same size twice
										VertexPair bi = new VertexPair(v1, v2, v);
										try {
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
						}
						return 0; 
					}
				});
			}
			for(int i = 0; i < futures.length; i++) {
				try {
					futures[i].get();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

		}
	
		long clusterLevelCost = 0;

		try {
			//(clusterResolutions).put(AbstractComputeMinCostTaskProducer.POISON_PILL);
			this.inference.getQueueClusterResolutions().put(clusterResolutions);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		for (VertexPair bi: clusterResolutions) {
			try {

				//if (bi == AbstractComputeMinCostTaskProducer.POISON_PILL) break;
				Vertex smallV = bi.cluster1;
				Vertex bigv = bi.cluster2;
				
				newMinCostTask(bigv).compute();
				newMinCostTask(smallV).compute();

				v._max_score = 0L;
				v._min_lc = smallV;
				v._min_rc = bigv;
				v._c = 0D ;

			} catch (Exception c) {
				c.printStackTrace();
				throw new RuntimeException("cannot resolve");
			}
		}
	
		
		if (v._done == getOtherDoneState())			
			v._done = 4;
		else
			v._done = getDoneState();
	
		return v._max_score;
	}

	@Override
	public Long getWeight(T t) {
		throw new RuntimeException("should not call");
	}

	

}
