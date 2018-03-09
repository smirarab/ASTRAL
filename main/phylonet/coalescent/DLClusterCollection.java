package phylonet.coalescent;

import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public class DLClusterCollection extends AbstractClusterCollection {
	
	protected HashMap<STITreeCluster, TreeSet<STBipartition>> geneTreeSTBByCluster;
	
	public DLClusterCollection(int len) {
		initialize(len);
		geneTreeSTBByCluster = new HashMap<STITreeCluster, TreeSet<STBipartition>>();
	}

	public void addGeneTreeSTB(STBipartition stb, int size) {
		if (! geneTreeSTBByCluster.containsKey(stb.c) ) {
			geneTreeSTBByCluster.put(stb.c, new TreeSet<STBipartition>());
		}
		geneTreeSTBByCluster.get(stb.c).add(stb);
		//System.err.println(geneTreeSTBByCluster);
	}
	
	@Override
	protected void addClusterToRet(Vertex vertex, int size, IClusterCollection ret) {
		super.addClusterToRet(vertex, size, ret);
		DLClusterCollection rett = (DLClusterCollection) ret;
		TreeSet<STBipartition> STBs = geneTreeSTBByCluster.get(vertex.getCluster());
		if (STBs != null) {
			rett.geneTreeSTBByCluster.put(vertex.getCluster(), STBs);
		}		
	}

	
	public Iterable<STBipartition> getContainedGeneTreeSTBs() {
		return new Iterable<STBipartition> () {

			@Override
			public Iterator<STBipartition> iterator() {
				return new Iterator<STBipartition>() {
					Iterator<TreeSet<STBipartition>> clustersIt = geneTreeSTBByCluster.values().iterator();
					Iterator<STBipartition> clusterSTBIt = null;
					//private STITreeCluster currentCluster;
					@Override
					public boolean hasNext() {
						if (clustersIt.hasNext()) {
							return true;
						}
						return clusterSTBIt !=null && clusterSTBIt.hasNext();
					}

					@Override
					public STBipartition next() {
						if (clusterSTBIt == null) {
							TreeSet<STBipartition> x = clustersIt.next();
							//System.out.println(x+  " is x " + geneTreeSTBByCluster);
							clusterSTBIt =x.iterator();
						}
						while (!clusterSTBIt.hasNext()) {
							clusterSTBIt = clustersIt.next().iterator();
						}
						
						return clusterSTBIt.next();
					}
					@Override
					public void remove() {
						throw new ConcurrentModificationException();
					}
				};
			}
			
		};
	}

	@Override
	public AbstractClusterCollection newInstance(int size) {
		return new DLClusterCollection(size);
	}
	
/*	class BasicSubClusterCollection extends BasicClusterCollection {

		public BasicSubClusterCollection(int len) {
			super(len);
		}

		@Override
		public boolean addCluster(Vertex nv, int size) {
			if (!clusters.containsKey(size)) {
				clusters.put(size, new TreeSet<Vertex>());
			}
			boolean added = clusters.get(size).add(nv);
			if (added) totalcount++;
			return added;
		}
		
		@Override
		public Vertex getVertexForCluster(STITreeCluster cluster1) {
			return BasicClusterCollection.this.getVertexForCluster(cluster1);
		}
	}
*/	
	/* The following code find all max-sub clusters (if ever needed)
	Map<Integer, TreeSet<Vertex>> maxSubClusters = 
		new HashMap<Integer, TreeSet<Vertex>>();
	
	for (int i = 1; i < clusterSize; i++) {
		TreeSet<Vertex> subClustersSizeI = containedVertecies.get(i);
		
		maxSubClusters.put(i, new TreeSet<Vertex>());

		if (subClustersSizeI == null) {
			continue;
		}

		TreeSet<Vertex> maxClustersSizeI = maxSubClusters.get(i);

		for (Vertex newMaxCluster : subClustersSizeI) {
			maxClustersSizeI.add(newMaxCluster);
		
			for (int j = i - 1; j > 0; j--) {
				
				List<Vertex> remove = new LinkedList<Vertex>();
				for (Vertex s : maxSubClusters.get(j)) {
					if (newMaxCluster._cluster.containsCluster(s._cluster)) {

						remove.add(s);
					}
				}
				maxSubClusters.get(j).removeAll(remove);
			}
		}

	}*/
	
}
