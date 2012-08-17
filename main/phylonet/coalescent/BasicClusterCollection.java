package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

public class BasicClusterCollection implements ClusterCollection {
	
	ArrayList<Set<Vertex>> clusters = new ArrayList<Set<Vertex>>();
	private HashMap<STITreeCluster, Vertex> clusterToVertx = new HashMap<STITreeCluster, Vertex>();
	int topClusterLength;
	int totalcount = 0;
	
	public BasicClusterCollection(int len) {
		this.topClusterLength = len;
		for (int i = 0; i <= len; i++) {
			clusters.add(new HashSet<Vertex>());
		}
	}

	@Override
	public Vertex getTopVertex() {
		Iterator<Vertex> it = clusters.get(topClusterLength).iterator();
		if (! it.hasNext()) {
			throw new NoSuchElementException();
		}
		return it.next();
	}

	@Override
	public int getClusterCount() {
		return totalcount;
	}
	int vertexIndex = 0;
	@Override
	public boolean addCluster(Vertex nv, int size) {
/*		if (!clusters.containsKey(size)) {
			clusters.put(size, new HashSet<Vertex>());
		}*/
		boolean added = clusters.get(size).add(nv);
		if (added) {
			nv.index=++vertexIndex;
			clusterToVertx.put(nv.getCluster(),nv);
		}
		if (added) totalcount++;
		return added;
	}

	@Override
	public boolean contains(Vertex vertex) {
		return clusters.get(vertex.getCluster().getClusterSize()).contains(vertex);
	}

	@Override
	public ClusterCollection getContainedClusters(STITreeCluster cluster) {
		//if (topClusterLength < 10)
		//	System.out.println("Contained: "+cluster+" "+clusterToVertx.keySet());
		ClusterCollection ret = new BasicClusterCollection(cluster.getClusterSize());
		int size = cluster.getClusterSize();
		ret.addCluster(getVertexForCluster(cluster), size);
		for (int i = size - 1 ; i > 0; i--) {
			Set<Vertex> sizeClusters = clusters.get(i);
			if (sizeClusters == null) continue;
			for (Vertex vertex : sizeClusters) {
				if (cluster.containsCluster(vertex.getCluster())) {
					ret.addCluster(vertex, i);
				}
			}
		}
		return ret;
	}

	@Override
	public Vertex getVertexForCluster(STITreeCluster cluster1) {
		return clusterToVertx.get(cluster1);
	}

	@Override
	public Collection<STBipartition> getClusterResolutions() {
		//System.out.println(topClusterLength+ " "+getTopVertex());
		STITreeCluster c = getTopVertex().getCluster();
		ArrayList<STBipartition> ret = new ArrayList<STBipartition>();
		int clusterSize = topClusterLength;
		for (int i = 1; i <= (clusterSize / 2); i++) {
			Set<Vertex> left = this.clusters.get(i);
			if (left == null || left.size() == 0) {
				continue;
			}
			Set<Vertex> right = this.clusters.get(clusterSize - i);
			if (right == null || right.size() == 0) {
				continue;
			}
			for (Vertex smallV : left) {
				
				for (Vertex bigv : right) {
					if (!smallV.getCluster().isDisjoint(bigv.getCluster())) {
						continue;
					}
					STBipartition bi = new STBipartition(
							smallV.getCluster(), bigv.getCluster(),
							c);
					ret.add(bi);
				}
			}
		}
		return ret;
	}

	@Override
	public Iterator<Set<Vertex>> getSubClusters() {
		return new Iterator<Set<Vertex>>() {
			int i = topClusterLength - 1;
			int next = topClusterLength ;
			@Override
			public boolean hasNext() {
				if (next > i) {
					next = i;
					while (clusters.get(next) == null || clusters.get(next).size() == 0) next--;
				}
				return next>0;
			}

			@Override
			public Set<Vertex> next() {
				if (! hasNext()) throw new NoSuchElementException();
				i = next;
				Set<Vertex> ret = clusters.get(i);
				i--;
				
				return ret;
			}

			@Override
			public void remove() {
				throw new NotImplementedException();
			}
		};
	}
	
/*	class BasicSubClusterCollection extends BasicClusterCollection {

		public BasicSubClusterCollection(int len) {
			super(len);
		}

		@Override
		public boolean addCluster(Vertex nv, int size) {
			if (!clusters.containsKey(size)) {
				clusters.put(size, new HashSet<Vertex>());
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
	Map<Integer, HashSet<Vertex>> maxSubClusters = 
		new HashMap<Integer, HashSet<Vertex>>();
	
	for (int i = 1; i < clusterSize; i++) {
		HashSet<Vertex> subClustersSizeI = containedVertecies.get(i);
		
		maxSubClusters.put(i, new HashSet<Vertex>());

		if (subClustersSizeI == null) {
			continue;
		}

		HashSet<Vertex> maxClustersSizeI = maxSubClusters.get(i);

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
