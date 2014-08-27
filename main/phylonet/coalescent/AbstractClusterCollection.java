package phylonet.coalescent;

import java.util.ArrayList;
import java.util.Collection;
import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractClusterCollection implements IClusterCollection{

	protected ArrayList<Set<Vertex>> clusters;
	HashMap<STITreeCluster, Vertex> globalVertexCash;// = new HashMap<STITreeCluster, STITreeCluster.Vertex>();
	protected int topClusterLength;
	int totalcount = 0;

	protected void initialize(int len, HashMap<STITreeCluster, Vertex> globalVertexCash) {
	    this.globalVertexCash = globalVertexCash;
		this.topClusterLength = len;
		 clusters = new ArrayList<Set<Vertex>>(len);
		for (int i = 0; i <= len; i++) {
			clusters.add(new HashSet<Vertex>());
			//geneTreeSTBBySize.add(new HashSet<STBipartition>());
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

	@Override
	public boolean addCluster(Vertex vertex, int size) {
		// See if this vertex is available in cash
		Vertex nv = globalVertexCash.get(vertex.getCluster());
		if (nv == null) {
			nv = vertex;
		}
		
		boolean added = clusters.get(size).add(nv);
		if (added) {
			//nv.index=++vertexIndex;
			globalVertexCash.put(nv.getCluster(),nv);
			totalcount++;
		}
		return added;
	}

	@Override
	public boolean contains(Vertex vertex) {
		return clusters.get(vertex.getCluster().getClusterSize()).contains(vertex);
	}

	@Override
	public IClusterCollection getContainedClusters(STITreeCluster cluster) {
		//if (topClusterLength < 10)
		//	System.out.println("Contained: "+cluster+" "+clusterToVertx.keySet());
		int size = cluster.getClusterSize();
		AbstractClusterCollection ret = newInstance(size);
		addClusterToRet(getVertexForCluster(cluster), size, ret);
		
		for (int i = size - 1 ; i > 0; i--) {
			Set<Vertex> sizeClusters = clusters.get(i);
			if (sizeClusters == null) continue;
			for (Vertex vertex : sizeClusters) {
				if (cluster.containsCluster(vertex.getCluster())) {
					addClusterToRet(vertex, i, ret);
				}
			}
		}
		return ret;
	}
	
	protected void addClusterToRet(Vertex vertex, int size, IClusterCollection ret) {
		ret.addCluster(vertex, size);	
	}

	@Override
	public Vertex getVertexForCluster(STITreeCluster cluster1) {
		return globalVertexCash.get(cluster1);
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
	
	public abstract AbstractClusterCollection newInstance(int size);

	@Override
	public Iterable<Set<Vertex>> getSubClusters() {
		return new Iterable<Set<Vertex>>() {
	
			@Override
			public Iterator<Set<Vertex>> iterator() {
	
				return new Iterator<Set<Vertex>>() {
					int i = topClusterLength - 1;
					int next = topClusterLength ;
					@Override
					public boolean hasNext() {
						if (next > i) {
							next = i;
							while (next>0 && (clusters.get(next) == null || clusters.get(next).size() == 0)) next--;
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
						throw new ConcurrentModificationException();
					}
				};
			}
		};
	}
	
	

}
