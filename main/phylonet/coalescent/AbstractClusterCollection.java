package phylonet.coalescent;

import java.util.ArrayList;
import java.util.ConcurrentModificationException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;

public abstract class AbstractClusterCollection implements IClusterCollection, Cloneable {

	public ArrayList<Set<Vertex>> clusters;
	protected int topClusterLength;
	protected int totalcount = 0;
	public Vertex topV;

	protected void initialize(int len) {
		this.topClusterLength = len;
		clusters = new ArrayList<Set<Vertex>>(len);
		for (int i = 0; i <= len; i++) {
			clusters.add(new HashSet<Vertex>());
		}
	}
	
	@Override
	public Vertex getTopVertex() {
		if (topV == null) {
			Iterator<Vertex> it = clusters.get(topClusterLength).iterator();
			if (! it.hasNext()) {
				throw new NoSuchElementException();
			}
				topV = it.next();
		}
		return topV;
	}

	@Override
	public int getClusterCount() {
		return totalcount;
	}

	@Override
	public boolean addCluster(Vertex vertex, int size) {
		
		boolean added = clusters.get(size).add(vertex);
		if (added) {
			totalcount++;
		}
		return added;
	}
	
	@Override	
	public boolean removeCluster(Vertex vertex, int size) {
		
		boolean removed = clusters.get(size).remove(vertex);
		if(removed){
			totalcount--;
		}
		return removed;
	}

	@Override
	public boolean contains(Vertex vertex) {
		return clusters.get(vertex.getCluster().getClusterSize()).contains(vertex);
	}

	@Override
	public IClusterCollection getContainedClusters(Vertex v) {
		STITreeCluster cluster = v.getCluster();
		int size = cluster.getClusterSize();
		AbstractClusterCollection ret = newInstance(size);
		addClusterToRet(v, size, ret);
		
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
	public Iterable<VertexPair> getClusterResolutions() {
		//System.out.println(topClusterLength+ " "+getTopVertex());
		//TODO: return an iterator directly instead of building a collection.
		ArrayList<VertexPair> ret = new ArrayList<VertexPair>();

		
		int clusterSize = topClusterLength;
		Vertex v = this.getTopVertex();
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
					VertexPair bi = new VertexPair(
							(VertexASTRAL3) smallV, (VertexASTRAL3) bigv, (VertexASTRAL3) v);
					ret.add(bi);
				}
			}
		}
		return ret;
	}
	
	public abstract AbstractClusterCollection newInstance(int size);

	@Override
	public Set<Vertex> getSubClusters(int size) {
		return this.clusters.get(size);
	}
	
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
	
	public AbstractClusterCollection clone() throws CloneNotSupportedException {
		AbstractClusterCollection clone = (AbstractClusterCollection) super.clone();
		clone.clusters = new ArrayList<Set<Vertex>>();
		
		for (Set<Vertex> vset : this.clusters) {
			HashSet<Vertex> nset = new HashSet<STITreeCluster.Vertex>();
			clone.clusters.add(nset);
			for (Vertex v: vset) {
				nset.add(v.getCluster().newVertex());
			}
		}
		
		return clone;
	}
	

}
