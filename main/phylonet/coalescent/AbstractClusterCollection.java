package phylonet.coalescent;

import java.util.ArrayList;
import java.util.ConcurrentModificationException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex;

public abstract class AbstractClusterCollection implements IClusterCollection, Cloneable {
	static int countResolutions = 0;
	static int countContained = 0; 
	static long timeResolutions;
	static long timeContained;
	
	protected ArrayList<Set<Vertex>> clusters;
	protected int topClusterLength;
	int totalcount = 0;
	
	protected void initialize(int len) {
		this.topClusterLength = len;
		clusters = new ArrayList<Set<Vertex>>(len);
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

	@Override
	public boolean addCluster(Vertex vertex, int size) {
		boolean added;
		synchronized(clusters.get(size)) {
			added = clusters.get(size).add(vertex);
		}
		if (added) {
			totalcount++;
		}
		return added;
	}

	@Override
	public boolean contains(Vertex vertex) {
		return clusters.get(vertex.getCluster().getClusterSize()).contains(vertex);
	}

/*	public IClusterCollection getContainedClusters2(Vertex v) {
		if(!queue3done) {
			AbstractClusterCollection ret;
			try{
				ret = queue3.take();
				if(ret == CommandLine.POISON_PILL_queue3) {
					queue3done = true;
				}
				else {
					return ret;
				}
			}
			catch (Exception e) {
				
			}
				
		}
		return getContainedClusters(v);
	}*/
	@Override
	public IClusterCollection getContainedClusters(Vertex v) {
		long start = System.nanoTime();
		STITreeCluster cluster = v.getCluster();
		int size = cluster.getClusterSize();
		AbstractClusterCollection ret = newInstance(size);
		addClusterToRet(v, size, ret);
		CountDownLatch latch = new CountDownLatch(size - 1);
		for (int i = size - 1 ; i > 0; i--) {
			CommandLine.eService.execute(new getContainedClustersLoop(ret, cluster, i, latch));
		}
		try {
			latch.await();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		/*
		timeContained += System.nanoTime() - start;
		countContained++;
		if(countContained >= 10000) {
			countContained -= 10000;
			System.out.println("The time for this round of 10000 getContainedClusters is: " + (double)timeContained/1000000000);
			timeContained = 0;
		}
		*/

		return ret;
	}
	public class getContainedClustersLoop implements Runnable{
		AbstractClusterCollection ret;
		int i;
		STITreeCluster cluster;
		CountDownLatch latch;
		public getContainedClustersLoop(AbstractClusterCollection ret, STITreeCluster cluster, int i, CountDownLatch latch) {
			this.ret = ret;
			this.i = i;
			this.cluster = cluster;
			this.latch = latch;
		}
		public void run() {
			Set<Vertex> sizeClusters = clusters.get(i);
			if (sizeClusters == null) {
				latch.countDown();
				return;
			}
			for (Vertex vertex : sizeClusters) {
				if (cluster.containsCluster(vertex.getCluster())) {
					addClusterToRet(vertex, i, ret);
				}
			}
			latch.countDown();
		}
	}
	protected void addClusterToRet(Vertex vertex, int size, IClusterCollection ret) {
		ret.addCluster(vertex, size);	
	}

	@Override
	public Iterable<VertexPair> getClusterResolutions() {
		long start = System.nanoTime();
		//System.out.println(topClusterLength+ " "+getTopVertex());
		//TODO: return an iterator directly instead of building a collection.
		ArrayList<VertexPair> ret = new ArrayList<VertexPair>();
		/*Iterable<VertexPair> r= new Iterable<IClusterCollection.VertexPair>() {
			
			@Override
			public Iterator<VertexPair> iterator() {
				
				return new Iterator<VertexPair>() {

					@Override
					public boolean hasNext() {
						// TODO Auto-generated method stub
						return false;
					}

					@Override
					public VertexPair next() {
						// TODO Auto-generated method stub
						return null;
					}

					@Override
					public void remove() {
						throw new UnsupportedOperationException();
					}
				};
			}
		};*/
		
		int clusterSize = topClusterLength;
		Vertex v = this.getTopVertex();
		Future<ArrayList<VertexPair>>[] futures = new Future[clusterSize / 2];
		for (int i = 1; i <= (clusterSize / 2); i++) {
			futures[i - 1] = CommandLine.eService.submit(new getClusterResolutionsLoop(this.clusters, i, v, clusterSize));
		}
		for(int i = 0; i < futures.length; i++) {
			try {
				ArrayList<VertexPair> partialRet = (ArrayList<VertexPair>) futures[i].get();
				if(partialRet != null)
					ret.addAll(partialRet);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		
		}
		/*
		timeResolutions += System.nanoTime() - start;
		countResolutions++;
		if(countResolutions >= 10000) {
			countResolutions -= 10000;
			System.err.println("The time for this round of 10000 getClusterResolutions is: " + (double)timeResolutions/1000000000);
			timeResolutions = 0;
		}
		*/

		return ret;
	}
	public class getClusterResolutionsLoop implements Callable{
		int i;
		ArrayList<Set<Vertex>> clusters;
		Vertex v;
		int clusterSize;
		public getClusterResolutionsLoop(ArrayList<Set<Vertex>> cluster, int i, Vertex v, int clusterSize) {
			this.i = i;
			this.clusters = cluster;
			this.v = v;
			this.clusterSize = clusterSize;
		}
		public ArrayList<VertexPair> call() {
			ArrayList<VertexPair> ret = new ArrayList<VertexPair>();

			Set<Vertex> left = this.clusters.get(i);
			if (left == null || left.size() == 0) {
				return null;
			}
			Set<Vertex> right = this.clusters.get(clusterSize - i);
			if (right == null || right.size() == 0) {
				return null;
			}
			for (Vertex smallV : left) {
				
				for (Vertex bigv : right) {
					if (!smallV.getCluster().isDisjoint(bigv.getCluster())) {
						continue;
					}
					VertexPair bi = new VertexPair(
							smallV, bigv, v);
					ret.add(bi);
				}
			}
			return ret;
		}
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
				nset.add(v.getCluster().new Vertex());
			}
		}
		
		return clone;
	}
	

}
