package phylonet.coalescent;

import java.util.Set;
import java.util.concurrent.CountDownLatch;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeCluster.Vertex; 

public abstract class AbstractClusterCollectionMP extends AbstractClusterCollection implements Cloneable {
	static int countResolutions = 0;
	static int countContained = 0; 
	static long timeResolutions;
	static long timeContained;


	@Override
	public int getClusterCount() {
		int s = 0;
		for (Set<Vertex> l:this.clusters) s += l.size();
		return s;
	}

	@Override
	public boolean addCluster(Vertex vertex, int size) {
		boolean added;
		synchronized(clusters.get(size)) {
			added = clusters.get(size).add((Vertex) vertex);
		}
		return added;
	}
	
	
	@Override
	public IClusterCollection getContainedClusters(Vertex v) {
		STITreeCluster cluster = v.getCluster();
		int size = cluster.getClusterSize();
		AbstractClusterCollection ret = newInstance(size);
		addClusterToRet(v, size, ret);
		CountDownLatch latch = new CountDownLatch(size - 1);
		for (int i = size - 1 ; i > 0; i--) {
			Threading.execute(newContainedClustersLoop(ret, (Vertex) v, i, latch));
		}
		try {
			latch.await();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}


		return ret;
	}
	
	public getContainedClustersLoop newContainedClustersLoop(AbstractClusterCollection ret, Vertex v, int i, CountDownLatch latch) {
		return new getContainedClustersLoop(ret, v, i, latch);
	}
	
	public class getContainedClustersLoop implements Runnable{
		AbstractClusterCollection ret;
		int i;
		//STITreeCluster cluster;
		CountDownLatch latch;
		Vertex v;
		
		public getContainedClustersLoop(AbstractClusterCollection ret, Vertex v, int i, CountDownLatch latch) {
			this.ret = ret;
			this.i = i;
			this.v = v;
			this.latch = latch;
		}
		public void run() {
			Set<Vertex> sizeClusters = clusters.get(i);
			if (sizeClusters == null) {
				latch.countDown();
				return;
			}
			STITreeCluster cluster = v.getCluster();
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



}
