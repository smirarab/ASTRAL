package phylonet.coalescent;

import java.util.Set;

import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.tree.model.sti.STITreeCluster.VertexASTRAL3;

/***
 * Contains a set of clusters (used to implement X in ASTRAL/DynaDup). 
 * @author smirarab
 *
 */
public interface IClusterCollection {

	/**
	 * Every cluster in the clustercollection is a subset
	 * of the top cluster (by construction). 
	 * @return
	 */
	Vertex getTopVertex();

	int getClusterCount();

	boolean addCluster(Vertex nv, int size);
	
	boolean removeCluster(Vertex nv, int size);

	boolean contains(Vertex reverse);

	/**
	 * Returns a new IClusterCollection where the top node is v
	 * @param v
	 * @return
	 */
	IClusterCollection getContainedClusters(Vertex v);

	/**
	 * Returns all ways of dividing the top cluster into two subsets
	 * such that the two subsets are both part of this IClusterColleciton
	 * @return
	 */
	Iterable<VertexPair> getClusterResolutions();

	/**
	 * Returns a collection of sets of clusters contained in 
	 * this IClusterCollection, sorted by size from ? to ?
	 * @return
	 */
	Iterable<Set<Vertex>> getSubClusters();
	
	/***
	 * Returns the subclusters of a certain size. 
	 * @param size
	 * @return
	 */
	Set<Vertex> getSubClusters(int size);
	
	class VertexPair implements Comparable<VertexPair> {
		public final Vertex cluster1;
		public final Vertex cluster2;
		public final Vertex both;
		public long weight = 0, upperbound = 0;
		
		
		public VertexPair(Vertex v1, Vertex v2, Vertex v) {
			cluster1 = v1;
			cluster2 =  v2;
			both = v;
		}

		@Override
		public int compareTo(VertexPair o) {
			if (upperbound > o.upperbound) return -1;
			if (upperbound == o.upperbound) return 0;
			return 1;
		}
	}

	//IClusterCollection clone() throws CloneNotSupportedException;

}
