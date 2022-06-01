package phylonet.tree.model.sti;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import phylonet.coalescent.Factory;
import phylonet.coalescent.Logging;
import phylonet.coalescent.TaxonIdentifier;
import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.util.BitSet;


/**
 * A subset of a set of taxa
 * @author smirarab
 *
 */
public class STITreeCluster implements Iterable<Integer>
{
	//protected String[] _taxa;
	private BitSet _cluster;

	private int hashCode = 0;

	/**
	 * This identifies the meaning of the biset set. 
	 * Bit number x in the bitset corresponds to taxon with ID x. 
	 */
	public TaxonIdentifier taxonIdentifier;


	/*  public STITreeCluster()
  {
    this.taxonIdentifier = GlobalMaps.taxonIdentifier;
    this._cluster = new BitSet(this.taxonIdentifier.taxonCount());
  }*/

	public STITreeCluster(STITreeCluster tc)
	{
		this.taxonIdentifier = tc.taxonIdentifier;
		this._cluster = new BitSet(this.taxonIdentifier.taxonCount());
		this._cluster.or(tc._cluster);    
	}

	public STITreeCluster(TaxonIdentifier taxonId)
	{
		this.taxonIdentifier = taxonId;
		this._cluster = new BitSet(this.taxonIdentifier.taxonCount());
	}

	public STITreeCluster(TaxonIdentifier taxonId, BitSet bitset)
	{
		this.taxonIdentifier = taxonId;
		this._cluster = bitset;
	}

	public void setCluster(BitSet c)
	{
			this._cluster = c;
	}

	public BitSet getBitSet()
	{
		return this._cluster;
	}

	public int getClusterSize() {
		return this._cluster.cardinality();
	}

	public String[] getClusterLeaves()
	{
		String[] cl = new String[this._cluster.cardinality()];
		int c = 0;
		for (int i = 0; i < this._cluster.length(); i++) {
			if (this._cluster.get(i)) {
				cl[(c++)] = this.taxonIdentifier.getTaxonName(i);
			}
		}

		return cl;
	}

	public void addLeaf(int i){
		if (i < this.taxonIdentifier.taxonCount())
			this._cluster.set(i);
		else
			throw new RuntimeException(i +" above the length");
	}

	public void removeLeaf(String l)
	{
		this._cluster.clear(this.taxonIdentifier.taxonId(l));
	}

	public boolean equals(Object o)
	{
		if (!(o instanceof STITreeCluster)) {
			return false;
		}

		STITreeCluster tc = (STITreeCluster)o;
		if ((tc == null) || (tc._cluster == null)) {
			Logging.log("Cluster is null. The function returns false.");
			return false;
		}
		return this._cluster.equals(tc._cluster);
	}

	//  /static HashMap<STITreeCluster,HashSet<STITreeCluster>> contains = new HashMap<STITreeCluster, HashSet<STITreeCluster>>();

	@Override
	public int hashCode()
	{ 
		if (hashCode == 0)
			hashCode =  this._cluster.hashCode() ;
		return hashCode;
	}

	public boolean isCompatible(STITreeCluster tc)
	{
		if ((tc == null) || (tc._cluster == null)) {
			Logging.log("Cluster is null. The function returns false.");
			return false;
		}

		BitSet temp = (BitSet)this._cluster.clone();
		temp.and(tc._cluster);
		if ((temp.equals(this._cluster)) || (temp.equals(tc._cluster)))
		{
			return true;
		}

		return temp.cardinality() == 0;
	}

	public boolean isDisjoint(BitSet tc) { 
		return ! this._cluster.intersects(tc);
	}

	public boolean isDisjoint(STITreeCluster tc)
	{
		if ((tc == null) || (tc._cluster == null)) {
			Logging.log("Cluster is null. The function returns false.");
			return false;
		}

		return isDisjoint(tc._cluster);
	}



	public boolean isComplementary(STITreeCluster tc)
	{

		if ((tc == null) || (tc._cluster == null)) {
			Logging.log("Cluster is null. The function returns false.");
			return false;
		}

		BitSet temp1 = (BitSet)this._cluster.clone();
		temp1.and(tc._cluster);

		BitSet temp2 = (BitSet)this._cluster.clone();
		temp2.or(tc._cluster);

		return (temp1.cardinality() == 0) && (temp2.cardinality() == this.taxonIdentifier.taxonCount());
	}

	public boolean containsLeaf(String l)
	{
		return this._cluster.get(this.taxonIdentifier.taxonId(l));
	}

	public boolean containsCluster(STITreeCluster tc)
	{   
		//  return containsCluster(tc._cluster);
		/*if (contains.containsKey(this) && contains.get(this).contains(tc)) {
		return true;
	}*/
		boolean ret = this._cluster.contains(tc._cluster);
		/*if (ret) {
    	HashSet<STITreeCluster> hashSet = contains.containsKey(this) ?
    			contains.get(this):
    				new HashSet<STITreeCluster>();
    	hashSet.add(tc);
    	contains.put(this, hashSet);
    }*/
		return ret;
	}

	public boolean containsCluster(BitSet bs)
	{
		BitSet temp = (BitSet)bs.clone();
		temp.and(this._cluster);

		return temp.equals(bs);
	}

	public STITreeCluster merge(STITreeCluster tc)
	{
		STITreeCluster temp = Factory.instance.newCluster(this);
		temp._cluster.or(tc._cluster);

		return temp;
	}

	public STITreeCluster complementaryCluster() {
		STITreeCluster cc = new STITreeCluster(this.taxonIdentifier);
		BitSet bs = (BitSet)this._cluster.clone();
		bs.flip(0,this.taxonIdentifier.taxonCount());
		/*    for (int i = 0; i < this._taxa.length; i++) {
      if (bs.get(i)) {
        bs.set(i, false);
      }
      else {
        bs.set(i, true);
      }
    }
		 */    
		cc.setCluster(bs);
		return cc;
	}

	public String toString2()
	{
		StringBuffer out = new StringBuffer();

		out.append("{");
		for (String s : getClusterLeaves()) {
			out.append(s + ", ");
		}
		out.delete(out.length() - 2, out.length());
		out.append("} ");

		out.append(this._cluster);
		out.append(" ");
		//out.append(this.taxonIdentifier.getTaxonList());
		out.append(" ");
		out.append(this.taxonIdentifier.taxonCount());
		return out.toString();
	}
	public String toString()
	{
		StringBuffer out = new StringBuffer();

		out.append("{");
		for (String s : getClusterLeaves()) {
			out.append(s + ", ");
		}
		out.delete(out.length() - 2, out.length());
		out.append("}");

		return out.toString();
	}

	public void updateHash() {
	}
	
	public long partionId() {
		return getBitSet().nextSetBit(0);
	}
	
	@Override
	public Iterator<Integer> iterator() {
		return new ClusterIterator();
	}
	
	public Vertex newVertex() {
		return new VertexASTRAL3();
	}
	
	public Vertex newVertex(int size) {
		return new VertexASTRAL3();
	}

	/**
	 * A node in the dynamic programming. Thus, it includes
	 *   -- A cluster that we are trying to divide (the instance of the outer class)
	 *   -- A best resolution of this cluster into two cluster (min_lc and min_rc)
	 *   -- The score of the best resolution (_max_score)
	 *   -- Whether this nodes has been processed (_done)
	 * @author smirarab
	 *
	 */
	public class Vertex {
		//public STITreeCluster _cluster = null;
		//public int _el_num = -1;
		//public int _min_cost = -1;
		public long _max_score = Integer.MIN_VALUE;
		public Vertex _min_lc = this._min_rc = null;
		public Vertex _min_rc;
	

		@Override
		public boolean equals(Object obj) {
			return ((Vertex) obj).getCluster().equals(STITreeCluster.this);
		}

		public STITreeCluster getCluster() {
			return STITreeCluster.this;
		}

		@Override
		public int hashCode() {
			return STITreeCluster.this.hashCode();
		}

		public String toString() {
			return STITreeCluster.this.toString() + "/" + this._max_score;
		}


	}
	
	public class VertexASTRAL3 extends Vertex {
		
		private long _estimated = Integer.MIN_VALUE;
		private long _upper_bound = Integer.MAX_VALUE;
		private double _c = 0;
		private byte _done = 0; // 0 for not, 1 for upper bound, 2 for estimated, 3 for yes
		private ArrayList<VertexPair> clusterResolutions = null;

		
		public long get_estimated() {
			return _estimated;
		}

		public void set_estimated(long _estimated) {
			this._estimated = _estimated;
		}

		public long get_upper_bound() {
			return _upper_bound;
		}

		public void set_upper_bound(long _upper_bound) {
			this._upper_bound = _upper_bound;
		}


		public double get_c() {
			return _c;
		}

		public void set_c(double _c) {
			this._c = _c;
		}

		public byte isDone() {
			return _done;
		}

		public void setDone(byte _done) {
			this._done = _done;
		}
		
		public void setDone(int _done) {
			this._done = (byte) _done;
		}

		public ArrayList<VertexPair> getClusterResolutions() {
			return clusterResolutions;
		}

		public void setClusterResolutions(ArrayList<VertexPair> clusterResolutions) {
			this.clusterResolutions = clusterResolutions;
		}
		
	}
	

	class ClusterIterator implements  Iterator{

		int i = STITreeCluster.this._cluster.nextSetBit(0);
		@Override
		public boolean hasNext() {
			return i >= 0;
		}

		@Override
		public Integer next() {
			int r = i;
			i = STITreeCluster.this._cluster.nextSetBit(i+1);
			return r;
		}

		@Override
		public void remove() {
			throw new RuntimeException("Not supported");

		}

	}

}
