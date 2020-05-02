package phylonet.tree.model.sti;

import java.util.Iterator;

import phylonet.coalescent.GlobalMaps;
import phylonet.coalescent.Logging;
import phylonet.coalescent.TaxonIdentifier;
import phylonet.util.BitSet;
import phylonet.util.BitSet.ImmutableBitSet;

/**
 * A subset of a set of taxa
 * @author smirarab
 *
 */
public class STITreeCluster implements Iterable<Integer>
{
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
  /**
   * A node in the dynamic programming. Thus, it includes
   *   -- A cluster that we are trying to divide (the instance of the outer class)
   *   -- A best resolution of this cluster into two cluster (min_lc and min_rc)
   *   -- The score of the best resolution (_max_score)
   *   -- Whether this nodes has been processed (_done)
   * @author smirarab
   *
   */
  public class Vertex  implements Comparable{

	  	public int clusterSize = 0;
		public long _max_score = Long.MIN_VALUE;

		public Vertex _min_lc = null;
		public Vertex _min_rc = null;
		public boolean _consDone = false;  
		public boolean _prodDone = false;
		
		public Vertex(int size) {
			super();
			STITreeCluster.this._cluster = STITreeCluster.this._cluster.new ImmutableBitSet();
			this.clusterSize = size;
		}
		
		@Override
		public int compareTo(Object arg0) {
			return this.getCluster().getBitSet().compareTo(((Vertex)arg0).getCluster().getBitSet());
		}
		
		/*public Vertex copy() {
			Vertex tmp = this.getCluster().new Vertex(this.clusterSize);
			//tmp.clusterSize = this.clusterSize;
			return tmp;
		}*/

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
  
  //protected String[] _taxa;
  private BitSet _cluster;
  
  public long hash1 = 0, hash2 = 0;
  
  /**
   * This identifies the meaning of the biset set. 
   * Bit number x in the bitset corresponds to taxon with ID x. 
   */
  TaxonIdentifier taxonIdentifier;

  /*public STITreeCluster()
  {
    this.taxonIdentifier = GlobalMaps.taxonIdentifier;
    this._cluster = new BitSet(this.taxonIdentifier.taxonCount());
  }*/
  
  public STITreeCluster(STITreeCluster tc)
  {
    this.taxonIdentifier = tc.taxonIdentifier;
    this.hash1 = tc.hash1;
    this.hash2 = tc.hash2;
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

  public void addLeaf(int i){
	  if (i < this.taxonIdentifier.taxonCount() && this._cluster.get(i) == false){
		  this._cluster.set(i);
		  if (hash1 != 0){
			  hash1 += GlobalMaps.taxonIdentifier.hash1[i];
			  hash2 += GlobalMaps.taxonIdentifier.hash2[i];
		  }
	  }
	    else
	    	throw new RuntimeException(i +" above the length");
  }

  
  public STITreeCluster complementaryCluster() {
    STITreeCluster cc = new STITreeCluster(this.taxonIdentifier);
    BitSet bs = (BitSet)this._cluster.copy();
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
    if (this._cluster instanceof ImmutableBitSet)
    	bs = bs.new ImmutableBitSet();
    cc.setCluster(bs);
    return cc;
  }

  public boolean containsCluster(BitSet bs)
  {
    BitSet temp = (BitSet)bs.clone();
    temp.and(this._cluster);

    return temp.equals(bs);
  }
  
//  public void addLeaf(String l)
//  {
//    addLeaf(GlobalMaps.taxonIdentifier.taxonId(l));
//  }
  
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
  
  public boolean containsLeaf(String l)
  {
    return this._cluster.get(this.taxonIdentifier.taxonId(l));
  }

  public boolean equals(Object o)
  {
    if (!(o instanceof STITreeCluster)) {
      return false;
    }

    STITreeCluster tc = (STITreeCluster)o;
    if ((tc == null) || (tc._cluster == null)) {
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }
    return this._cluster.equals(tc._cluster);
  }

//  /static HashMap<STITreeCluster,HashSet<STITreeCluster>> contains = new HashMap<STITreeCluster, HashSet<STITreeCluster>>();
  
  public BitSet getBitSet()
  {
    return this._cluster;
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

  public int getClusterSize() {
    return this._cluster.cardinality();
  }

  @Override
  public int hashCode()
  { 
	  return _cluster.hashCode();
  }

  public boolean isCompatible(STITreeCluster tc)
  {
    if ((tc == null) || (tc._cluster == null)) {
      System.err.println("Cluster is null. The function returns false.");
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

  public boolean isComplementary(STITreeCluster tc)
  {

    if ((tc == null) || (tc._cluster == null)) {
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }

    BitSet temp1 = (BitSet)this._cluster.clone();
    temp1.and(tc._cluster);

    BitSet temp2 = (BitSet)this._cluster.clone();
    temp2.or(tc._cluster);

    return (temp1.cardinality() == 0) && (temp2.cardinality() == this.taxonIdentifier.taxonCount());
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

  @Override
public Iterator<Integer> iterator() {
	return new ClusterIterator();
}

  public STITreeCluster merge(STITreeCluster tc)
  {
    STITreeCluster temp = new STITreeCluster(this);
    temp._cluster.or(tc._cluster);
    temp.hash1 = 0;
    temp.hash2 = 0;
    return temp;
  }

  public void removeLeaf(String l)
  {
      int i = this.taxonIdentifier.taxonId(l);
      if (this._cluster.get(i) == false) return;
	  this._cluster.clear(i);
	  if (hash1 != 0){
		  hash1 -= GlobalMaps.taxonIdentifier.hash1[i];
		  hash2 -= GlobalMaps.taxonIdentifier.hash2[i];
	  }
  }
  public void setCluster(BitSet c)
  {
    if (c != null) {
      this._cluster = c;
      this.hash1 = 0;
      this.hash2 = 0;
    }
    else
    	Logging.log("Null bit set.");
  }
  
  @Override
  public STITreeCluster clone() {
    return new STITreeCluster(this);
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
	//  public String[] getTaxa() {
	//    return this._taxa;
	//  }
	  public void updateHash()
	  {
		if (hash1 != 0) return;
		hash1 = 0;
		hash2 = 0;
		BitSet b = getBitSet();	
		for (int k = b.nextSetBit(0); k >= 0; k = b.nextSetBit(k + 1)) {
		  hash1 += GlobalMaps.taxonIdentifier.hash1[k];
		  hash2 += GlobalMaps.taxonIdentifier.hash2[k];
		}
	  }
}
