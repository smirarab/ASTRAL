package phylonet.tree.model.sti;

import phylonet.coalescent.GlobalMaps;
import phylonet.coalescent.TaxonIdentifier;
import phylonet.coalescent.IClusterCollection.VertexPair;
import phylonet.util.BitSet;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * A subset of a set of taxa
 * @author smirarab
 *
 */
public class STITreeCluster implements Iterable<Integer>
{
  //protected String[] _taxa;
  protected BitSet _cluster;
  private int hashCode = 0;
  
  /**
   * This identifies the meaning of the biset set. 
   * Bit number x in the bitset corresponds to taxon with ID x. 
   */
  private TaxonIdentifier taxonIdentifier;


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
//  public String[] getTaxa() {
//    return this._taxa;
//  }

  public void setCluster(BitSet c)
  {
    if (c != null) {
      this._cluster = c;
    }
    else
      System.err.println("Null bit set.");
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
  
//  public void addLeaf(String l)
//  {
//    addLeaf(GlobalMaps.taxonIdentifier.taxonId(l));
//  }
  
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
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }
    return this._cluster.equals(tc._cluster);
  }

//  /static HashMap<STITreeCluster,HashSet<STITreeCluster>> contains = new HashMap<STITreeCluster, HashSet<STITreeCluster>>();
  
  public int hashCode()
  { 
	  if (hashCode == 0)
		  hashCode =  this._cluster.hashCode() ;
	  return hashCode;
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

  public boolean isDisjoint(STITreeCluster tc)
  {
    if ((tc == null) || (tc._cluster == null)) {
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }

    return isDisjoint(tc._cluster);
  }

  public boolean isDisjoint(BitSet tc) { 
    return ! this._cluster.intersects(tc);
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
    STITreeCluster temp = new STITreeCluster(this);
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
	  	public long hash1 = 0, hash2 = 0;
	  	public int clusterSize = 0;
		public double _max_score = Integer.MIN_VALUE;
		public double _estimated = Integer.MIN_VALUE;
		public double _upper_bound = Integer.MAX_VALUE;
		public double _c = 0;
		public Vertex _min_lc = this._min_rc = null;
		public Vertex _min_rc;
		public List<Vertex> _subcl = null;	 // Don't matter	
		public byte _done = 0; // 0 for not, 1 for upper bound, 2 for estimated, 3 for yes
		public ArrayList<VertexPair> clusterResolutions = null;
		
		public Vertex() {
			super();
		}
		
		public String toString() {
			return STITreeCluster.this.toString() + "/" + this._max_score;
		}
		
		public STITreeCluster getCluster() {
			return STITreeCluster.this;
		}

		@Override
		public boolean equals(Object obj) {
			return ((Vertex) obj).getCluster().equals(STITreeCluster.this);
		}

		@Override
		public int hashCode() {
			return STITreeCluster.this.hashCode();
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
	@Override
	public Iterator<Integer> iterator() {
		return new ClusterIterator();
	}
}
