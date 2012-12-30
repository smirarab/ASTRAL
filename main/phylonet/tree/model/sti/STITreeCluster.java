package phylonet.tree.model.sti;

import phylonet.util.BitSet;
import java.util.Arrays;
import java.util.List;


public class STITreeCluster
{
  protected String[] _taxa;
  protected BitSet _cluster;

  public STITreeCluster(String[] taxa)
  {
    if ((taxa == null) || (taxa.length == 0)) {
      System.err.println("Invalid cluster");

      this._taxa = null;
      this._cluster = null;
      return;
    }

    this._taxa = taxa;
    this._cluster = new BitSet(taxa.length);
  }

  public STITreeCluster(STITreeCluster tc)
  {
    assert ((tc._taxa != null) && (tc._taxa.length > 0));

    this._taxa = tc._taxa;
    this._cluster = new BitSet(this._taxa.length);
    this._cluster.or(tc._cluster);    
  }

  public String[] getTaxa() {
    return this._taxa;
  }

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
    assert ((this._taxa != null) && (this._taxa.length > 0));

    String[] cl = new String[this._cluster.cardinality()];
    int c = 0;
    for (int i = 0; i < this._cluster.length(); i++) {
      if (this._cluster.get(i)) {
        cl[(c++)] = this._taxa[i];
      }
    }

    return cl;
  }

  public void addLeaf(String l)
  {
    int i = -1;
    for (i = 0; i < this._taxa.length; i++) {
      if (l.equals(this._taxa[i]))
      {
        break;
      }
    }
    if (i < this._taxa.length)
      this._cluster.set(i);
    else
    	throw new RuntimeException(l +" not found in the taxon list");
  }
  
  public void removeLeaf(String l)
  {
    int i = -1;
    for (i = 0; i < this._taxa.length; i++) {
      if (l.equals(this._taxa[i]))
      {
        break;
      }
    }
    if (i < this._taxa.length)
      this._cluster.clear(i);
    else
    	throw new RuntimeException(l +" not found in the taxon list");
  }

  public boolean equals(Object o)
  {
    assert ((this._taxa != null) && (this._taxa.length > 0));
    
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
    return this._cluster.hashCode() + this._taxa.hashCode();
  }

  public boolean isCompatible(STITreeCluster tc)
  {
    assert ((this._taxa != null) && (this._taxa.length > 0));

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
    assert ((this._taxa != null) && (this._taxa.length > 0));

    if ((tc == null) || (tc._cluster == null)) {
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }

    return isDisjoint(tc._cluster);
  }

  public boolean isDisjoint(BitSet tc) {
    assert ((this._taxa != null) && (this._taxa.length > 0));

    if (tc == null) {
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }

    BitSet temp = (BitSet)this._cluster.clone();
    temp.and(tc);

    return temp.cardinality() == 0;
  }

  public boolean isComplementary(STITreeCluster tc)
  {
    assert ((this._taxa != null) && (this._taxa.length > 0));

    if ((tc == null) || (tc._cluster == null)) {
      System.err.println("Cluster is null. The function returns false.");
      return false;
    }

    BitSet temp1 = (BitSet)this._cluster.clone();
    temp1.and(tc._cluster);

    BitSet temp2 = (BitSet)this._cluster.clone();
    temp2.or(tc._cluster);

    return (temp1.cardinality() == 0) && (temp2.cardinality() == this._taxa.length);
  }

  public boolean containsLeaf(String l)
  {
	int i = 0;
    for (; i < this._taxa.length; i++) {
      if (this._taxa[i].equals(l))
      {
        break;
      }
    }
    return this._cluster.get(i);
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
    STITreeCluster cc = new STITreeCluster(this._taxa);
    BitSet bs = (BitSet)this._cluster.clone();
    bs.flip(0,this._taxa.length);
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
    out.append(this._taxa);
    out.append(" ");
    out.append(this._taxa.length);
    out.append(" ");
    out.append(Arrays.toString(this._taxa));
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
  
  public class Vertex {
		//public STITreeCluster _cluster = null;
		//public int _el_num = -1;
		//public int _min_cost = -1;
		public double _max_score = Integer.MIN_VALUE;
		public double _c = 0;
		public Vertex _min_lc = this._min_rc = null;
		public Vertex _min_rc;
		public List<Vertex> _subcl = null;		
		public byte _done = 0; // 0 for not, 1 for yes, 2 for failed
		
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
}
