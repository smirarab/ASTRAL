package phylonet.tree.model.sti;

import java.util.Iterator;

import phylonet.coalescent.GlobalMaps;
import phylonet.coalescent.TaxonIdentifier;
import phylonet.coalescent.TaxonIdentifierMP;
import phylonet.tree.model.sti.STITreeCluster.Vertex;
import phylonet.util.BitSet;
import phylonet.util.BitSet.ImmutableBitSet;

/**
 * A subset of a set of taxa
 * @author smirarab
 *
 */
public class STITreeClusterMP extends STITreeCluster
{

	TaxonIdentifierMP taxonid =  (TaxonIdentifierMP) GlobalMaps.taxonIdentifier;
	public long hash1 = 0, hash2 = 0;

	/**
	 * This identifies the meaning of the biset set. 
	 * Bit number x in the bitset corresponds to taxon with ID x. 
	 */
	//TaxonIdentifierMP taxonIdentifier;

	/*public STITreeCluster()
  {
    this.taxonIdentifier = GlobalMaps.taxonIdentifier;
    this._cluster = new BitSet(this.taxonIdentifier.taxonCount());
  }*/

	public STITreeClusterMP(STITreeClusterMP tc)
	{
		super(tc);
		this.hash1 = tc.hash1;
		this.hash2 = tc.hash2; 
	}

	public STITreeClusterMP(TaxonIdentifier taxonId)
	{
		super(taxonId);
	}


	public STITreeClusterMP(TaxonIdentifier taxonId, BitSet bitset)
	{
		super(taxonId,bitset);
	}

	@Override
	public void setCluster(BitSet c)
	{
		super.setCluster(c);
		if (c != null) {
			this.hash1 = 0;
			this.hash2 = 0;
		}
	}


	@Override
	public void addLeaf(int i){
		try {
		if (i < this.taxonIdentifier.taxonCount() && this._cluster.get(i) == false){
			this._cluster.set(i);
			if (hash1 != 0){
				hash1 += taxonid.hash1[i];
				hash2 += taxonid.hash2[i];
			}
		}
		else
			throw new RuntimeException(i +" above the length");
		} catch (Exception e) {
			System.err.println(e);
		}
	}

	@Override
	public void removeLeaf(String l)
	{
		int i = this.taxonIdentifier.taxonId(l);
		if (this._cluster.get(i) == false) return;
		this._cluster.clear(i);
		if (hash1 != 0){
			hash1 -= taxonid.hash1[i];
			hash2 -= taxonid.hash2[i];
		}
	}

	/*public boolean equals(Object o)
	{
		if (!(o instanceof STITreeClusterMP)) {
			return false;
		}

		STITreeClusterMP tc = (STITreeClusterMP)o;
		if ((tc == null) || (tc._cluster == null)) {
			Logging.log("Cluster is null. The function returns false.");
			return false;
		}
		return this._cluster.equals(tc._cluster);
	}*/

	//  /static HashMap<STITreeCluster,HashSet<STITreeCluster>> contains = new HashMap<STITreeCluster, HashSet<STITreeCluster>>();

	@Override
	public int hashCode()
	{ 
		return _cluster.hashCode();
	}


	public STITreeClusterMP merge(STITreeClusterMP tc)
	{
		STITreeClusterMP temp = new STITreeClusterMP(this);
		temp._cluster.or(tc._cluster);
		temp.hash1 = 0;
		temp.hash2 = 0;
		return temp;
	}

	public STITreeClusterMP complementaryCluster() {
		STITreeClusterMP cc = new STITreeClusterMP( (TaxonIdentifierMP) this.taxonIdentifier);
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


	@Override
	public STITreeClusterMP clone() {
		return new STITreeClusterMP(this);
	}

	@Override
	public void updateHash()
	{
		if (hash1 != 0) return;
		hash1 = 0;
		hash2 = 0;
		BitSet b = getBitSet();	
		for (int k = b.nextSetBit(0); k >= 0; k = b.nextSetBit(k + 1)) {
			hash1 += taxonid.hash1[k];
			hash2 += taxonid.hash2[k];
		}
	}
	
	@Override
	public long partionId() {
		return this.hash1;
	}
	
	public Vertex newVertex() {
		throw new RuntimeException("invaild call"); // This is not to be used for ASTRAL-MP
	}
	
	@Override
	public Vertex newVertex(int size) {
		return new VertexMP(size);
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
	public class VertexMP  extends STITreeCluster.Vertex implements Comparable {

		public int clusterSize = 0;
		//long _max_score = Long.MIN_VALUE;
		
		//VertexMP _min_lc = null;
		//VertexMP _min_rc = null;

		private boolean _consDone = false;  
		private boolean _prodDone = false;

		private VertexMP(int size) {
			super();
			STITreeClusterMP.this._cluster = STITreeClusterMP.this._cluster.new ImmutableBitSet();
			this.clusterSize = size;
		}

		@Override
		public int compareTo(Object arg0) {
			return this.getCluster().getBitSet().compareTo(((VertexMP)arg0).getCluster().getBitSet());
		}

		/*public Vertex copy() {
				Vertex tmp = this.getCluster().new Vertex(this.clusterSize);
				//tmp.clusterSize = this.clusterSize;
				return tmp;
			}*/

		/*@Override
		public boolean equals(Object obj) {
			return ((Vertex) obj).getCluster().equals(STITreeClusterMP.this);
		}*/

		public STITreeClusterMP getCluster() {
			return STITreeClusterMP.this;
		}

		/*@Override
		public int hashCode() {
			return STITreeClusterMP.this.hashCode();
		}*/

		public String toString() {
			return super.toString() + " "+ (this.isProdDone()?"1":"0")+ (this.isConsDone()?"1":"0");
		}

		public boolean isConsDone() {
			return _consDone;
		}

		public void setConsDone(boolean _consDone) {
			this._consDone = _consDone;
		}

		public boolean isProdDone() {
			return _prodDone;
		}

		public void setProdDone(boolean _prodDone) {
			this._prodDone = _prodDone;
		}
	}
	
}
