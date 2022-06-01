package phylonet.coalescent;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.tree.model.sti.STITreeClusterMP;
import phylonet.util.BitSet;

public class HashOnlyTreeCluster extends STITreeClusterMP {

	static TaxonIdentifierMP taxid = (TaxonIdentifierMP) GlobalMaps.taxonIdentifier;
	public HashOnlyTreeCluster(){
		super(taxid);
		this.setCluster(null);
	}

	public HashOnlyTreeCluster(int i){
		super(taxid);
		this.hash1 = taxid.hash1[i];
		this.hash2 = taxid.hash2[i];
		this.setCluster(null);
	}

	public HashOnlyTreeCluster(long h1, long h2){
		super(taxid);
		this.hash1 = h1;
		this.hash2 = h2;
		this.setCluster(null);
	}

	public HashOnlyTreeCluster(STITreeClusterMP c){
		super((TaxonIdentifierMP) c.taxonIdentifier);
		c.updateHash();
		this.hash1 = c.hash1;
		this.hash2 = c.hash2;
		this.setCluster(null);
	}

	@Override
	public void updateHash(){}

	@Override
	public void setCluster(BitSet b)
	{
		if (b != null) {
			STITreeClusterMP c = new STITreeClusterMP(taxid, b);
			c.updateHash();
			this.hash1 = c.hash1;
			this.hash2 = c.hash2;
		}
		else
			super.setCluster(null);
	}

	@Override
	public final BitSet getBitSet()
	{
		Logging.log("Not applicable.");
		return null;
	}

	@Override
	public int getClusterSize() {
		Logging.log("Not applicable.");
		return 0;
	}

	@Override
	public String[] getClusterLeaves()
	{
		Logging.log("Not applicable.");
		return null;
	}

	//  public void addLeaf(String l)
	//  {
	//	    addLeaf(taxid.taxonId(l));
	//  }

	@Override
	public void addLeaf(int i){
		if (i < taxid.taxonCount()){
			hash1 += taxid.hash1[i];
			hash2 += taxid.hash2[i];
		}
		else
			throw new RuntimeException(i +" above the length");
	}

	@Override
	public void removeLeaf(String l)
	{
		int i = taxid.taxonId(l);
		hash1 -= taxid.hash1[i];
		hash2 -= taxid.hash2[i];
	}

	@Override
	public boolean equals(Object o)
	{
		if (!(o instanceof HashOnlyTreeCluster)) {
			return false;
		}

		HashOnlyTreeCluster tc = (HashOnlyTreeCluster)o;
		return this.hash1 == tc.hash1 && this.hash2 == tc.hash2;
	}

	//  /static HashMap<STITreeCluster,HashSet<STITreeCluster>> contains = new HashMap<STITreeCluster, HashSet<STITreeCluster>>();

	@Override
	public int hashCode()
	{
		return (int) hash1;
	}

	@Override
	public boolean isCompatible(STITreeCluster tc)
	{
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public boolean isDisjoint(STITreeCluster tc)
	{
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public boolean isDisjoint(BitSet tc) { 
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public boolean isComplementary(STITreeCluster tc)
	{
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public boolean containsLeaf(String l)
	{
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public boolean containsCluster(STITreeCluster tc)
	{   
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public boolean containsCluster(BitSet bs)
	{
		Logging.log("Not applicable.");
		return false;
	}

	@Override
	public STITreeClusterMP merge(STITreeClusterMP tc)
	{
		Logging.log("Not applicable.");
		return null;
	}

	@Override
	public STITreeClusterMP complementaryCluster() {
		Logging.log("Not applicable.");
		return null;
	}

	@Override
	public String toString2()
	{
		return "HashOnlyTreeCluster: " + this.toString();
	}

	@Override
	public String toString()
	{
		return this.hash1 + " | " + this.hash2;
	}

	@Override
	public HashOnlyTreeCluster clone() {
		return new HashOnlyTreeCluster(this);
	}

	public HashOnlyTreeCluster disjointClusterMerge(HashOnlyTreeCluster c){
		return new HashOnlyTreeCluster(this.hash1 + c.hash1, this.hash2 + c.hash2);
	}

	public HashOnlyTreeCluster subclusterComplement(HashOnlyTreeCluster c){
		return new HashOnlyTreeCluster(this.hash1 - c.hash1, this.hash2 - c.hash2);
	}
}
