package phylonet.tree.model.sti;

import phylonet.tree.model.sti.STITreeCluster;
import phylonet.coalescent.GlobalMaps;
import phylonet.coalescent.TaxonIdentifier;
import phylonet.util.BitSet;

public class HashOnlyTreeCluster extends STITreeCluster {
	  public HashOnlyTreeCluster(){
		  super(GlobalMaps.taxonIdentifier);
	  }
	  
	  public HashOnlyTreeCluster(int i){
		super(GlobalMaps.taxonIdentifier);
		this.hash1 = GlobalMaps.taxonIdentifier.hash1[i];
		this.hash2 = GlobalMaps.taxonIdentifier.hash2[i];
	  }
	  
	  public HashOnlyTreeCluster(long h1, long h2){
		super(GlobalMaps.taxonIdentifier);
		this.hash1 = h1;
		this.hash2 = h2;
	  }
	  
	  public HashOnlyTreeCluster(STITreeCluster c){
		super(c.taxonIdentifier);
		c.updateHash();
		this.hash1 = c.hash1;
		this.hash2 = c.hash2;
	  }
	  
	  @Override
	  public void updateHash(){}

	  @Override
	  public void setCluster(BitSet b)
	  {
	    if (b != null) {
	      STITreeCluster c = new STITreeCluster(GlobalMaps.taxonIdentifier, b);
	      c.updateHash();
	      this.hash1 = c.hash1;
		  this.hash2 = c.hash2;
	    }
	    else
	      System.err.println("Null bit set.");
	  }

	  @Override
	  public final BitSet getBitSet()
	  {
		System.err.println("Not applicable.");
	    return null;
	  }

	  @Override
	  public int getClusterSize() {
		System.err.println("Not applicable.");
	    return 0;
	  }

	  @Override
	  public String[] getClusterLeaves()
	  {
		System.err.println("Not applicable.");
		return null;
	  }
	  
	//  public void addLeaf(String l)
	//  {
//	    addLeaf(GlobalMaps.taxonIdentifier.taxonId(l));
	//  }
	  
	  @Override
	  public void addLeaf(int i){
		if (i < GlobalMaps.taxonIdentifier.taxonCount()){
	      hash1 += GlobalMaps.taxonIdentifier.hash1[i];
		  hash2 += GlobalMaps.taxonIdentifier.hash2[i];
		}
		else
		  throw new RuntimeException(i +" above the length");
	  }
	  
	  @Override
	  public void removeLeaf(String l)
	  {
	    int i = GlobalMaps.taxonIdentifier.taxonId(l);
		hash1 -= GlobalMaps.taxonIdentifier.hash1[i];
		hash2 -= GlobalMaps.taxonIdentifier.hash2[i];
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
		System.err.println("Not applicable.");
		return false;
	  }
	  
	  @Override
	  public boolean isDisjoint(STITreeCluster tc)
	  {
		System.err.println("Not applicable.");
		return false;
	  }

	  @Override
	  public boolean isDisjoint(BitSet tc) { 
		System.err.println("Not applicable.");
		return false;
	  }
	  
	  @Override
	  public boolean isComplementary(STITreeCluster tc)
	  {
		System.err.println("Not applicable.");
		return false;
	  }
	  
	  @Override
	  public boolean containsLeaf(String l)
	  {
		System.err.println("Not applicable.");
		return false;
	  }
	  
	  @Override
	  public boolean containsCluster(STITreeCluster tc)
	  {   
		System.err.println("Not applicable.");
		return false;
	  }
	  
	  @Override
	  public boolean containsCluster(BitSet bs)
	  {
		System.err.println("Not applicable.");
		return false;
	  }

	  @Override
	  public STITreeCluster merge(STITreeCluster tc)
	  {
		System.err.println("Not applicable.");
		return null;
	  }

	  @Override
	  public STITreeCluster complementaryCluster() {
		System.err.println("Not applicable.");
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
