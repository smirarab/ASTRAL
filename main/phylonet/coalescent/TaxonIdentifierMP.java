package phylonet.coalescent;

import java.util.Random;

public class TaxonIdentifierMP extends TaxonIdentifier {

	/**
	 * Hash values for each taxon
	 */
	public long[] hash1, hash2;


	public void lock() {
		super.lock();
		hash1 = new long[taxonCount]; 
		hash2 = new long[taxonCount];
		Random rnd = GlobalMaps.random;
		for (int i = 0; i < taxonCount; i++) {
			hash1[i] = rnd.nextLong() ^ rnd.nextLong() ^ rnd.nextLong(); 
			hash2[i] = ( rnd.nextLong() ^ rnd.nextLong() ) +  ( rnd.nextLong() ^ rnd.nextLong() );
			//System.err.println(idToName.get(i) +" "+hash1[i] +" "+ hash2[i]);
		}
	}
}
