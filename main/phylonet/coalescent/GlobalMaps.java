package phylonet.coalescent;

import java.util.Random;

/**
 * Global singelton classes
 * @author smirarab
 *
 */
public class GlobalMaps{

	/***
	 * Maps gene names to gene IDs (and vice versa)
	 */
	public static TaxonIdentifier taxonIdentifier = new TaxonIdentifier();
	/**
	 * Manages naming between gene and species names
	 */
	public static TaxonNameMap taxonNameMap;
	/**
	 * Random number generator
	 */
	public static Random random;

}