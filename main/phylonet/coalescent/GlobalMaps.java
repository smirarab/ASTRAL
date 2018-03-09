package phylonet.coalescent;

import java.util.HashMap;
import java.util.Random;
import java.util.concurrent.ExecutorService;

import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;

import phylonet.tree.model.sti.STITreeCluster.Vertex;

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
	/**
	 * Hash values for each taxon
	 */
	public static long[] hash1, hash2;
	
	static ExecutorService eService;
	static cl_context_properties contextProperties;
	static cl_context context;
	static cl_device_id[] usedDevices;
	static int numThreads = Runtime.getRuntime().availableProcessors();
	public static void logTimeMessage(String message) {
		if (GlobalMaps.timerOn) {
			System.err.println("TIME TOOK FROM LAST NOTICE " + message 
					+" "+ (double) (System.currentTimeMillis() - GlobalMaps.timer) / 1000);
			GlobalMaps.timer = System.currentTimeMillis();
		}
	}
	public static long timer;
	public static final boolean timerOn = false;
	public static void generateHashValues() {
		if (hash1 == null){
			Random rnd = GlobalMaps.random;
			int n = GlobalMaps.taxonIdentifier.taxonCount();
			long[] hash1 = new long[n], hash2 = new long[n];
			GlobalMaps.hash1 = hash1;
			GlobalMaps.hash2 = hash2;
			for (int i = 0; i < n; i++) {
				hash1[i] = rnd.nextLong();
				hash2[i] = rnd.nextLong();
			}
		}
	}
}