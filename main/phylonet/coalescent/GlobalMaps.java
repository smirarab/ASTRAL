package phylonet.coalescent;

import java.util.Random;
import java.util.concurrent.ExecutorService;

import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;

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
}