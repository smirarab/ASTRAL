package phylonet.coalescent;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;

public class Threading {

	private static ExecutorService eService;
	static cl_context_properties contextProperties;
	static cl_context [] context;
	static cl_device_id[] usedDevices;
	static String[] deviceVendors;
	private static int numThreads = Runtime.getRuntime().availableProcessors();
	private static int distMatrixChunks = -1;
	
	public static Future submit(Callable c) {
		return Threading.eService.submit(c);
	}
	
	public static void execute(Runnable c) {
		Threading.eService.execute(c);
	}

	public static void shutdown() {
		Logging.log("Shutting down threading");
		Logging.log(Logging.ENDMESSAGE);
		Threading.eService.shutdown();
	}

	public static void startThreading(int t) {
		Threading.numThreads = t;
		if (Threading.numThreads<2) {
			throw new RuntimeException("Sorry, at least two threads are needed. Switch to normal ASTRAL if you have only 1 thread"); 
		}
		Threading.eService = Executors.newFixedThreadPool(Threading.numThreads);
		Logging.log("There are " + Threading.getNumThreads() + " threads used to run.");
		
	}

	public static int getNumThreads() {
		return numThreads;
	}
	
	public static int getDistMatrixChunkSize() {
		return distMatrixChunks;
	}

	static void setDistMatrixChunkSize(int chunks) {
		distMatrixChunks  = chunks;
	}
}
