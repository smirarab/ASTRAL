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
	static cl_context context;
	static cl_device_id[] usedDevices;
	private static int numThreads = Runtime.getRuntime().availableProcessors();
	
	public static Future submit(Callable c) {
		return Threading.eService.submit(c);
	}
	
	public static void execute(Runnable c) {
		Threading.eService.execute(c);
	}

	public static void shutdown() {
		Threading.eService.shutdown();
		
	}

	public static void startThreading(int t) {
		Threading.numThreads = t;
		Threading.eService = Executors.newFixedThreadPool(Threading.numThreads);
		System.err.println("There are " + Threading.getNumThreads() + " threads used to run.");
		
	}

	public static int getNumThreads() {
		return numThreads;
	}

}
