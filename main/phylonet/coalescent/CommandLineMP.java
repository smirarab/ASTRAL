package phylonet.coalescent;

import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_NAME;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_DEVICE_VENDOR;
import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetDeviceInfo;
import static org.jocl.CL.clGetPlatformIDs;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_platform_id;

import com.martiansoftware.jsap.FlaggedOption;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import com.martiansoftware.jsap.Switch;

import phylonet.tree.model.Tree;


public class CommandLineMP extends CommandLine {

	protected CommandLineMP() {
		Factory.instance = new FactoryAstralMP();
		this.ASTRAL = "ASTRAL-MP";
	}

	protected Parameter[] getAstralParameters() {
		List<Parameter> parameters = new ArrayList(Arrays.asList(super.getAstralParameters()));
		List<Parameter> mpParameters = Arrays.asList(
				new Parameter[] { 

						new Switch("cpu only", 'C', "cpu-only", "Do not use GPUs."),
						new FlaggedOption("cpu threads", JSAP.INTEGER_PARSER, "-1", JSAP.NOT_REQUIRED, 'T',
								"cpu-threads", "Number of threads to use. Has to be at least 2. "),
						new FlaggedOption("GPU", JSAP.STRING_PARSER, null, JSAP.NOT_REQUIRED, 'G', "GPU",
								" the index of GPUs to be used, provided as a comma-separated list. If missing, all GPUs are used. (default)"),


						new FlaggedOption("matrixcount", JSAP.INTEGER_PARSER, null, JSAP.NOT_REQUIRED, JSAP.NO_SHORTFLAG,
								"matrixcount",
								"Sets the number of concurrent threads used for similarity matrix calculation. "
										+ " Reducing this can help reducing memory footprint. "),
				});
		parameters.addAll(mpParameters);
		return parameters.toArray(new Parameter[]{});
	}

	protected void exitWithErr(String extraMessage) {
		Logging.log("");
		Logging.log(extraMessage);
		Logging.log("");
		Logging.log("Usage: java -jar astralmp."+_version+".jar "+ jsap.getUsage());
		Logging.log("");
		Logging.log(jsap.getHelp());
		Threading.shutdown();
		System.exit( 1 );
	}
	
	@Override
	protected void setupComputing(JSAPResult config) {
		if (!config.getBoolean("cpu only")) {
			try {
				final int platformIndex = 0;
				final long deviceType = CL_DEVICE_TYPE_ALL;
				final int deviceIndex = 0;

				// Enable exceptions and subsequently omit error checks in this
				CL.setExceptionsEnabled(true);

				// Obtain the number of platforms
				int numPlatformsArray[] = new int[1];
				clGetPlatformIDs(0, null, numPlatformsArray);
				int numPlatforms = numPlatformsArray[0];
				// Obtain a platform ID
				cl_platform_id platforms[] = new cl_platform_id[numPlatforms];
				clGetPlatformIDs(platforms.length, platforms, null);
				cl_platform_id platform = platforms[platformIndex];
				// Initialize the context properties
				Threading.contextProperties = new cl_context_properties();
				Threading.contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);
				// Obtain the number of devices for the platform
				int numDevicesArray[] = new int[1];
				clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
				int numDevices = numDevicesArray[0];
				// Obtain a device ID
				cl_device_id devices[] = new cl_device_id[numDevices];
				clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
				Arrays.sort(devices, new Comparator<cl_device_id>() {

					@Override
					public int compare(cl_device_id arg0, cl_device_id arg1) {

						return getString(arg0, CL_DEVICE_NAME).compareTo(getString(arg1, CL_DEVICE_NAME));
					}
				});

				Logging.log("Detected GPU devices: ");
				for (int i = 0; i < numDevices; i++) {
					String deviceName = getString(devices[i], CL_DEVICE_NAME);
					String deviceVendor = getString(devices[i], CL_DEVICE_VENDOR);
					Logging.log("Device " + (i + 1) + " of " + numDevices + ": " + deviceName + " " + devices[i]
							+ " Vendor: " + deviceVendor);
				}

				ArrayList<cl_device_id> usedGPUs = new ArrayList<cl_device_id>();
				if (config.getString("GPU") != null) {
					try {
						for (String  si : config.getString("GPU").split(",")) {
							if (!usedGPUs.contains(devices[Integer.parseInt(si)-1]))
								usedGPUs.add(devices[Integer.parseInt(si)-1]);
						}
					} catch (Exception e) {
						exitWithErr("Could not parse GPU selection '" +config.getString("GPU")+"'. "
								+ "This should be comma-delimited list of integers between 1 and "+ numDevices+"\n"+e.toString());

					}
				} else {
					for (int i = 0; i < devices.length; i++) {
						usedGPUs.add(devices[i]);
					}
				}

				ArrayList<cl_device_id> usedDevicesAL = new ArrayList<cl_device_id>();
				ArrayList<String> deviceVendorsAL = new ArrayList<String>();
				for (cl_device_id d : usedGPUs) {
					deviceVendorsAL.add(getString(d, CL_DEVICE_VENDOR));	
					usedDevicesAL.add(d);
					Logging.log("Will use Device : " + getString(d, CL_DEVICE_NAME));
				}
				// testing only
				// deviceVendorsAL.add(getString(devices[0], CL_DEVICE_VENDOR));
				// usedDevicesAL.add(devices[0]);
				// usedDevicesAL.add(devices[1]);
				// usedDevicesAL.add(devices[2]);
				// usedDevicesAL.add(devices[3]);
				Threading.usedDevices = new cl_device_id[usedDevicesAL.size()];
				Threading.usedDevices = usedDevicesAL.toArray(Threading.usedDevices);
				Threading.deviceVendors = new String[deviceVendorsAL.size()];
				Threading.deviceVendors = deviceVendorsAL.toArray(Threading.deviceVendors);
				Threading.context = new  cl_context [Threading.usedDevices.length];
				for (int c = 0; c <  Threading.usedDevices.length; c++) {
					Threading.context[c] = clCreateContext(Threading.contextProperties, 1,
							new cl_device_id[] {Threading.usedDevices[c]}, null, null, null);
				}
				// johng23 end
			} catch (Exception e) {
				Logging.log("Warning:\n\n Problem using GPU. Proceeding without GPU\n"+e);
				Logging.log(e.toString());
				Threading.usedDevices = null;
			} catch (Error e) {
				Logging.log("Warning:\n\n Problem using GPU. Proceeding without GPU\n"+e);
				Logging.log(e.toString());
				Threading.usedDevices = null;
			}
		}
		int numThreads = config.getInt("cpu threads");
		if (numThreads == -1) {
			numThreads = Runtime.getRuntime().availableProcessors();
		}
		Logging.log("Starting " + numThreads + " threads ...");
		Threading.startThreading(numThreads);

	}

	protected Options readOptions( boolean rooted, boolean extrarooted, double wh,
			JSAPResult config, List<Tree> mainTrees, List<List<String>> bootstrapInputSets) 
					throws JSAPException, IOException {
		Options ret = super.readOptions(rooted, extrarooted, wh, config, mainTrees, bootstrapInputSets);

		if (config.contains("matrixcount"))
			Threading.setDistMatrixChunkSize(config.getInt("matrixcount"));
		else
			Threading.setDistMatrixChunkSize( Math.min(Threading.getNumThreads(), (10^9/GlobalMaps.taxonIdentifier.taxonCount()^2)/2 ) );
		return ret;
	}

	public void process(String[] args) throws Exception {
		try {
			super.process(args);
			Threading.shutdown();
		} catch (Exception e) {
			Threading.shutdown();
			throw (e);
		}
	}

	@Override
	protected void checkLibraries() {
		try {
			System.loadLibrary("Astral");
			Logging.log("Using native AVX batch computing.");
		} catch (Throwable e) {
			// e.printStackTrace();
			Logging.log("Warning: \n Fail to load native library " + System.mapLibraryName("Astral")
			+ "; use Java default computing method without AVX2, which is 4X slower. \n"
			+ " Make sure you are using the correct Djava.library.path (to the `lib` directory under ASTRAL where "
			+ System.mapLibraryName("Astral") + " can be found). \n"
			+ " Trying running make.sh. For mode debugging, run: java -Djava.library.path=lib/ -jar native_library_tester.jar");
		}
	}



	protected String getString(cl_device_id device, int paramName) {
		long size[] = new long[1];
		clGetDeviceInfo(device, paramName, 0, null, size);
		byte buffer[] = new byte[(int) size[0]];
		clGetDeviceInfo(device, paramName, buffer.length, Pointer.to(buffer), null);
		return new String(buffer, 0, buffer.length - 1);
	}

	public static void main(String[] args) throws Exception{
		CommandLineMP cm = new CommandLineMP();
		cm.process(args);
	}
}
