import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_MEM_COPY_HOST_PTR;
import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_WRITE_ONLY;
import static org.jocl.CL.CL_TRUE;
import static org.jocl.CL.clBuildProgram;
import static org.jocl.CL.clCreateBuffer;
import static org.jocl.CL.clCreateCommandQueue;
import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clCreateKernel;
import static org.jocl.CL.clCreateProgramWithSource;
import static org.jocl.CL.clEnqueueNDRangeKernel;
import static org.jocl.CL.clEnqueueReadBuffer;
import static org.jocl.CL.clEnqueueWriteBuffer;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetPlatformIDs;
import static org.jocl.CL.clSetKernelArg;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_platform_id;
import org.jocl.cl_program;

public class GPUCall {
	int[] tripartitions;
	int[] geneTreesAsInts;
	int[] allArray;
	long[] weightArray;
	
	public long WORK_GROUP_SIZE = 1L << 12;
	private cl_context context;
	private cl_command_queue commandQueue;
	private cl_kernel kernel;
	private cl_mem d_geneTreesAsInts;
	private cl_mem d_tripartitions;
	private cl_mem d_allArray;
	private cl_mem d_weightArray;
	
	public GPUCall (int[] geneTreesAsInts, int[] all, int[] trip) {
		this.geneTreesAsInts = geneTreesAsInts;
		allArray = all;
		tripartitions = trip;
		weightArray = new long[trip.length / 3 / Main.SPECIES_LENGTH];
		initCL();
	}
	private void initCL() {
		final int platformIndex = 0;
		final long deviceType = CL_DEVICE_TYPE_ALL;
		final int deviceIndex = 0;

		// Enable exceptions and subsequently omit error checks in this sample
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
		cl_context_properties contextProperties = new cl_context_properties();
		contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);

		// Obtain the number of devices for the platform
		int numDevicesArray[] = new int[1];
		clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
		int numDevices = numDevicesArray[0];

		// Obtain a device ID
		cl_device_id devices[] = new cl_device_id[numDevices];
		clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
		cl_device_id device = devices[deviceIndex];

		// Create a context for the selected device
		context = clCreateContext(contextProperties, 1, new cl_device_id[] { device }, null, null, null);

		// Create a command-queue for the selected device
		commandQueue = clCreateCommandQueue(context, device, 0, null);

		// Program Setup
		String source = readFile("main/gputest/calculateWeight.cl");
		source = source.replaceAll("SPECIES_LENGTH", "200");
		
		// Create the program
		cl_program cpProgram = clCreateProgramWithSource(context, 1, new String[] { source }, null, null);

		// Build the program
		clBuildProgram(cpProgram, 0, null, "-cl-mad-enable", null, null);

		// Create the kernel
		kernel = clCreateKernel(cpProgram, "calcWeight", null);

		// Create the memory object which will be filled with the
		// pixel data

		d_geneTreesAsInts = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY,
				Sizeof.cl_int * geneTreesAsInts.length, Pointer.to(geneTreesAsInts), null);
		d_allArray = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_int * allArray.length,
				Pointer.to(allArray), null);
		d_tripartitions = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_int * tripartitions.length,
				Pointer.to(tripartitions), null);
		d_weightArray = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_long * weightArray.length,
				Pointer.to(weightArray), null);



	}

	private String readFile(String fileName) {
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
			StringBuffer sb = new StringBuffer();
			String line = null;
			while (true) {
				line = br.readLine();
				if (line == null) {
					break;
				}
				sb.append(line).append("\n");
			}
			return sb.toString();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
			return null;
		}
	}

	public void update() {
		
		// Set work size and execute the kernel
		

		long globalWorkSize[] = new long[1];
		globalWorkSize[0] = WORK_GROUP_SIZE;
		
		clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(d_geneTreesAsInts));
		clSetKernelArg(kernel, 1, Sizeof.cl_int, Pointer.to(new int[] {geneTreesAsInts.length}));
		clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(d_allArray));
		clSetKernelArg(kernel, 3, Sizeof.cl_mem, Pointer.to(d_tripartitions));
		clSetKernelArg(kernel, 4, Sizeof.cl_mem, Pointer.to(d_weightArray));
		clEnqueueWriteBuffer(commandQueue, d_tripartitions, CL_TRUE, 0L, Sizeof.cl_int * tripartitions.length, Pointer.to(tripartitions), 0,
				null, null);
		
		int i = 0;
		for(; i < tripartitions.length / Main.SPECIES_LENGTH / 3 / WORK_GROUP_SIZE; i++) {
			clEnqueueNDRangeKernel(commandQueue, kernel, 1, new long[]{WORK_GROUP_SIZE * i}, globalWorkSize, null, 0, null, null);	
		}
		clEnqueueNDRangeKernel(commandQueue, kernel, 1, new long[]{WORK_GROUP_SIZE * i}, new long[] {(tripartitions.length / Main.SPECIES_LENGTH / 3) % WORK_GROUP_SIZE}, null, 0, null, null);	
 
		clEnqueueReadBuffer(commandQueue, d_weightArray, CL_TRUE, 0L, Sizeof.cl_long * weightArray.length, Pointer.to(weightArray), 0, null, null);
	}
}
