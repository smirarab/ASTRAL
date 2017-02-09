package phylonet.coalescent;

import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_LOCAL_MEM_SIZE;
import static org.jocl.CL.CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE;
import static org.jocl.CL.CL_DEVICE_NAME;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_MEM_COPY_HOST_PTR;
import static org.jocl.CL.CL_MEM_OBJECT_IMAGE2D;
import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_READ_WRITE;
import static org.jocl.CL.CL_MEM_WRITE_ONLY;
import static org.jocl.CL.CL_RGBA;
import static org.jocl.CL.CL_SIGNED_INT16;
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
import static org.jocl.CL.clGetDeviceInfo;
import static org.jocl.CL.clGetPlatformIDs;
import static org.jocl.CL.clSetKernelArg;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.concurrent.LinkedBlockingQueue;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_image_desc;
import org.jocl.cl_image_format;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_platform_id;
import org.jocl.cl_program;

public class TurnTaskToScores implements Runnable {
    public static final long workGroupSize = 1L<<13;
	public static final Object POISON_PILL = new Object();
    public static String clFile = "calculateWeight.cl";

	public LinkedBlockingQueue<Long> queue2;
	public LinkedBlockingQueue<Tripartition> queue1;
	public AbstractInference inference;
	public long[] tripartition1;
	public long[] tripartition2;
	public long[] tripartition3;
	
	public long[] all;
	public int tripCounter = 0;
	public boolean done = false;
	public GPUCall gpu;
	
	public final boolean pjohng23 = false;
	public int speciesWordLength;
	public TurnTaskToScores(AbstractInference inf, LinkedBlockingQueue<Tripartition> queue1, LinkedBlockingQueue<Long> queue2, int[] geneTreeAsInts, long[] all, int speciesWordLength) {
		this.inference = inf;
		this.queue1 = queue1;
		this.queue2 = queue2;
		this.all = all;
		this.speciesWordLength = speciesWordLength;
		tripartition1 = new long[(int)(speciesWordLength * workGroupSize)];
		tripartition2 = new long[(int)(speciesWordLength * workGroupSize)];
		tripartition3 = new long[(int)(speciesWordLength * workGroupSize)];
		System.err.println("global work group size is : " + workGroupSize);

		gpu = new GPUCall(geneTreeAsInts, all, tripartition1, tripartition2, tripartition3, inference, pjohng23);
		
	}

	public void run() {
		int writeCounter = 0;

		long timer = System.nanoTime();
		Object taken = null;
		Tripartition task = null;
		((WQWeightCalculator)inference.weightCalculator).lastTime = System.currentTimeMillis();
		while(true) {
			try {
				taken = queue1.take();
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			if(taken == POISON_PILL) {
				break;
			}
			task = (Tripartition)taken;
			for(int i = speciesWordLength - 1; i >= 0; i--)
				//tripartition1[tripCounter * SPECIES_WORD_LENGTH + SPECIES_WORD_LENGTH - i - 1] = task.trip.cluster1.getBitSet().words[i];
				tripartition1[(speciesWordLength - i - 1) * (int)workGroupSize + tripCounter]=task.cluster1.getBitSet().words[i];
			for(int i = speciesWordLength - 1; i >= 0; i--)
				//tripartition2[tripCounter * SPECIES_WORD_LENGTH + SPECIES_WORD_LENGTH - i - 1] = task.trip.cluster2.getBitSet().words[i];
				tripartition2[(speciesWordLength - i - 1) * (int)workGroupSize + tripCounter]=task.cluster2.getBitSet().words[i];
			for(int i = speciesWordLength - 1; i >= 0; i--)
				//tripartition3[tripCounter * SPECIES_WORD_LENGTH + SPECIES_WORD_LENGTH - i - 1] = task.trip.cluster3.getBitSet().words[i];
				tripartition3[(speciesWordLength - i - 1) * (int)workGroupSize + tripCounter]=task.cluster3.getBitSet().words[i];
			tripCounter++;
			if(tripCounter == workGroupSize) {
				gpu.compute(workGroupSize);
				tripCounter = 0;
				try {
					for(int i = 0; i < gpu.weightArray.length; i++) {
						queue2.put(gpu.weightArray[i]);
						writeCounter++;
					}
				}
				catch (Exception e) {
					e.printStackTrace();
				}
				timer = System.nanoTime();
			}
/*
			if((System.nanoTime() - timer) >= 15000000000L && tripCounter >= 128){
				
				gpu.compute(tripCounter);
				try {
					for(int i = 0; i < tripCounter; i++) {
						queue2.put(gpu.weightArray[i]);
					}
				}
				tripCounter = 0;
				catch (Exception e) {
					e.printStackTrace();
				}

			}
*/			
		}
		gpu.compute(tripCounter);
		try {
			for(int i = 0; i < tripCounter; i++) {
				queue2.put(gpu.weightArray[i]);
				writeCounter++;
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		System.err.println("THE WRITECOUNTER IS: " + writeCounter);

		if(pjohng23) {
			clEnqueueReadBuffer(gpu.commandQueue, gpu.d_profile, CL_TRUE, 0L, Sizeof.cl_long * gpu.profile.length, Pointer.to(gpu.profile), 0, null, null);
			long total = 0;	
			for(int i = 0; i < 4; i++) {
				total += gpu.profile[i];
			}
			System.out.println("intersecting with the all arrays takes: " + 4*(double)gpu.profile[0]/1000000000);
			System.out.println("adding numbers to the stack takes: " + 4*(double)gpu.profile[1]/1000000000);
			System.out.println("calculating weight of a tripartition takes: " + 4*(double)gpu.profile[2]/1000000000);
			System.out.println("calculating weight of a polytomy takes: " + 4*(double)gpu.profile[3]/1000000000);
		}
		try {
			queue2.put(-23L); //random specific number used as a "poison pill" for AbstractWeightCalculator
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		if(CommandLine.timerOn) {
	       	System.err.println("TIME TOOK FROM LAST NOTICE: " + (double)(System.nanoTime()-timer)/1000000000);
			timer = System.nanoTime();
		}
	}
	
	public class GPUCall {
		public long[] tripartitions1;
		public long[] tripartitions2;
		public long[] tripartitions3;
		
		public int[] geneTreesAsInts;
		public short[] geneTreesAsShorts;
		public long[] allArray;
		public long[] weightArray;
		public short[] stack;
		public long[] profile;
		public AbstractInference inference;
		
		public boolean p;
		private cl_context context;
		private cl_command_queue commandQueue;
		private cl_kernel kernel;
		private cl_mem d_geneTreesAsInts;
		private cl_mem[] d_geneTreesAsIntsConst;
		private cl_mem d_tripartitions1;
		private cl_mem d_tripartitions2;
		private cl_mem d_tripartitions3;
		private cl_mem d_allArray;
		private cl_mem d_weightArray;
		private cl_mem d_stack;
		private cl_mem d_profile;
		public GPUCall (int[] geneTreesAsInts, long[] all, long[] trip1, long[] trip2, long[] trip3, AbstractInference inference, boolean p) {
		this.p = p;
			this.geneTreesAsInts = geneTreesAsInts;
			geneTreesAsShorts = new short[(geneTreesAsInts.length/4 + 1) * 4];
			for(int i = 0; i < geneTreesAsInts.length; i++) {
				geneTreesAsShorts[i] = (short)geneTreesAsInts[i];
				if(geneTreesAsInts[i] == Integer.MIN_VALUE)
					geneTreesAsShorts[i] = Short.MIN_VALUE;
			}
			allArray = all;
			tripartitions1 = trip1;
			tripartitions2 = trip2;
			tripartitions3 = trip3;
			weightArray = new long[trip1.length / speciesWordLength];
			this.inference = inference;
			initCL();

			prepare();

		}
		public void initCL() {
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
			for (int i=0; i<numDevices; i++)
	        {
	            String deviceName = getString(devices[i], CL_DEVICE_NAME);
	            System.out.println("Device "+i+" of "+numDevices+": "+deviceName);
	        }
			cl_device_id device = devices[deviceIndex];

			// Create a context for the selected device
			context = clCreateContext(contextProperties, 1, new cl_device_id[] { device }, null, null, null);

			// Create a command-queue for the selected device
			commandQueue = clCreateCommandQueue(context, device, 0, null);

			//moved to here to edit source
			cl_image_format geneTreesImageFormat = new cl_image_format();
			geneTreesImageFormat.image_channel_data_type = CL_SIGNED_INT16;
			geneTreesImageFormat.image_channel_order = CL_RGBA;
			
			cl_image_desc geneTreesImageDesc = new cl_image_desc();
			geneTreesImageDesc.image_type = CL_MEM_OBJECT_IMAGE2D;
			geneTreesImageDesc.image_width = 4*((int)Math.sqrt((double)geneTreesAsShorts.length)/16);
			geneTreesImageDesc.image_height = geneTreesAsShorts.length/geneTreesImageDesc.image_width+1;
			geneTreesImageDesc.image_depth = 0;
			geneTreesImageDesc.image_array_size = 0;
			geneTreesImageDesc.image_row_pitch = 0;
			geneTreesImageDesc.image_slice_pitch = 0;
			geneTreesImageDesc.num_mip_levels = 0;
			geneTreesImageDesc.num_samples = 0;
			geneTreesImageDesc.buffer = null;
			
			// getting the tree height
			int treeheight = ((WQWeightCalculator)inference.weightCalculator).maxHeight();
			System.out.println("TREE HEIGHT IS: " + treeheight);
			// Program Setup
			String source = readFile(getClass().getResourceAsStream(clFile));
			// String source = readFile(clFile);
			source = source.replaceAll("SPECIES_WORD_LENGTH - 1", Long.toString(speciesWordLength - 1));
			source = source.replaceAll("SPECIES_WORD_LENGTH", Long.toString(speciesWordLength));
			source = source.replaceAll("LONG_BIT_LENGTH", "64");
//			source = source.replaceAll("(STACK_SIZE + 1) * 3", Integer.toString((treeheight + 1) * 3));
//			source = source.replaceAll("(STACK_SIZE + 2) * 3", Integer.toString((treeheight + 2) * 3));
//			source = source.replaceAll("(STACK_SIZE + 1) * 2", Integer.toString((treeheight + 1) * 2));
//			source = source.replaceAll("(STACK_SIZE + 2) * 2", Integer.toString((treeheight + 2) * 2));
//			source = source.replaceAll("(STACK_SIZE + 1)", Integer.toString(treeheight + 1));
//			source = source.replaceAll("(STACK_SIZE + 2)", Integer.toString(treeheight + 2));
			source = source.replaceAll("STACK_SIZE", Integer.toString(treeheight + 2));
			source = source.replaceAll("TAXON_SIZE", Integer.toString(GlobalMaps.taxonIdentifier.taxonCount()));
			source = source.replaceAll("INT_MIN", "SHRT_MIN");
			source = source.replaceAll("WORK_GROUP_SIZE", Long.toString(workGroupSize));
			source = source.replaceAll("IMAGE_WIDTH", Long.toString(geneTreesImageDesc.image_width));
			
			long[] localmemsize = new long[1];
			long[] constmemsize = new long[1];
			clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, 8, Pointer.to(localmemsize), null);
			clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, 8, Pointer.to(constmemsize), null);
			System.out.println("LOCAL MEMORY SIZE IS: " + localmemsize[0]);
			System.out.println("CONSTANT MEMORY SIZE IS: " + constmemsize[0]);
			// Create the program
			cl_program cpProgram = clCreateProgramWithSource(context, 1, new String[] { source }, null, null);

			// Build the program
			if(p)
				clBuildProgram(cpProgram, 0, null, "-cl-opt-disable", null, null);
			else
				clBuildProgram(cpProgram, 0, null, "-cl-mad-enable -cl-strict-aliasing", null, null);

			// Create the kernel
			kernel = clCreateKernel(cpProgram, "calcWeight", null);

			// Create the memory object which will be filled with the
			// pixel data
			
			d_geneTreesAsInts = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_short * geneTreesAsInts.length, Pointer.to(geneTreesAsShorts), null);
			d_allArray = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * allArray.length,
					Pointer.to(allArray), null);
			d_tripartitions1 = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * tripartitions1.length,
					Pointer.to(tripartitions1), null);
			d_tripartitions2 = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * tripartitions2.length,
					Pointer.to(tripartitions2), null);
			d_tripartitions3 = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * tripartitions3.length,
					Pointer.to(tripartitions3), null);
			
			//stack = new short[(int) (Sizeof.cl_ushort * 3 * (2 + treeheight) * workGroupSize)];
			d_stack = clCreateBuffer(context, CL_MEM_READ_WRITE , Sizeof.cl_long * 3 * (2 + treeheight) * workGroupSize,
					null, null);
			d_weightArray = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_long * weightArray.length,
					Pointer.to(weightArray), null);
		if(p){
				profile = new long[20];
				d_profile = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, Sizeof.cl_long * profile.length,
					Pointer.to(profile), null);
			}
		//d_c = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, Sizeof.cl_int,
			//		Pointer.to(c), null);

		}


		private String readFile(InputStream inputStream) {
			try {
				BufferedReader br = new BufferedReader(new InputStreamReader(inputStream));
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

		public void prepare() {

			
			clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(d_geneTreesAsInts));
			clSetKernelArg(kernel, 1, Sizeof.cl_int, Pointer.to(new int[] {geneTreesAsInts.length}));
			clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(d_allArray));
			clSetKernelArg(kernel, 3, Sizeof.cl_mem, Pointer.to(d_tripartitions1));
			clSetKernelArg(kernel, 4, Sizeof.cl_mem, Pointer.to(d_tripartitions2));
			clSetKernelArg(kernel, 5, Sizeof.cl_mem, Pointer.to(d_tripartitions3));
	   		clSetKernelArg(kernel, 6, Sizeof.cl_mem, Pointer.to(d_weightArray));

//		clSetKernelArg(kernel, 7, Sizeof.cl_mem, Pointer.to(d_stack));
		
		if(p)
	   			clSetKernelArg(kernel, 8, Sizeof.cl_mem, Pointer.to(d_profile));

//	 		clSetKernelArg(kernel, 8, Sizeof.cl_mem, Pointer.to(d_c));
	   		
//			clSetKernelArg(kernel, 7, Sizeof.cl_long * workGroupSize * SPECIES_WORD_LENGTH, null);
//			clSetKernelArg(kernel, 8, Sizeof.cl_long * workGroupSize * SPECIES_WORD_LENGTH, null);
//			clSetKernelArg(kernel, 9, Sizeof.cl_long * workGroupSize * SPECIES_WORD_LENGTH, null);

		}
		
		public void compute(long workSize) {
			
			// Set work size and execute the kernel
			
			clEnqueueWriteBuffer(commandQueue, d_tripartitions1, CL_TRUE, 0L, Sizeof.cl_long * tripartitions1.length, Pointer.to(tripartitions1), 0,
					null, null);
			clEnqueueWriteBuffer(commandQueue, d_tripartitions2, CL_TRUE, 0L, Sizeof.cl_long * tripartitions2.length, Pointer.to(tripartitions2), 0,
					null, null);
			clEnqueueWriteBuffer(commandQueue, d_tripartitions3, CL_TRUE, 0L, Sizeof.cl_long * tripartitions3.length, Pointer.to(tripartitions3), 0,
					null, null);
//			if(workSize >= 32 && workSize % 32 == 0) {
//				clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, new long[]{workSize}, new long[]{32L}, 0, null, null);
//			}
//			else	
				clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, new long[]{workSize}, null, 0, null, null);	

			clEnqueueReadBuffer(commandQueue, d_weightArray, CL_TRUE, 0L, Sizeof.cl_long * weightArray.length, Pointer.to(weightArray), 0, null, null);
			
		}
	}

	private static String getString(cl_device_id device, int paramName)
	{
	    long size[] = new long[1];
	    clGetDeviceInfo(device, paramName, 0, null, size);
	    byte buffer[] = new byte[(int)size[0]];
	    clGetDeviceInfo(device, paramName, 
	        buffer.length, Pointer.to(buffer), null);
	    return new String(buffer, 0, buffer.length-1);
	}	
}

