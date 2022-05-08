package phylonet.coalescent;

import static org.jocl.CL.CL_MEM_COPY_HOST_PTR;
import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_WRITE_ONLY;
import static org.jocl.CL.CL_TRUE;
import static org.jocl.CL.clBuildProgram;
import static org.jocl.CL.clCreateBuffer;
import static org.jocl.CL.clCreateCommandQueue;
import static org.jocl.CL.clCreateKernel;
import static org.jocl.CL.clCreateProgramWithSource;
import static org.jocl.CL.clEnqueueNDRangeKernel;
import static org.jocl.CL.clEnqueueReadBuffer;
import static org.jocl.CL.clEnqueueWriteBuffer;
import static org.jocl.CL.clSetKernelArg;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_program;

public class TurnTaskToScores implements Runnable {
	static final long THEEND = -23L;
	private static final int LOG_FREQ = 100000;
	private static final long workGroupSize = 1L << 14;
	private static final String clFileNVidia = "calculateWeightNVidia.cl";
	private static final String clFileAMD = "calculateWeightAMD.cl";
	
	private static final int cpuChunkSize = 64;
	public static final Object POISON_PILL = new Object();
	
	private static int nextLog = 0;
	private static long timer3;

	final PriorityBlockingQueue<ComparablePair<Long, Integer>> queue2Helper;
	final BlockingQueue<Long> queueWeightResults;
	final BlockingQueue<Tripartition> tripartitionsQueue;
	final AbstractWeightCalculatorMP<Tripartition> wqWeightCalculator;
	final WQInferenceConsumerMP inference;
	final WQDataCollectionMP dataCollection;
	
	AtomicInteger threadCount = new AtomicInteger(0);
	Object positionOutLock = new Object();
	
	private int positionOut = 0;

	
	int tripCounter = 0;
	//boolean done = false;
	
	final GPUCall gpuRunner;
	
	public final int speciesWordLength;
	

	public TurnTaskToScores(WQInferenceConsumerMP inf, BlockingQueue<Tripartition> queue1) {
		this.inference = inf;
		this.wqWeightCalculator = (AbstractWeightCalculatorMP<Tripartition>) inf.weightCalculator;
		this.dataCollection = (WQDataCollectionMP) inf.dataCollection;
		this.tripartitionsQueue = queue1;
		this.queueWeightResults = GlobalQueues.instance.getQueueWeightResults();
		
		this.queue2Helper = new PriorityBlockingQueue<ComparablePair<Long, Integer>>();
		this.speciesWordLength = (GlobalMaps.taxonIdentifier.taxonCount() / 64 + 1);
		
		//Logging.log("global work group size is : " + workGroupSize);
		if (Threading.usedDevices == null || Threading.usedDevices.length == 0) {
			this.gpuRunner = null;
		} else {
			this.gpuRunner = new GPUCall(((WQWeightCalculatorMP)this.wqWeightCalculator).geneTreesAsInts(), 
					((WQDataCollectionMP)inf.dataCollection).getAllArray(),
					inference, Threading.usedDevices, 
					Threading.context, Threading.contextProperties);
		}

	}

	private boolean noGPU() {
		return this.gpuRunner == null;
	}
	
	public void run() {
		int positionIn = 0;
		//long timer2 = System.nanoTime();
		//long timeWait = 0;
		Tripartition task = null;
		((WQWeightCalculatorMP) inference.weightCalculator).lastTime = System.currentTimeMillis();
		int currentGPU = 0;
		//long timer3 = System.nanoTime();
		Tripartition[] tripsForCPU = new Tripartition[cpuChunkSize];
		int[] tripsForCPULabel = new int[cpuChunkSize];
		int tripsForCPUCounter = 0;
		Logging.logTimeMessage(" TurnTaskToScores:92: ");

		Object taken = null; 
		while (true) {
			try{
				taken = this.tripartitionsQueue.take(); //tasks,(int)(2*workGroupSize));
			 } catch (Exception e) {
				 e.printStackTrace();
			}
			if (taken == POISON_PILL) {
				break;
			}
			task = (Tripartition) taken;
			if (noGPU() || (currentGPU == -1)) {
				tripsForCPULabel[tripsForCPUCounter] = positionIn++;
				tripsForCPU[tripsForCPUCounter++] = task;
				if (tripsForCPUCounter == cpuChunkSize) {
					Threading.execute(new CPUCalculationThread(tripsForCPU, tripsForCPULabel));
					tripsForCPU = new Tripartition[cpuChunkSize];
					tripsForCPULabel = new int[cpuChunkSize];
					tripsForCPUCounter = 0;
				}
			} 
			if (!noGPU()) {
				// Logging.logTimeMessage(" TurnTaskToScores:133: current GPU " + currentGPU);
				if (currentGPU == -1) {
					synchronized (gpuRunner.gpuLock) {
						for (int i = 0; i < Threading.usedDevices.length; i++) {
							if (gpuRunner.available[i].get()) {
								currentGPU = i;
								break;
							}
						}
						if (currentGPU == -1) {
							for (int i = 0; i < Threading.usedDevices.length; i++) {
								if (gpuRunner.storageAvailable[i].get()) {
									currentGPU = i;
									break;
								}
							}
						}
					}
				} else {
					gpuRunner.label[gpuRunner.currentLabel[currentGPU]][currentGPU][tripCounter] = positionIn++;

					for (int i = speciesWordLength - 1; i >= 0; i--){
							gpuRunner.tripartitions1[currentGPU][(speciesWordLength - i - 1) * (int) workGroupSize
						                               + tripCounter] = task.cluster1.getBitSet().words[i];
					}
					for (int i = speciesWordLength - 1; i >= 0; i--){
							gpuRunner.tripartitions2[currentGPU][(speciesWordLength - i - 1) * (int) workGroupSize
						                               + tripCounter] = task.cluster2.getBitSet().words[i];
					}
					for (int i = speciesWordLength - 1; i >= 0; i--){
							gpuRunner.tripartitions3[currentGPU][(speciesWordLength - i - 1) * (int) workGroupSize
						                               + tripCounter] = task.cluster3.getBitSet().words[i];
					}
					tripCounter++;
					//if (tripCounter % (workGroupSize/4) == 0) {
					//	Logging.logTimeMessage(" TurnTaskToScores:194: filling " + currentGPU + " "+tripCounter);
					//}
					if (tripCounter == workGroupSize) {
						Logging.logTimeMessage(" TurnTaskToScores:194: filled " + currentGPU + " "+ tripartitionsQueue.size() );
						synchronized (gpuRunner.gpuLock) {
							gpuRunner.storageAvailable[currentGPU].set(false);
							if (gpuRunner.available[currentGPU].get()) {
								gpuRunner.compute(workGroupSize, currentGPU);
								gpuRunner.available[currentGPU].set(false);
							}	
							tripCounter = 0;
							currentGPU = -1;
							for (int i = 0; i < Threading.usedDevices.length; i++) {
								if (gpuRunner.available[i].get()) {
									currentGPU = i;
									break;
								}
							}
							if (currentGPU == -1) {
								for (int i = 0; i < Threading.usedDevices.length; i++) {
									if (gpuRunner.storageAvailable[i].get()) {
										currentGPU = i;
										break;
									}
								}
							}
						}
						//Logging.logTimeMessage(" TurnTaskToScores:194: now to " + currentGPU );
					}									
				}
			}
		}
		//}
		// System.out.println("Time used to wait was: " + timeWait/1000);
		//TODO: Check if the call below is asynchronous. If not, it should be. 
		if (currentGPU != -1 && !noGPU()) {
			gpuRunner.compute(tripCounter, currentGPU);
		}
		if (tripsForCPUCounter != 0) {
			Threading.execute(new CPUCalculationThread(tripsForCPU, tripsForCPULabel, tripsForCPUCounter));
		}

		try {
			queue2Helper.offer(new ComparablePair<Long, Integer>(THEEND, positionIn++));
			// random  specific number used as a "poison pill" for AbstractWeightCalculator
			Logging.log("Finishing up:" + positionIn + " " + positionOut + " " + queue2Helper.peek().value.intValue());
		} catch (Exception e) {
			e.printStackTrace();
		}
		synchronized (positionOutLock) {
			while (!queue2Helper.isEmpty() && positionOut == queue2Helper.peek().value.intValue()
					&& queueWeightResults.add(queue2Helper.remove().key)) {
				positionOut++;
			}

			logWeights();
		}

		//GlobalMaps.logTimeMessage("Time used to wait on queue1.take() with at least one gpu available: "
		//		+ (double) (timeWait) / 1000000000+"\nTurnTaskToScores:199: ");

	}

	public void logWeights() {
		if (nextLog == 0) {
			timer3 = System.currentTimeMillis();
			nextLog += LOG_FREQ;
		}
		else if(positionOut > nextLog) {
			nextLog += LOG_FREQ;
			Logging.log(positionOut + " weights calculated " + 
					((double)System.currentTimeMillis() - timer3)/1000);
			timer3 = System.currentTimeMillis();
		}
	}

	public class GPUCall {
		public AtomicBoolean[] available;
		public AtomicBoolean[] storageAvailable;

		public long[][] tripartitions1;
		public long[][] tripartitions2;
		public long[][] tripartitions3;
		public long[][] weightArray;
		public int[][][] label;
		public int[] currentLabel;

		public int[] geneTreesAsInts;
		public short[] geneTreesAsShorts;
		public long[] allArray;
		public short[] stack;
		public long[] profile;
		public WQInferenceConsumerMP inference;

		public boolean DEBUG = false;

		private cl_context[] context;
		//private cl_context_properties contextProperties;
		private cl_kernel[] kernels;

		
		private cl_mem[] d_geneTreesAsInts;
		private cl_mem[] d_tripartitions1;
		private cl_mem[] d_tripartitions2;
		private cl_mem[] d_tripartitions3;
		private cl_mem[] d_allArray;
		private cl_mem[] d_weightArray;
		private cl_command_queue[] commandQueues;
		private cl_device_id[] devices;

		public Object gpuLock = new Object();
		
		public GPUCall(int[] geneTreesAsInts, long[] all, WQInferenceConsumerMP inference,
				cl_device_id[] devices, cl_context[] context, cl_context_properties contextProperties) {
			this.geneTreesAsInts = geneTreesAsInts;
			geneTreesAsShorts = new short[(geneTreesAsInts.length / 4 + 1) * 4];
			for (int i = 0; i < geneTreesAsInts.length; i++) {
				geneTreesAsShorts[i] = (short) geneTreesAsInts[i];
				if (geneTreesAsInts[i] == Integer.MIN_VALUE)
					geneTreesAsShorts[i] = Short.MIN_VALUE;
			}
			allArray = all;

			weightArray = new long[devices.length][(int) workGroupSize];
			tripartitions1 = new long[devices.length][(int) (speciesWordLength * workGroupSize)];
			tripartitions2 = new long[devices.length][(int) (speciesWordLength * workGroupSize)];
			tripartitions3 = new long[devices.length][(int) (speciesWordLength * workGroupSize)];
			label = new int[2][devices.length][(int) workGroupSize];
			currentLabel = new int[devices.length];
			Logging.log("device length is: " + devices.length);
			available = new AtomicBoolean[devices.length];
			storageAvailable = new AtomicBoolean[devices.length];

			for (int i = 0; i < devices.length; i++) {
				available[i] = new AtomicBoolean(true);
				storageAvailable[i] = new AtomicBoolean(false);
			}
			this.inference = inference;
			this.devices = devices;
			this.context = context;
			//this.contextProperties = contextProperties;

			initCL();

		}

		public void initCL() {
			int treeheight = ((WQWeightCalculatorMP) inference.weightCalculator).maxHeight();
			Logging.log("TREE HEIGHT IS: " + treeheight);
			
			//boolean NVidia = false;
			//boolean AMD = false;
			
			kernels = new cl_kernel[devices.length];
			for(int i = 0; i < devices.length; i++) {
				if(Threading.deviceVendors[i].toLowerCase().contains("nvidia")) {
					buildKernel(treeheight, true, i);
				}
				else {
					buildKernel(treeheight, false, i);
				}
			}
			
			// Program Setup. We want to avoid compiling twice. The true is for compiling nvidia code. The false if for compiling amd code.


			d_tripartitions1 = new cl_mem[devices.length];
			d_tripartitions2 = new cl_mem[devices.length];
			d_tripartitions3 = new cl_mem[devices.length];
			d_weightArray = new cl_mem[devices.length];
			commandQueues = new cl_command_queue[devices.length];
			d_geneTreesAsInts = new  cl_mem[devices.length];
			d_allArray = new  cl_mem[devices.length];

			for (int i = 0; i < devices.length; i++) {

				d_geneTreesAsInts[i] = clCreateBuffer(context[i], CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
					Sizeof.cl_short * geneTreesAsInts.length, Pointer.to(geneTreesAsShorts), null);
				d_allArray[i] = clCreateBuffer(context[i], CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY,
					Sizeof.cl_long * allArray.length, Pointer.to(allArray), null);
				d_tripartitions1[i] = clCreateBuffer(context[i], CL_MEM_READ_ONLY,
						Sizeof.cl_long * tripartitions1[i].length, null, null);
				d_tripartitions2[i] = clCreateBuffer(context[i], CL_MEM_READ_ONLY,
						Sizeof.cl_long * tripartitions2[i].length, null, null);
				d_tripartitions3[i] = clCreateBuffer(context[i], CL_MEM_READ_ONLY,
						Sizeof.cl_long * tripartitions3[i].length, null, null);
				d_weightArray[i] = clCreateBuffer(context[i], CL_MEM_WRITE_ONLY, Sizeof.cl_long * weightArray[i].length,
						null, null);
				commandQueues[i] = clCreateCommandQueue(context[i], devices[i], 0, null);
				
				clSetKernelArg(kernels[i], 0, Sizeof.cl_mem, Pointer.to(d_geneTreesAsInts[i]));
				clSetKernelArg(kernels[i], 1, Sizeof.cl_int, Pointer.to(new int[] { geneTreesAsInts.length }));
				clSetKernelArg(kernels[i], 2, Sizeof.cl_mem, Pointer.to(d_allArray[i]));

				
			}
		}

		public void buildKernel(int treeheight, boolean NVidia, int deviceIndex) {
			String source = "";
			if(NVidia) {
				source = readFile(getClass().getResourceAsStream(clFileNVidia));
			}
			else {
				source = readFile(getClass().getResourceAsStream(clFileAMD));
			}
						
			source = setupGPUSourceFile(treeheight, source);

			cl_program cpProgram = clCreateProgramWithSource(context[deviceIndex], 1, new String[] { source }, null, null);

			// Build the program
			if (DEBUG)
				clBuildProgram(cpProgram, 0, null, "-cl-opt-disable", null, null);
			else
				clBuildProgram(cpProgram, 0, null, "-cl-mad-enable -cl-strict-aliasing", null, null);

			// Create the kernel
			if(NVidia)
				kernels[deviceIndex] = clCreateKernel(cpProgram, "calcWeight", null);
			else
				kernels[deviceIndex] = clCreateKernel(cpProgram, "calcWeight", null);
			
		}

		public String setupGPUSourceFile(int treeheight, String source) {
			source = source.replaceAll("SPECIES_WORD_LENGTH - 1", Long.toString(speciesWordLength - 1));
			source = source.replaceAll("SPECIES_WORD_LENGTH", Long.toString(speciesWordLength));
			source = source.replaceAll("LONG_BIT_LENGTH", "64");
			source = source.replaceAll("STACK_SIZE", Integer.toString(treeheight + 2));
			source = source.replaceAll("TAXON_SIZE", Integer.toString(GlobalMaps.taxonIdentifier.taxonCount()));
			source = source.replaceAll("INT_MIN", "SHRT_MIN");
			source = source.replaceAll("WORK_GROUP_SIZE", Long.toString(workGroupSize));
			return source;
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

		public class ComputeThread implements Runnable {
			private long workSize;
			private int deviceIndex;

			public ComputeThread(long workSize, int deviceIndex) {
				this.workSize = workSize;
				this.deviceIndex = deviceIndex;
			}

			public void run() {
				// System.out.println("I am gone!!");

				do {
					// threadCount.incrementAndGet();
					Logging.logTimeMessage(" TurnTaskToScores:456: started " + deviceIndex );
					int labelIndex = currentLabel[deviceIndex];
					clEnqueueWriteBuffer(commandQueues[deviceIndex], d_tripartitions1[deviceIndex], CL_TRUE, 0L,
							Sizeof.cl_long * tripartitions1[deviceIndex].length,
							Pointer.to(tripartitions1[deviceIndex]), 0, null, null);
					//Logging.logTimeMessage(" TurnTaskToScores:542: started * " + deviceIndex );
					clEnqueueWriteBuffer(commandQueues[deviceIndex], d_tripartitions2[deviceIndex], CL_TRUE, 0L,
							Sizeof.cl_long * tripartitions2[deviceIndex].length,
							Pointer.to(tripartitions2[deviceIndex]), 0, null, null);
					clEnqueueWriteBuffer(commandQueues[deviceIndex], d_tripartitions3[deviceIndex], CL_TRUE, 0L,
							Sizeof.cl_long * tripartitions3[deviceIndex].length,
							Pointer.to(tripartitions3[deviceIndex]), 0, null, null);
					currentLabel[deviceIndex] = -labelIndex + 1;
					synchronized (gpuLock) {
						storageAvailable[deviceIndex].set(true);
					}
					cl_kernel  kernel = kernels[deviceIndex];
					synchronized (kernel) {
						clSetKernelArg(kernel, 3, Sizeof.cl_mem, Pointer.to(d_tripartitions1[deviceIndex]));
						clSetKernelArg(kernel, 4, Sizeof.cl_mem, Pointer.to(d_tripartitions2[deviceIndex]));
						clSetKernelArg(kernel, 5, Sizeof.cl_mem, Pointer.to(d_tripartitions3[deviceIndex]));
						clSetKernelArg(kernel, 6, Sizeof.cl_mem, Pointer.to(d_weightArray[deviceIndex]));
			
						clEnqueueNDRangeKernel(commandQueues[deviceIndex], kernel, 1, null, new long[] { workSize },
								null, 0, null, null);
					}
					clEnqueueReadBuffer(commandQueues[deviceIndex], d_weightArray[deviceIndex], CL_TRUE, 0L,
							Sizeof.cl_long * weightArray[deviceIndex].length, Pointer.to(weightArray[deviceIndex]), 0,
							null, null);
					for (int i = 0; i < workSize; i++) {
						queue2Helper.offer(
								new ComparablePair<Long, Integer>(weightArray[deviceIndex][i], label[labelIndex][deviceIndex][i]));
						// System.out.println("I have been here!" +
						// queue2Helper.peek().value + " " + positionOut);

					}
					synchronized (positionOutLock) {
						while (!queue2Helper.isEmpty() && positionOut == queue2Helper.peek().value.intValue()
								&& queueWeightResults.add(queue2Helper.remove().key)) {
							positionOut++;
						}
					}
					logWeights();
					synchronized (gpuLock) {
						Logging.logTimeMessage(" TurnTaskToScores:505:  finished " +deviceIndex + " "+ storageAvailable[deviceIndex].get() );
						if (storageAvailable[deviceIndex].get()) {
							available[deviceIndex].set(true);
							break;
						}
					}
				} while (true);
				// System.out.println(threadCount.decrementAndGet());
				// System.out.println("I am back!!!");
			}

		}

		public void compute(long workSize, int deviceNum) {

			Threading.execute(new ComputeThread(workSize, deviceNum));
		}
	}

	public class CPUCalculationThread implements Runnable {
		int[][] stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 2][3];

		int[][] overlap = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		int[][] overlapind = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		Tripartition[] trips;
		int[] positions;
		int numRuns = cpuChunkSize;
		
		Collection<Tripartition> nulls = Arrays.asList(new Tripartition[]{null});

		CPUCalculationThread(Tripartition [] trips, int[] positions) {
			this.trips = trips;
			this.positions = positions;
		}

		CPUCalculationThread(Tripartition [] trips, int[] positions, int numRuns) {

			this.positions = positions;
			this.numRuns = numRuns;
			this.trips = Arrays.copyOf(trips,numRuns);

		}

		public void run(){
			
			threadCount.incrementAndGet();
			
			Long[] weights = TurnTaskToScores.this.wqWeightCalculator.calculateWeight(trips);
				
			for(int i = 0; i < numRuns; i++) {
				queue2Helper.add(new ComparablePair<Long, Integer>(weights[i], positions[i]));
			}
			synchronized(positionOutLock) {
				while(!queue2Helper.isEmpty() && positionOut == queue2Helper.peek().value.intValue() && queueWeightResults.add(queue2Helper.remove().key)) {
					positionOut++;
				}
			}
			
			logWeights();

			threadCount.decrementAndGet();

		}

	}

}
