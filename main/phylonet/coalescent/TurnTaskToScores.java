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
import java.util.Iterator;
import java.util.concurrent.LinkedBlockingQueue;
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

import phylonet.tree.model.sti.STITreeCluster;

public class TurnTaskToScores implements Runnable {
	public static final long workGroupSize = 1L << 13;
	public static final Object POISON_PILL = new Object();
	public static String clFile = "calculateWeight.cl";
	public static final int cpuChunkSize = 10;

	public PriorityBlockingQueue<ComparablePair<Long, Integer>> queue2Helper;
	public LinkedBlockingQueue<Long> queue2;
	AtomicInteger threadCount = new AtomicInteger(0);
	Object positionOutLock = new Object();
	int positionOut = 0;
	int positionIn = 0;
	public LinkedBlockingQueue<Tripartition> queue1;
	public AbstractInference inference;
	public long[] all;
	public long timer3;
	public int tripCounter = 0;
	public boolean done = false;
	public GPUCall gpu;
	public cl_device_id[] devices;
	public int[] geneTreesAsInts;
	public WQDataCollection dataCollection;
	public final boolean pjohng23 = false;
	public int speciesWordLength;
	public boolean noGPU = false;
	public Object gpuLock = new Object();

	public TurnTaskToScores(AbstractInference inf, LinkedBlockingQueue<Tripartition> queue1,
			LinkedBlockingQueue<Long> queue2, int[] geneTreeAsInts, long[] all, int speciesWordLength,
			cl_device_id[] devices, cl_context context, cl_context_properties contextProperties) {
		this.inference = inf;
		this.queue1 = queue1;
		this.queue2 = queue2;
		this.geneTreesAsInts = geneTreeAsInts;
		queue2Helper = new PriorityBlockingQueue<ComparablePair<Long, Integer>>();
		this.all = all;
		this.speciesWordLength = speciesWordLength;
		System.err.println("global work group size is : " + workGroupSize);
		if (devices == null || devices.length == 0)
			noGPU = true;
		else {
			this.devices = devices;
			gpu = new GPUCall(geneTreeAsInts, all, inference, pjohng23, devices, context, contextProperties);
		}
		this.dataCollection = (WQDataCollection) inf.dataCollection;
	}

	public void run() {
		long timer2 = System.nanoTime();
		long timeWait = 0;
		Object taken = null;
		Tripartition task = null;
		((WQWeightCalculator) inference.weightCalculator).lastTime = System.currentTimeMillis();
		int currentGPU = 0;
		long timer3 = System.nanoTime();
		Tripartition[] tripsForCPU = new Tripartition[cpuChunkSize];
		int[] tripsForCPULabel = new int[cpuChunkSize];
		int tripsForCPUCounter = 0;
		if (CommandLine.timerOn) {
			System.err.println("TIME TOOK FROM LAST NOTICE TurnTaskToScores:92: "
					+ (double) (System.nanoTime() - CommandLine.timer) / 1000000000);
			CommandLine.timer = System.nanoTime();
		}
		while (true) {

			try {
				taken = queue1.take();
			} catch (Exception e) {
				e.printStackTrace();
			}
			if (taken == POISON_PILL) {
				break;
			}
			task = (Tripartition) taken;
			if (noGPU) {
				tripsForCPULabel[tripsForCPUCounter] = positionIn++;
				tripsForCPU[tripsForCPUCounter++] = task;
				if (tripsForCPUCounter == cpuChunkSize) {
					CommandLine.eService.execute(new CPUCalculationThread(tripsForCPU, tripsForCPULabel));
					tripsForCPU = new Tripartition[cpuChunkSize];
					tripsForCPULabel = new int[cpuChunkSize];
					tripsForCPUCounter = 0;
				}
			} else if (currentGPU == -1) {
				tripsForCPULabel[tripsForCPUCounter] = positionIn++;
				tripsForCPU[tripsForCPUCounter++] = task;
				if (tripsForCPUCounter == cpuChunkSize) {
					CommandLine.eService.execute(new CPUCalculationThread(tripsForCPU, tripsForCPULabel));
					tripsForCPU = new Tripartition[cpuChunkSize];
					tripsForCPULabel = new int[cpuChunkSize];
					tripsForCPUCounter = 0;
				}
				synchronized (gpuLock) {
					for (int i = 0; i < devices.length; i++) {
						if (gpu.available[i].get()) {
							currentGPU = i;
							break;
						}
					}
					if (currentGPU == -1) {
						for (int i = 0; i < devices.length; i++) {
							if (gpu.storageAvailable[i].get()) {
								currentGPU = i;
								break;
							}
						}
					}
				}
			} else {
				gpu.label[gpu.currentLabel[currentGPU]][currentGPU][tripCounter] = positionIn++;

				for (int i = speciesWordLength - 1; i >= 0; i--)
					gpu.tripartitions1[currentGPU][(speciesWordLength - i - 1) * (int) workGroupSize
							+ tripCounter] = task.cluster1.getBitSet().words[i];
				for (int i = speciesWordLength - 1; i >= 0; i--)
					gpu.tripartitions2[currentGPU][(speciesWordLength - i - 1) * (int) workGroupSize
							+ tripCounter] = task.cluster2.getBitSet().words[i];
				for (int i = speciesWordLength - 1; i >= 0; i--)
					gpu.tripartitions3[currentGPU][(speciesWordLength - i - 1) * (int) workGroupSize
							+ tripCounter] = task.cluster3.getBitSet().words[i];
				tripCounter++;
				if (tripCounter == workGroupSize) {
					synchronized (gpuLock) {
						if (gpu.available[currentGPU].get()) {
							gpu.compute(workGroupSize, currentGPU);
							gpu.available[currentGPU].set(false);
							gpu.storageAvailable[currentGPU].set(false);
						} else {
							gpu.storageAvailable[currentGPU].set(false);
						}
						tripCounter = 0;
						currentGPU = -1;
						for (int i = 0; i < devices.length; i++) {
							if (gpu.available[i].get()) {
								currentGPU = i;
								break;
							}
						}
						if (currentGPU == -1) {
							for (int i = 0; i < devices.length; i++) {
								if (gpu.storageAvailable[i].get()) {
									currentGPU = i;
									break;
								}
							}
						}
					}
				}
			}

		}
		// System.out.println("Time used to wait was: " + timeWait/1000000000);
		if (currentGPU != -1 && !noGPU) {
			gpu.compute(tripCounter, currentGPU);
		}
		if (tripsForCPUCounter != 0) {
			CommandLine.eService.execute(new CPUCalculationThread(tripsForCPU, tripsForCPULabel, tripsForCPUCounter));
		}

		try {
			queue2Helper.offer(new ComparablePair(-23L, positionIn++)); // random
																		// specific
																		// number
																		// used
																		// as a
																		// "poison
																		// pill"
																		// for
																		// AbstractWeightCalculator
			System.out.println(positionIn + " " + positionOut + " " + queue2Helper.peek().value.intValue());
		} catch (Exception e) {
			e.printStackTrace();
		}
		synchronized (positionOutLock) {
			while (!queue2Helper.isEmpty() && positionOut == queue2Helper.peek().value.intValue()
					&& queue2.add(queue2Helper.remove().key)) {
				positionOut++;
				if (positionOut % 100000 == 0) {
					System.out.println(positionOut + " weights calculated in turntasktoscores "
							+ ((double) System.nanoTime() - timer3) / 1000000000);
					timer3 = System.nanoTime();
				}
			}
		}
		// long timer2 = System.nanoTime();
		// System.out.println("WAAAAIIIT A MINUTE!");
		// System.err.println("THE WRITECOUNTER IS: " + positionOut);
		//
		// while(positionOut != positionIn) {
		// if(System.nanoTime() - timer2 > 1000000000) {
		// timer2 += 1000000000;
		// System.out.println("num threads: " + threadCount.get());
		// System.out.println("The queue is empty??" + queue2Helper.isEmpty());
		// System.out.println("Then lets check that the front of the queue
		// is..." + queue2Helper.peek().value + " and that positionout is" +
		// positionOut);
		// }
		// }

		if (CommandLine.timerOn) {
			System.err.println("Time used to wait on queue1.take() with at least one gpu available: "
					+ (double) (timeWait) / 1000000000);
			System.err.println("TIME TOOK FROM LAST NOTICE TurnTaskToScores:199: "
					+ (double) (System.nanoTime() - CommandLine.timer) / 1000000000);
			CommandLine.timer = System.nanoTime();
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
		public AbstractInference inference;

		public boolean p;

		private cl_context context;
		private cl_context_properties contextProperties;
		private cl_kernel kernel;
		private cl_mem d_geneTreesAsInts;
		private cl_mem[] d_tripartitions1;
		private cl_mem[] d_tripartitions2;
		private cl_mem[] d_tripartitions3;
		private cl_mem d_allArray;
		private cl_mem[] d_weightArray;
		private cl_command_queue[] commandQueues;
		private cl_device_id[] devices;

		public GPUCall(int[] geneTreesAsInts, long[] all, AbstractInference inference, boolean p,
				cl_device_id[] devices, cl_context context, cl_context_properties contextProperties) {
			this.p = p;
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
			System.out.println("device length is: " + devices.length);
			available = new AtomicBoolean[devices.length];
			storageAvailable = new AtomicBoolean[devices.length];

			for (int i = 0; i < devices.length; i++) {
				available[i] = new AtomicBoolean(true);
				storageAvailable[i] = new AtomicBoolean(false);
			}
			this.inference = inference;
			this.devices = devices;
			this.context = context;
			this.contextProperties = contextProperties;

			initCL();

		}

		public void initCL() {
			int treeheight = ((WQWeightCalculator) inference.weightCalculator).maxHeight();
			System.out.println("TREE HEIGHT IS: " + treeheight);
			// Program Setup
			String source = readFile(getClass().getResourceAsStream(clFile));
			// String source = readFile(clFile);
			source = source.replaceAll("SPECIES_WORD_LENGTH - 1", Long.toString(speciesWordLength - 1));
			source = source.replaceAll("SPECIES_WORD_LENGTH", Long.toString(speciesWordLength));
			source = source.replaceAll("LONG_BIT_LENGTH", "64");
			source = source.replaceAll("STACK_SIZE", Integer.toString(treeheight + 2));
			source = source.replaceAll("TAXON_SIZE", Integer.toString(GlobalMaps.taxonIdentifier.taxonCount()));
			source = source.replaceAll("INT_MIN", "SHRT_MIN");
			source = source.replaceAll("WORK_GROUP_SIZE", Long.toString(workGroupSize));

			cl_program cpProgram = clCreateProgramWithSource(context, 1, new String[] { source }, null, null);

			// Build the program
			if (p)
				clBuildProgram(cpProgram, 0, null, "-cl-opt-disable", null, null);
			else
				clBuildProgram(cpProgram, 0, null, "-cl-mad-enable -cl-strict-aliasing", null, null);

			// Create the kernel
			kernel = clCreateKernel(cpProgram, "calcWeight", null);

			d_tripartitions1 = new cl_mem[devices.length];
			d_tripartitions2 = new cl_mem[devices.length];
			d_tripartitions3 = new cl_mem[devices.length];
			d_weightArray = new cl_mem[devices.length];
			commandQueues = new cl_command_queue[devices.length];

			d_geneTreesAsInts = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
					Sizeof.cl_short * geneTreesAsInts.length, Pointer.to(geneTreesAsShorts), null);
			d_allArray = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY,
					Sizeof.cl_long * allArray.length, Pointer.to(allArray), null);
			for (int i = 0; i < devices.length; i++) {

				d_tripartitions1[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
						Sizeof.cl_long * tripartitions1[i].length, null, null);
				d_tripartitions2[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
						Sizeof.cl_long * tripartitions2[i].length, null, null);
				d_tripartitions3[i] = clCreateBuffer(context, CL_MEM_READ_ONLY,
						Sizeof.cl_long * tripartitions3[i].length, null, null);
				d_weightArray[i] = clCreateBuffer(context, CL_MEM_WRITE_ONLY, Sizeof.cl_long * weightArray[i].length,
						null, null);
				commandQueues[i] = clCreateCommandQueue(context, devices[i], 0, null);
				clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(d_geneTreesAsInts));
				clSetKernelArg(kernel, 1, Sizeof.cl_int, Pointer.to(new int[] { geneTreesAsInts.length }));
				clSetKernelArg(kernel, 2, Sizeof.cl_mem, Pointer.to(d_allArray));
			}
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
					int labelIndex = currentLabel[deviceIndex];
					clEnqueueWriteBuffer(commandQueues[deviceIndex], d_tripartitions1[deviceIndex], CL_TRUE, 0L,
							Sizeof.cl_long * tripartitions1[deviceIndex].length,
							Pointer.to(tripartitions1[deviceIndex]), 0, null, null);
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
								new ComparablePair(weightArray[deviceIndex][i], label[labelIndex][deviceIndex][i]));
						// System.out.println("I have been here!" +
						// queue2Helper.peek().value + " " + positionOut);

					}
					synchronized (positionOutLock) {
						while (!queue2Helper.isEmpty() && positionOut == queue2Helper.peek().value.intValue()
								&& queue2.add(queue2Helper.remove().key)) {
							positionOut++;
							if (positionOut % 100000 == 0) {
								System.out.println(positionOut + " weights calculated in turntasktoscores "
										+ ((double) System.nanoTime() - timer3) / 1000000000);
							}
							timer3 = System.nanoTime();
						}
					}
					synchronized (gpuLock) {
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
			CommandLine.eService.execute(new ComputeThread(workSize, deviceNum));
		}
	}

	public class CPUCalculationThread implements Runnable {
		int[][] stack = new int[GlobalMaps.taxonIdentifier.taxonCount() + 2][3];

		int[][] overlap = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		int[][] overlapind = new int[GlobalMaps.taxonIdentifier.taxonCount() + 1][3];
		Tripartition[] trips;
		int[] positions;
		int numRuns = cpuChunkSize;

		CPUCalculationThread(Tripartition[] trips, int[] positions) {
			this.trips = trips;
			this.positions = positions;
		}

		CPUCalculationThread(Tripartition[] trips, int[] positions, int numRuns) {
			this.trips = trips;
			this.positions = positions;
			this.numRuns = numRuns;
		}

		long F(int a, int b, int c) {
			if (a < 0 || b < 0 || c < 0) {
				throw new RuntimeException("negative side not expected: " + a + " " + b + " " + c);
			}
			long ret = (a + b + c - 3);
			ret *= a * b * c;
			return ret;
		}

		public void run(){

			threadCount.incrementAndGet();
			int[] allsides = null;
			long[] weights = new long[numRuns];
			for(int ii = 0; ii < numRuns; ii++) {
				Iterator<STITreeCluster> tit = dataCollection.treeAllClusters.iterator();
				Tripartition trip = trips[ii];
				boolean newTree = true;
				int top = 0; // The first empty place on stack (generally)
				for (Integer gtb : geneTreesAsInts) {
					if (newTree) {
						STITreeCluster all = tit.next();
						allsides = new int[] {
								trip.cluster1.getBitSet().intersectionSize(all.getBitSet()),
								trip.cluster2.getBitSet().intersectionSize(all.getBitSet()),
								trip.cluster3.getBitSet().intersectionSize(all.getBitSet())};
						newTree = false;
					}
					if (gtb >= 0) { // Leaf nodes
						if (trip.cluster1.getBitSet().get(gtb)) {
							stack[top][0] = 1;
							stack[top][1] = 0;
							stack[top][2] = 0;
						} else if (trip.cluster2.getBitSet().get(gtb)) {
							stack[top][0] = 0;
							stack[top][1] = 1;
							stack[top][2] = 0;
						} else if (trip.cluster3.getBitSet().get(gtb)) {
							stack[top][0] = 0;
							stack[top][1] = 0;
							stack[top][2] = 1;
						} else { // This can happen due to missing data
							stack[top][0] = 0;
							stack[top][1] = 0;
							stack[top][2] = 0;
						}
						top++;
					} else if (gtb == Integer.MIN_VALUE) { // delimiter between trees
						top = 0;
						newTree = true;
					} else if (gtb == -2) { // Internal nodes
		
						top--;
						int newSides0 = stack[top][0] + stack[top - 1][0];
						int newSides1 = stack[top][1] + stack[top - 1][1];
						int newSides2 = stack[top][2] + stack[top - 1][2];
		
						int side3s0 = allsides[0] - newSides0;
						int side3s1 = allsides[1] - newSides1;
						int side3s2 = allsides[2] - newSides2;
		
						weights[ii] += F(stack[top][0], stack[top - 1][1], side3s2)
								+ F(stack[top][0], stack[top - 1][2], side3s1)
								+ F(stack[top][1], stack[top - 1][0], side3s2)
								+ F(stack[top][1], stack[top - 1][2], side3s0)
								+ F(stack[top][2], stack[top - 1][0], side3s1)
								+ F(stack[top][2], stack[top - 1][1], side3s0);
		
						stack[top - 1][0] = newSides0;
						stack[top - 1][1] = newSides1;
						stack[top - 1][2] = newSides2;
					} else { // The following case is relevant only for polytomies.
						int[] sxy = new int[3];
						long tempWeight = 0;

                        int[] newSides = new int[3];

                        for(int side = 0; side < 3; side++) {
                                stack[top][side] = allsides[side];
                                for(int i = top - 1; i>= top + gtb; i--) {
                                        stack[top][side] -= stack[i][side];
                                }
                        }
                        for(int i = top; i>= top + gtb; i--) {
                        	for(int side = 0; side < 3; side++) {
                                newSides[side] = stack[i][side];
                            }
                            sxy[0] += newSides[1] * newSides[2];
                            sxy[1] += newSides[0] * newSides[2];
                            sxy[2] += newSides[0] * newSides[1];
                        }
                        for(int i = top; i >= top + gtb; i--) {
                            for(int side = 0; side < 3; side++) {
                                newSides[side] = stack[i][side];
                            }
                            tempWeight += ((allsides[1] - newSides[1]) * (allsides[2] - newSides[2]) - sxy[0] + newSides[1] * newSides[2]) * newSides[0] * (newSides[0] - 1) +
                                            ((allsides[0] - newSides[0]) * (allsides[2] - newSides[2]) - sxy[1] + newSides[0] * newSides[2]) * newSides[1] * (newSides[1] - 1) +
                                            ((allsides[0] - newSides[0]) * (allsides[1] - newSides[1]) - sxy[2] + newSides[0] * newSides[1]) * newSides[2] * (newSides[2] - 1);
                        }

                        for(int side = 0; side < 3; side++) {
                            stack[top + gtb][side] = allsides[side] - stack[top][side];
                        }
                        weights[ii] += tempWeight;
                        top = top + gtb + 1;
		
					} // End of polytomy section
		
				}
				
			}
			for(int i = 0; i < numRuns; i++) {
				queue2Helper.add(new ComparablePair<Long, Integer>(weights[i], positions[i]));
			}
			synchronized(positionOutLock) {
				while(!queue2Helper.isEmpty() && positionOut == queue2Helper.peek().value.intValue() && queue2.add(queue2Helper.remove().key)) {
					positionOut++;
					if(positionOut % 100000 == 0) {
						System.out.println(positionOut + " weights calculated in turntasktoscores " + ((double)System.nanoTime() - timer3)/1000000000);
						timer3 = System.nanoTime();
					}
				}
			}

			threadCount.decrementAndGet();

		}

	}
}
