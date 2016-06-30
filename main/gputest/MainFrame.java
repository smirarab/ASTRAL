
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
import static org.jocl.CL.CL_CONTEXT_PLATFORM;

import java.awt.Dimension;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

import javax.swing.JFrame;

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

public class MainFrame extends JFrame {

	public static final int PREFERRED_WIDTH = 256;
	public static final int PREFERRED_LENGTH = 256;
	public static final int PRECISION = 3;
	public static final int MAX_ITERATIONS = 255;
	public static final long SIGN_MASK = 0x8000000000000000L;
	public static final long UINT_MAX_VALUE = 0xffffffffL;
	public static final long SCALE_MASK = 0xffffffffL;
	public static final long ZERO_MASK = 0x4000000000000000L;
	public static final long SET_ZERO_MASK = 0x4000000000000000L;
	public static final long X_SCROLL_SPEED = 10;
	public static final long Y_SCROLL_SPEED = 10;
	public static final long ZOOMIN_SPEED = 3435974000L;
	public static final double ZOOMOUT_SPEED = 1.2;
	public static final long FLIP_FIRST_MASK = 0x8000000000000000L;
	public static final long FLIP_MASK = 0x00000000ffffffffL;
	public static final int WORK_GROUP_SIZE = 64;
	
	public static final int check1 = 62;
	public static final int check2 = 0;

	public long[] zoom = new long[PRECISION];
	public long[] centerR = new long[PRECISION];
	public long[] centerI = new long[PRECISION];
	public int[] dankness = new int[PREFERRED_WIDTH * PREFERRED_LENGTH];
	public MainPanel MP = new MainPanel(dankness, PREFERRED_WIDTH, PREFERRED_LENGTH, MAX_ITERATIONS);

	private cl_context context;
	private cl_command_queue commandQueue;
	private cl_kernel kernel;
	private cl_mem d_dankness;
	private cl_mem d_centerR;
	private cl_mem d_centerI;
	private cl_mem d_zoom;

	public MainFrame() {
		zoom[1] = 3435974000L/64;
		zoom[0] = 0xffffffffL;
		
		centerR[0] = ZERO_MASK;
		centerI[0] = ZERO_MASK;
		
		add(MP);

		initCL();

		setSize(new Dimension(PREFERRED_WIDTH, PREFERRED_LENGTH));
		setVisible(true);
		addKeyListener(new KeyPressed());

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

	}

	public void normalize(long[] in) {

		// converting (8 16) to 96 basically
		for (int i = PRECISION - 1; i > 1; i--) {
			in[i - 1] += in[i] >>> 32;
			in[i] = in[i] & 0x00000000ffffffffL;
		}
		// in case in[1] is over uintmaxvalue, then shift everything over
		if (in[1] >>> 32 != 0) {
			System.arraycopy(in, 1, in, 2, PRECISION - 2);
			in[1] = in[2] >>> 32;
			in[2] = in[2] & UINT_MAX_VALUE;
			in[0]++;
			in[0] = in[0] & 0xff000000ffffffffL;
		}
		// get rid of leading 0s
		else {
			int counter = 1;
			while (in[counter] == 0) {
				counter++;
			}
			if (counter != 1) {
				System.arraycopy(in, counter, in, 1, PRECISION - counter);
			}
			for (int i = PRECISION - counter + 1; i < PRECISION; i++) {
				in[i] = 0;
			}
			addScale(in, -(counter - 1));
		}

	}

	public void addScale(long[] in, int augend) {
		long topIn = in[0] & 0xffffffff00000000L;
		long botIn = in[0] & 0x00000000ffffffffL;
		botIn = (in[0] + augend) & 0x00000000ffffffffL;
		topIn += botIn;
		in[0] = topIn;
	}

	public long[] multiply(long[] in, long multiplicand) {

		long[] out = new long[PRECISION];
		out[0] = in[0];
		if (multiplicand == 0) {
			out[0] |= 0x4000000000000000L;
		}
		if (multiplicand < 0) {
			multiplicand = -multiplicand;
			out[0] ^= 0x8000000000000000L;
		}
		for (int i = PRECISION - 1; i > 0; i--) {
			out[i] = (in[i] & UINT_MAX_VALUE) * multiplicand;
		}
		normalize(out);
		return out;
	}

	public long[] multiply(long[] in, double multiplicand) {
		long[] out = new long[PRECISION];
		out[0] = in[0];

		if (multiplicand == 0) {
			out[0] |= 0x4000000000000000L;
		}
		if (multiplicand < 0) {
			multiplicand = -multiplicand;
			out[0] ^= 0x8000000000000000L;
		}
		for (int i = PRECISION - 1; i > 0; i--) {
			out[i] = (long) ((in[i] & UINT_MAX_VALUE) * multiplicand);
		}
		normalize(out);
		return out;
	}

	private void zoomIn() {
		zoom = multiply(zoom, ZOOMIN_SPEED);
		addScale(zoom, -1);
	}

	private void zoomOut() {
		zoom = multiply(zoom, ZOOMOUT_SPEED);
	}

	public long[] add(long[] num1, long[] num2) {
		long carry = 0;
		long sum;
		int outScale;
		long outSign;
		long[] out = new long[PRECISION];
		int scaleDifference = (int) (num1[0] & SCALE_MASK) - (int) (num2[0] & SCALE_MASK);
		if ((num1[0] & ZERO_MASK) != 0) {
			for (int i = 0; i < PRECISION; i++) {
				out[i] = num2[i];
			}
			return out;
		}
		if ((num2[0] & ZERO_MASK) != 0) {
			for (int i = 0; i < PRECISION; i++) {
				out[i] = num1[i];
			}
			return out;
		}

		// check if we have to subtract
		int subtract = 0;
		if (((num1[0] & SIGN_MASK) ^ (num2[0] & SIGN_MASK)) != 0) {
			subtract = 1;
		}
		if (scaleDifference >= 0) {
			outScale = (int) (num1[0]);
			outSign = num1[0] & SIGN_MASK;
			// copies first few digits of num1
			for (int i = 1; i < scaleDifference + 1 && i < PRECISION; i++) {
				out[i] = num1[i];
			}

			// adds
			if (subtract == 0) {
				carry = 0;

				// actual adding
				for (int i = PRECISION - 1; i > scaleDifference; i--) {
					out[i] = (num1[i] & UINT_MAX_VALUE) + (num2[i - scaleDifference] & UINT_MAX_VALUE);
					// out[i] = sum & UINT_MAX_VALUE;
					// carry = sum >>> 32;

				}

				// if theres a carry, shift all the uints and copy carry to the
				// first one

			}

			// subtracts
			else {

				// actual subtraction
				for (int i = PRECISION - 1; i > scaleDifference + 1; i--) {

					if ((num1[i] & UINT_MAX_VALUE) + carry >= (num2[i - scaleDifference] & UINT_MAX_VALUE)) {
						out[i] = (num1[i] & UINT_MAX_VALUE) - (num2[i - scaleDifference] & UINT_MAX_VALUE) + carry;
						carry = 0;

					} else {
						out[i] = UINT_MAX_VALUE
								- ((num2[i - scaleDifference] & UINT_MAX_VALUE) - (num1[i] & UINT_MAX_VALUE)) + carry
								+ 1; // +1
						// because
						// UINT_MAX_VALUE
						// is
						// analogous
						// to
						// 9
						// in
						// decimal
						carry = -1;
					}
				}

				// if num1 and num2 are same size and num2 is bigger, flip sign
				// since it was assigned num1's sign
				if (scaleDifference == 0) {
					if ((num1[1] & UINT_MAX_VALUE) + carry >= (num2[1] & UINT_MAX_VALUE)) {
						out[1] = (num1[1] & UINT_MAX_VALUE) - (num2[1] & UINT_MAX_VALUE) + carry;
					} else {
						out[1] = UINT_MAX_VALUE - (num2[1] - num1[1]);
						for (int i = 1; i < PRECISION; i++) {
							out[i] = out[i] ^ FLIP_MASK;
						}
						outSign = outSign ^ FLIP_FIRST_MASK;
					}
				}

				// subtracts uppermost digit
				else {
					// just need to subtract uppermost digit, no hassle
					if (num1[scaleDifference + 1] + carry >= num2[1]) {
						out[scaleDifference + 1] = (num1[scaleDifference + 1] & UINT_MAX_VALUE)
								- (num2[1] & UINT_MAX_VALUE) + carry;
					}

					// uppermost digit of num1 is less, have to borrow
					else {
						out[scaleDifference + 1] = UINT_MAX_VALUE
								- ((num2[1] & UINT_MAX_VALUE) - (num1[scaleDifference + 1] & UINT_MAX_VALUE)) + carry
								+ 1;
						int counter = scaleDifference;

						// does borrowing
						while (out[counter] == 0 && counter > 0) {
							out[counter] = UINT_MAX_VALUE;
							counter--;
						}
						out[counter]--;
					}
				}
			}
		}

		// num2 scale is bigger for sure
		else {
			scaleDifference = -scaleDifference;
			outScale = (int) num2[0];
			outSign = num2[0] & SIGN_MASK;

			// copies first few digits of num2
			for (int i = 1; i < scaleDifference + 1 && i < PRECISION; i++) {
				out[i] = num2[i];
			}

			// adds
			if (subtract == 0) {
				carry = 0;

				// sum the last few terms of num2 and num1 - scaleDifference
				for (int i = PRECISION - 1; i > scaleDifference; i--) {
					out[i] = (num1[i - scaleDifference] & UINT_MAX_VALUE) + (num2[i] & UINT_MAX_VALUE);
				}

				// sum the carry and keeps carrying in case it is needed
				for (int i = scaleDifference; carry != 0 && i > 0; i--) {
					sum = (num2[i] & UINT_MAX_VALUE) + carry;
					out[i] = sum & UINT_MAX_VALUE;
					carry = sum >>> 32;
				}
				// if there's an overflow, shift everything over and add in the
				// carry
				if (carry != 0) {

					for (int i = PRECISION - 1; i > 1; i--) {
						out[i] = out[i - 1];
					}
					out[1] = carry;
					outScale++; // increments scale

				}
			}

			// subtracts
			else {
				// actual subtraction
				for (int i = PRECISION - 1; i > scaleDifference + 1; i--) {

					if ((num2[i] & UINT_MAX_VALUE) + carry >= (num1[i - scaleDifference] & UINT_MAX_VALUE)) {
						out[i] = (num2[i] & UINT_MAX_VALUE) - (num1[i - scaleDifference] & UINT_MAX_VALUE) + carry;
						carry = 0;

					} else {
						out[i] = UINT_MAX_VALUE
								- ((num1[i - scaleDifference] & UINT_MAX_VALUE) - (num2[i] & UINT_MAX_VALUE)) + carry
								+ 1; // +1
						// because
						// UINT_MAX_VALUE
						// is
						// analogous
						// to
						// 9
						// in
						// decimal
						carry = -1;
					}
				}

				// just need to subtract uppermost digit, no hassle
				if ((num2[scaleDifference + 1] & UINT_MAX_VALUE) + carry >= (num1[1] & UINT_MAX_VALUE)) {
					out[scaleDifference + 1] = (num2[scaleDifference + 1] & UINT_MAX_VALUE) - (num1[1] & UINT_MAX_VALUE)
							+ carry;
				}

				// uppermost digit of num1 is less, have to borrow
				else {
					out[scaleDifference + 1] = UINT_MAX_VALUE
							- ((num1[1] & UINT_MAX_VALUE) - (num2[scaleDifference + 1] & UINT_MAX_VALUE)) + carry + 1;
					int counter = scaleDifference;

					// does borrowing
					while (out[counter] == 0) {
						out[counter] = UINT_MAX_VALUE;
						counter--;
					}
					out[counter]--;
				}

			}
		}
		int zero = 1;
		for (int i = 1; i < PRECISION; i++) {
			if (out[i] != 0) {
				zero = 0;
				break;
			}
		}

		if (zero == 1) {
			out[0] |= SET_ZERO_MASK;
		} else {
			out[0] |= outSign;
			out[0] += outScale & 0x00000000ffffffffL;
			normalize(out);
		}
		return out;
	}

	public void moveUp() {
		centerI = add(centerI, multiply(zoom, -Y_SCROLL_SPEED));
	}

	public void moveDown() {
		centerI = add(centerI, multiply(zoom, Y_SCROLL_SPEED));
	}

	public void moveRight() {
		centerR = add(centerR, multiply(zoom, X_SCROLL_SPEED));
	}

	public void moveLeft() {
		centerR = add(centerR, multiply(zoom, -X_SCROLL_SPEED));
	}

	private class KeyPressed implements KeyListener {

		@Override
		public void keyPressed(KeyEvent keyEvent) {
			int keyCode = keyEvent.getKeyCode();

			if (keyCode == KeyEvent.VK_LEFT) {
				moveLeft();
				update();
			}
			if (keyCode == KeyEvent.VK_RIGHT) {
				moveRight();
				update();
			}
			if (keyCode == KeyEvent.VK_DOWN) {
				moveDown();
				update();
			}
			if (keyCode == KeyEvent.VK_UP) {
				moveUp();
				update();
			}
			if (keyCode == KeyEvent.VK_1) {
				zoomIn();
				update();
			}
			if (keyCode == KeyEvent.VK_2) {
				zoomOut();
				update();
			}
			/*System.out.print(Long.toHexString(zoom[0]) + " ");
			for (int j = 1; j < zoom.length; j++) {

				System.out.print(((zoom[j] & 0x00000000ffffffffL)) + " ");

			}
			System.out.println();
			System.out.print(Long.toHexString(centerR[0]) + " ");
			for (int j = 1; j < centerR.length; j++) {

				System.out.print(((centerR[j] & 0x00000000ffffffffL)) + " ");

			}
			System.out.println();
			System.out.print(Long.toHexString(centerI[0]) + " ");
			for (int j = 1; j < centerI.length; j++) {

				System.out.print(((centerI[j] & 0x00000000ffffffffL)) + " ");

			}
			System.out.println();
			System.out.println();
			System.out.println();
			*/
		}

		@Override
		public void keyReleased(KeyEvent arg0) {
			// TODO Auto-generated method stub

		}

		@Override
		public void keyTyped(KeyEvent arg0) {
			// TODO Auto-generated method stub

		}

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
		String source = readFile("toDank.cl");
		source = source.replaceAll("ELEMENT_LENGTH", Integer.toString(PRECISION));
		source = source.replaceAll("SCALE_MASK", "0x3fffffffffffffff");
		source = source.replaceAll("SIGN_MASK", "0x8000000000000000");
		source = source.replaceAll("SET_ZERO_MASK", "0x4000000000000000");
		source = source.replaceAll("ZERO_MASK", "0x4000000000000000");
		source = source.replaceAll("FLIP_FIRST_MASK", "0x8000000000000000");
		source = source.replaceAll("UINT_MAX_VALUE", "0x00000000ffffffff");
		source = source.replaceAll("FLIP_MASK", "0x00000000ffffffff");
		source = source.replaceAll("UNSIGNEDLONG_TOP_MASK", "0xffffffff00000000");
		source = source.replaceAll("UNSIGNEDLONG_BOTTOM_MASK", "0x00000000ffffffff");
		source = source.replaceAll("check1", Integer.toString(check1));
		source = source.replaceAll("check2", Integer.toString(check2));
		
		// Create the program
		cl_program cpProgram = clCreateProgramWithSource(context, 1, new String[] { source }, null, null);

		// Build the program
		clBuildProgram(cpProgram, 0, null, "-cl-mad-enable", null, null);

		// Create the kernel
		kernel = clCreateKernel(cpProgram, "fillDankness", null);

		// Create the memory object which will be filled with the
		// pixel data

		d_dankness = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_WRITE_ONLY,
				Sizeof.cl_int * PREFERRED_WIDTH * PREFERRED_LENGTH, Pointer.to(dankness), null);
		d_centerR = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * PRECISION,
				Pointer.to(centerR), null);
		d_centerI = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * PRECISION,
				Pointer.to(centerI), null);
		d_zoom = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY, Sizeof.cl_long * PRECISION,
				Pointer.to(zoom), null);



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

	private void update() {
		
		// Set work size and execute the kernel
		

		long globalWorkSize[] = new long[2];
		globalWorkSize[0] = WORK_GROUP_SIZE;
		globalWorkSize[1] = WORK_GROUP_SIZE;

		clEnqueueWriteBuffer(commandQueue, d_centerR, CL_TRUE, 0L, Sizeof.cl_long * PRECISION, Pointer.to(centerR), 0,
				null, null);
		clEnqueueWriteBuffer(commandQueue, d_centerI, CL_TRUE, 0L, Sizeof.cl_long * PRECISION, Pointer.to(centerI), 0,
				null, null);
		clEnqueueWriteBuffer(commandQueue, d_zoom, CL_TRUE, 0L, Sizeof.cl_long * PRECISION, Pointer.to(zoom), 0, null,
				null);
		
		clSetKernelArg(kernel, 0, Sizeof.cl_mem, Pointer.to(d_dankness));
		clSetKernelArg(kernel, 1, Sizeof.cl_int, Pointer.to(new int[] { PREFERRED_LENGTH }));
		clSetKernelArg(kernel, 2, Sizeof.cl_int, Pointer.to(new int[] { PREFERRED_WIDTH }));
		clSetKernelArg(kernel, 3, Sizeof.cl_mem, Pointer.to(d_centerR));
		clSetKernelArg(kernel, 4, Sizeof.cl_mem, Pointer.to(d_centerI));
		clSetKernelArg(kernel, 5, Sizeof.cl_mem, Pointer.to(d_zoom));
		clSetKernelArg(kernel, 6, Sizeof.cl_int, Pointer.to(new int[]{MAX_ITERATIONS}));
//		System.out.println(Arrays.toString(dankness));
//		System.out.println(Arrays.toString(centerR));
//		System.out.println(Arrays.toString(centerI));
		System.out.println(Arrays.toString(zoom));
		
		for(int i = 0; i < PREFERRED_WIDTH / WORK_GROUP_SIZE; i++) {
			for(int j = 0; j < PREFERRED_LENGTH / WORK_GROUP_SIZE; j++) {
				clEnqueueNDRangeKernel(commandQueue, kernel, 2, new long[]{WORK_GROUP_SIZE * i, WORK_GROUP_SIZE * j}, globalWorkSize, null, 0, null, null);		
			}
		}
		// Read the pixel data into the BufferedImage
        clEnqueueReadBuffer(commandQueue, d_dankness, CL_TRUE, 0L, Sizeof.cl_int * PREFERRED_LENGTH * PREFERRED_WIDTH, Pointer.to(dankness), 0, null, null);

		repaint();
		
	}
}
