void addScale(
	long* in,
	int augend
) 
{
	long topIn = in[0] & 0xffffffff00000000L;
	long botIn = in[0] & 0x00000000ffffffffL;
	botIn = (in[0] + augend) & 0x00000000ffffffffL;
	topIn += botIn;
	in[0] = topIn;

}
void truncate( //right is exclusive. Truncates num1 into a base uint number
	long* in,
	int inLength
)
{


	for (int i = inLength - 1; i > 1; i--) {
		in[i - 1] += (in[i] >> 32) & 0x00000000ffffffffL;
		in[i] = in[i] & 0x00000000ffffffffL;
	}


	if (in[1] >> 32 != 0) {
		
		for(int j = inLength - 1; j > 0; j--) {
			in[j] = in[j-1];
		}
		
		in[1] = (in[2] >> 32) & 0x00000000ffffffffL;
		in[2] = in[2] & UNSIGNEDLONG_BOTTOM_MASK;

		in[0] = (in[0] + 1) & 0xff000000ffffffffL;

	}
	
	else {
	
	
		int counter = 1;
		while(in[counter] == 0) {
			counter++;
		}
		if(counter != 1) {

			for(int j = 1; j < inLength - counter + 1; j++) {
				in[j] = in[j + counter - 1];
			}
			for(int i = inLength - counter + 1; i < inLength; i++) {
				in[i] = 0;
			}
			addScale(in, -(counter - 1));

		}

		
	}

}

void addGlobal(
	 __global long* num1,
	 long* num2,
	 long* out
)
{
	
	long carry = 0;
	long sum;
	int outScale;
	long outSign;
	int scaleDifference = (int) (num1[0] & SCALE_MASK) - (int) (num2[0] & SCALE_MASK);
	if ((num1[0] & ZERO_MASK) != 0 || scaleDifference < -ELEMENT_LENGTH) {
		for (int i = 0; i < ELEMENT_LENGTH; i++) {
			out[i] = num2[i];
		}
		return;
	}
	if ((num2[0] & ZERO_MASK) != 0 || scaleDifference > ELEMENT_LENGTH) {
		for (int i = 0; i < ELEMENT_LENGTH; i++) {
			out[i] = num1[i];
		}
		return;
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
		for (int i = 1; i < scaleDifference + 1 && i < ELEMENT_LENGTH; i++) {
			out[i] = num1[i];
		}
		// adds
		if (subtract == 0) {
			carry = 0;
			// actual adding
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference; i--) {
				out[i] = (num1[i] & UNSIGNEDLONG_BOTTOM_MASK) + (num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK);
//					out[i] = sum & UINT_MAX_VALUE;
//					carry = sum >> 32;
				}
			// if theres a carry, shift all the uints and copy carry to the
			// first one
		}
			// subtracts
		else {
			// actual subtraction
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference + 1; i--) {
				if ((num1[i] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK)) {
					out[i] = (num1[i] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
					carry = 0;
					} else {
					out[i] = UINT_MAX_VALUE - ((num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[i] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1; // +1 because UINT_MAX_VALUE is analogous to 9 in decimal
					carry = -1;
				}
				
			}
			// if num1 and num2 are same size and num2 is bigger, flip sign
			// since it was assigned num1's sign
			if (scaleDifference == 0) {
				if ((num1[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num2[1] & UNSIGNEDLONG_BOTTOM_MASK)) {
					out[1] = (num1[1] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
				} else {
					out[1] = UINT_MAX_VALUE - (num2[1] - num1[1]);
					for (int i = 1; i < ELEMENT_LENGTH; i++) {
						out[i] = out[i] ^ FLIP_MASK;
					}
					outSign = outSign ^ FLIP_FIRST_MASK;
				}
			}
			// subtracts uppermost digit
			else {
				// just need to subtract uppermost digit, no hassle
				if (num1[scaleDifference + 1] + carry >= num2[1]) {
					out[scaleDifference + 1] = (num1[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
				}
					// uppermost digit of num1 is less, have to borrow
				else {
					out[scaleDifference + 1] = UINT_MAX_VALUE - ((num2[1] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1;
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
		for (int i = 1; i < scaleDifference + 1 && i < ELEMENT_LENGTH; i++) {
			out[i] = num2[i];
		}
		// adds
		if (subtract == 0) {
			carry = 0;
			// sum the last few terms of num2 and num1 - scaleDifference
				
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference; i--) {
			
				out[i] = (num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) + (num2[i] & UNSIGNEDLONG_BOTTOM_MASK);
			}

		}
		// subtracts
		else {
			// actual subtraction
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference + 1; i--) {
				if ((num2[i] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK)) {
					out[i] = (num2[i] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
					carry = 0;
				} 
				else {
					out[i] = UINT_MAX_VALUE - ((num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[i] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1; // +1 because UINT_MAX_VALUE is analogous to 9 in decimal
					carry = -1;
				}
			}
			
			// just need to subtract uppermost digit, no hassle
			if ((num2[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num1[1] & UNSIGNEDLONG_BOTTOM_MASK)) {
				out[scaleDifference + 1] = (num2[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
			}
			// uppermost digit of num1 is less, have to borrow
			else {
				out[scaleDifference + 1] = UINT_MAX_VALUE - ((num1[1] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1;
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
	for (int i = 1; i < ELEMENT_LENGTH; i++) {
		if (out[i] != 0) {
			zero = 0;
			break;
		}
	}
	if (zero == 1) {
		out[0] |= SET_ZERO_MASK;
	}
	else {
		out[0] = outSign;
		out[0] |= outScale & 0x00000000ffffffffL;	
	
		truncate(out, ELEMENT_LENGTH);
	}
}
void add(
	 long* num1,
	 long* num2,
	 long* out
)
{
	long carry = 0;
	long sum;
	int outScale;
	long outSign;
	int scaleDifference = (int) (num1[0] & SCALE_MASK) - (int) (num2[0] & SCALE_MASK);
	if ((num1[0] & ZERO_MASK) != 0 || scaleDifference < -ELEMENT_LENGTH) {
		for (int i = 0; i < ELEMENT_LENGTH; i++) {
			out[i] = num2[i];
		}
		return;
	}
	if ((num2[0] & ZERO_MASK) != 0 || scaleDifference > ELEMENT_LENGTH) {
		for (int i = 0; i < ELEMENT_LENGTH; i++) {
			out[i] = num1[i];
		}
		return;
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

		for (int i = 1; i < scaleDifference + 1 && i < ELEMENT_LENGTH; i++) {
			out[i] = num1[i];
		}
		// adds
		if (subtract == 0) {
			carry = 0;
			// actual adding
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference; i--) {
				out[i] = (num1[i] & UNSIGNEDLONG_BOTTOM_MASK) + (num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK);

			}
		}
			// subtracts
		else {
			// actual subtraction
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference + 1; i--) {
				if ((num1[i] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK)) {
					out[i] = (num1[i] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
					carry = 0;
					} else {
					out[i] = UINT_MAX_VALUE - ((num2[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[i] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1; // +1 because UINT_MAX_VALUE is analogous to 9 in decimal
					carry = -1;
				}
				
			}
			// if num1 and num2 are same size and num2 is bigger, flip sign
			// since it was assigned num1's sign
			if (scaleDifference == 0) {
				if ((num1[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num2[1] & UNSIGNEDLONG_BOTTOM_MASK)) {
					out[1] = (num1[1] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
				} else {
					out[1] = UINT_MAX_VALUE - (num2[1] - num1[1]);
					for (int i = 1; i < ELEMENT_LENGTH; i++) {
						out[i] = out[i] ^ FLIP_MASK;
					}
					outSign = outSign ^ FLIP_FIRST_MASK;
				}
			}
			// subtracts uppermost digit
			else {
				// just need to subtract uppermost digit, no hassle
				if (num1[scaleDifference + 1] + carry >= num2[1]) {
					out[scaleDifference + 1] = (num1[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
				}
					// uppermost digit of num1 is less, have to borrow
				else {
					out[scaleDifference + 1] = UINT_MAX_VALUE - ((num2[1] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1;
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
		for (int i = 1; i < scaleDifference + 1 && i < ELEMENT_LENGTH; i++) {
			out[i] = num2[i];
		}
		// adds
		if (subtract == 0) {
			carry = 0;
			// sum the last few terms of num2 and num1 - scaleDifference
				
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference; i--) {
				out[i] = (num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) + (num2[i] & UNSIGNEDLONG_BOTTOM_MASK);
			}

		}
		// subtracts
		else {
			// actual subtraction
			for (int i = ELEMENT_LENGTH - 1; i > scaleDifference + 1; i--) {
				if ((num2[i] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK)) {
					out[i] = (num2[i] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
					carry = 0;
				} 
				else {
					out[i] = UINT_MAX_VALUE - ((num1[i - scaleDifference] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[i] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1; // +1 because UINT_MAX_VALUE is analogous to 9 in decimal
					carry = -1;
				}
			}
			
			// just need to subtract uppermost digit, no hassle
			if ((num2[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK) + carry >= (num1[1] & UNSIGNEDLONG_BOTTOM_MASK)) {
				out[scaleDifference + 1] = (num2[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK) - (num1[1] & UNSIGNEDLONG_BOTTOM_MASK) + carry;
			}
			// uppermost digit of num1 is less, have to borrow
			else {
				out[scaleDifference + 1] = UINT_MAX_VALUE - ((num1[1] & UNSIGNEDLONG_BOTTOM_MASK) - (num2[scaleDifference + 1] & UNSIGNEDLONG_BOTTOM_MASK)) + carry + 1;
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
	for (int i = 1; i < ELEMENT_LENGTH; i++) {
		if (out[i] != 0) {
			zero = 0;
			break;
		}
	}
	if (zero == 1) {
		out[0] |= SET_ZERO_MASK;
	}
	else {
		out[0] = outSign;
		out[0] |= outScale & 0x00000000ffffffffL;	
		
		truncate(out, ELEMENT_LENGTH);
	
	}
	
}

void multiplyLongGlobal( //num1, num2, and out all have same length and ELEMENT_LENGTH. They also all end with a single integer which is divided into 2 bit sign and 6 bit exponent SSEEEEEE 
	__global long* in,
	long multiplicand,
	long* out
)
{
	out[0] = in[0];
	if(multiplicand == 0 || in[0] & ZERO_MASK) {
		out[0] |= 0x4000000000000000L;
		return;
	}

	if(multiplicand < 0) {
		multiplicand = -multiplicand;
		out[0] ^= 0x8000000000000000L;
	}
	
	for (int i = ELEMENT_LENGTH - 1; i > 0; i--) {
		out[i] = (in[i] & UINT_MAX_VALUE) * multiplicand;
	}

	truncate(out, ELEMENT_LENGTH);

}
void multiplyLong( //num1, num2, and out all have same length and ELEMENT_LENGTH. They also all end with a single integer which is divided into 2 bit sign and 6 bit exponent SSEEEEEE 
	long* in,
	long multiplicand,
	long* out
)
{
	out[0] = in[0];
	if(multiplicand == 0 || in[0] & ZERO_MASK) {
		out[0] |= 0x4000000000000000L;
		return;
	}
	if(multiplicand < 0) {
		multiplicand = -multiplicand;
		out[0] ^= 0x8000000000000000L;
	}
	for (int i = ELEMENT_LENGTH - 1; i > 0; i--) {
		out[i] = (in[i] & UINT_MAX_VALUE) * multiplicand;
	}
	
	truncate(out, ELEMENT_LENGTH);

}
void multiplyArray(
	long * num1,
	long * num2,
	long * out
)
{
	ulong product;
	int outScale = (num1[0] & SCALE_MASK) + (num2[0] & SCALE_MASK);
	//to make sure the scale doesn't overflow negatively
	if(outScale > 50000) {
		out[0] = out[0] & SET_ZERO_MASK;
		return;
	}
	long outSign = (num1[0] & SIGN_MASK) ^ (num2[0] & SIGN_MASK);
	long outPrep[(ELEMENT_LENGTH - 1) * 2];
	
	if(num1[0] & ZERO_MASK || num2[0] & ZERO_MASK) {
		out[0] = out[0] | SET_ZERO_MASK;
		return;
	}

	for(int i = 0; i < ELEMENT_LENGTH * 2 - 2; i++) {
		outPrep[i] = 0;
	}
	//actual multiplication
	for(int i = ELEMENT_LENGTH - 1; i > 0; i--) {
		for(int j = ELEMENT_LENGTH - 1; j > 0; j--) {
			product = (ulong)num1[i] * num2[j];
			outPrep[i + j - 2] += ((product & UNSIGNEDLONG_TOP_MASK) >> 32) & UNSIGNEDLONG_BOTTOM_MASK;
			outPrep[i + j - 1] += product & UNSIGNEDLONG_BOTTOM_MASK;

		}
	}

	//turns outPrep into out by truncating and then setting bits
	//truncate(outPrep, ELEMENT_LENGTH * 2 - 1);
	for (int i = ELEMENT_LENGTH * 2 - 3; i > 0; i--) {
		outPrep[i - 1] += (outPrep[i] >> 32) & 0x00000000ffffffffL;
		outPrep[i] = outPrep[i] & 0x00000000ffffffffL;
	}
	//calculate how much outScale should be to make sure it doesnt have leading zeros
	int counter = 0;
	while(!outPrep[counter] && counter < ELEMENT_LENGTH * 5) {
		counter++;
		outScale--;
	}
	//copying the good part from outPrep
	for(int i = 1; i < ELEMENT_LENGTH && counter + i - 1 < (ELEMENT_LENGTH - 1) * 2; i++) {

		out[i] = outPrep[counter + i - 1];
	}
	outScale++;
	out[0] = outSign;
	out[0] |= (outScale & UNSIGNEDLONG_BOTTOM_MASK);

}


bool isMagnitudeLessThanTwo(long * r, long * i) { //r^2 + i^2

	long test[ELEMENT_LENGTH];

	add(r, i, test);

	if((int)(test[0] & UNSIGNEDLONG_BOTTOM_MASK) > 0){
		return false;
	}
	if((int)(test[0] & UNSIGNEDLONG_BOTTOM_MASK) == 0 && (int)test[1] >= 4) {
		return false;
	}
	return true;

}


int calcMandelbrot(long * r, long * i, int maxIteration) {
	long rCalc[ELEMENT_LENGTH];
	long iCalc[ELEMENT_LENGTH];
	long rr[ELEMENT_LENGTH];
	long ii[ELEMENT_LENGTH];
	long ri[ELEMENT_LENGTH];
	long rCalcTemp[ELEMENT_LENGTH];
	long iCalcTemp[ELEMENT_LENGTH];
	int counter = 0;
	
	multiplyArray(r, r, rr);
	multiplyArray(i, i, ii);
	
	
	if(isMagnitudeLessThanTwo(rr, ii)) {
		multiplyArray(r, i, ri);
		
		ii[0] ^= SIGN_MASK;
		add(rr, ii, rCalc);
		multiplyLong(ri, 2, iCalc);

		add(rCalc, r, rCalcTemp);
		add(iCalc, i, iCalcTemp);
		
		counter++;

		multiplyArray(rCalcTemp, rCalcTemp, rr);
		multiplyArray(iCalcTemp, iCalcTemp, ii);
		while(isMagnitudeLessThanTwo(rr, ii) && counter < maxIteration) { //mandelbrot algorithm

			multiplyArray(rCalcTemp, iCalcTemp, ri);
			ii[0] ^= SIGN_MASK;
			add(rr, ii, rCalc);
			multiplyLong(ri, 2, iCalc);
			
			add(rCalc, r, rCalcTemp);
			add(iCalc, i, iCalcTemp);

			multiplyArray(rCalcTemp, rCalcTemp, rr);
			multiplyArray(iCalcTemp, iCalcTemp, ii);
			
			counter++;
		}
	}
	return counter;
}

__kernel void test(__global long * num1, __global long * num2, __global long * num3) {
	long test1[ELEMENT_LENGTH];
	long test2[ELEMENT_LENGTH];
	long test3[ELEMENT_LENGTH];
	for(int i = 0; i < ELEMENT_LENGTH; i++) {
		test1[i] = num1[i];
		test2[i] = num2[i];
	}
	multiplyArray(test1, test1, test3);
	for(int i = 0; i < ELEMENT_LENGTH; i++) {
		num3[i] = test3[i];
	}
}	


__kernel void fillDankness(
	__global int* out,
	int preferredLength,
	int preferredWidth,
	__global long* centerR,
	__global long* centerI,
	__global long* zoom, //going to be zoom * (idx - preferred width) = xshift
	int maxIteration
)
{

	int idx = get_global_id(0);
	int idy = get_global_id(1);
	
	long r[ELEMENT_LENGTH];
	long i[ELEMENT_LENGTH];
	long rTemp[ELEMENT_LENGTH];
	long iTemp[ELEMENT_LENGTH];
	multiplyLongGlobal(zoom, (idx - preferredWidth / 2), rTemp);
	multiplyLongGlobal(zoom, (idy - preferredLength / 2), iTemp);

	addGlobal(centerR, rTemp, r);
	addGlobal(centerI, iTemp, i);
	
	int num = calcMandelbrot(r, i, maxIteration);

	out[idx * preferredLength + idy] = num;
	
	
}
