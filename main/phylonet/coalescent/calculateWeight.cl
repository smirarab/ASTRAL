__constant sampler_t sampler =
      CLK_NORMALIZED_COORDS_FALSE
    | CLK_ADDRESS_CLAMP_TO_EDGE
    | CLK_FILTER_NEAREST;
inline long F(long a, long b, long c) {
	return ((a + b + c - 3))*a*b*c;
}
inline long FF(long a1, long a2, long a3, long b1, long b2, long b3, long c1, long c2, long c3) {
	return a1*((a1 + b2 + c3 - 3)*b2*c3 + (a1 + b3 + c2 - 3)*b3*c2) + a2*((a2 + b1 + c3 - 3)*b1*c3 + (a2 + b3 + c1 - 3)*b3*c1) + a3*((a3 + b1 + c2 - 3)*b1*c2 + (a3 + b2 + c1 - 3)*b2*c1);
}
struct cl_Tripartition {
	__global const long * cluster1;
	__global const long * cluster2;
	__global const long * cluster3;
};
inline uint popcnt(const ulong i) {
        uint n;
        asm("popc.b64 %0, %1;" : "=r"(n) : "l" (i));
        return n;
}
inline int bitIntersectionSize(__global const long input1[SPECIES_WORD_LENGTH], __global const long input2[SPECIES_WORD_LENGTH]) {
	int out = 0;
	for (int i = 0; i < SPECIES_WORD_LENGTH; i++) {
		out += popcnt(input1[i]&input2[i * WORK_GROUP_SIZE]);
	}
	return out;
}
__kernel void calcWeight(
	__global const short* geneTreesAsInts,
	int geneTreesAsIntsLength,
	__global const long* allArray,
	__global const long* tripartitions1glob,
	__global const long* tripartitions2glob,
	__global const long* tripartitions3glob,
	__global long* weightArray
){
	long weight = 0;
	struct cl_Tripartition trip;
	int idx = get_global_id(0);
	
	int allsides[3];

	int newTree = 1;
	int counter = 0;
	int treeCounter = 0;

	long overlap [(TAXON_SIZE + 1) * 3];
	long overlapind [(TAXON_SIZE + 1) * 3];
	long stack[STACK_SIZE * 3];

	int top = 0;
	short geneInt = 0;
	while(counter < geneTreesAsIntsLength){
		geneInt = geneTreesAsInts[counter];
		counter++;
		
		if (newTree) {
			newTree = 0;

			allsides[0] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], tripartitions1glob + idx);
			allsides[1] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], tripartitions2glob + idx);
			allsides[2] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], tripartitions3glob + idx);

			treeCounter++;
		}
		if (geneInt >= 0) {
			stack[top * 3] = ((tripartitions1glob[idx + (SPECIES_WORD_LENGTH - 1 - geneInt / LONG_BIT_LENGTH) * WORK_GROUP_SIZE])>>(geneInt % LONG_BIT_LENGTH)) & 1;
			stack[top * 3 + 1] = ((tripartitions2glob[idx + (SPECIES_WORD_LENGTH - 1 - geneInt / LONG_BIT_LENGTH) * WORK_GROUP_SIZE])>>(geneInt % LONG_BIT_LENGTH)) & 1;
			stack[top * 3 + 2] = ((tripartitions3glob[idx + (SPECIES_WORD_LENGTH - 1 - geneInt / LONG_BIT_LENGTH) * WORK_GROUP_SIZE])>>(geneInt % LONG_BIT_LENGTH)) & 1;
			top++;
		}
		else if (geneInt == INT_MIN) {
			top = 0;
			newTree = 1;
		}
		else if (geneInt == -2) {
			top--;
			int topminus1 = top - 1;
			long topa[3];
			long topminus1a[3];
			
			topa[0] = stack[top * 3];
			topa[1] = stack[top * 3 + 1];
			topa[2] = stack[top * 3 + 2];

			topminus1a[0] = stack[topminus1 * 3];
			topminus1a[1] = stack[topminus1 * 3 + 1];
			topminus1a[2] = stack[topminus1 * 3 + 2];


			long newSides0 = topa[0] + topminus1a[0];
			long newSides1 = topa[1] + topminus1a[1];
			long newSides2 = topa[2] + topminus1a[2];
			long side3s0 = allsides[0] - newSides0;
			long side3s1 = allsides[1] - newSides1;
			long side3s2 = allsides[2] - newSides2;

			weight += FF(topa[0], topa[1], topa[2], topminus1a[0], topminus1a[1], topminus1a[2], side3s0, side3s1, side3s2);	

			stack[topminus1 * 3] = newSides0;
			stack[topminus1 * 3 + 1] = newSides1;
			stack[topminus1 * 3 + 2] = newSides2;

		}
		else { //for polytomies
		
			int nzc[3];
			nzc[0] = nzc[1] = nzc[2] = 0;
			int newSides[3];
			newSides[0] = newSides[1] = newSides[2] = 0;
			
			for(int side = 0; side < 3; side++) {
				for(int i = top - 1; i >= top + geneInt; i--) {
					if(stack[i * 3 + side] > 0) {
						newSides[side] += stack[i * 3 + side];
						overlap[nzc[side]+ side * (TAXON_SIZE + 1)] = stack[i * 3 + side];
						overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = i;
						nzc[side]++;
					}
				}
				
				stack[top * 3 + side] = allsides[side] - newSides[side];
				
				if(stack[top * 3 + side] > 0) {
					overlap[nzc[side] + side * (TAXON_SIZE + 1)] = stack[top * 3 + side];
					overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = top;
					nzc[side]++;					
				}
				stack[(top + geneInt) * 3 + side] = newSides[side];
			}
			
			for(int i = nzc[0] - 1; i >= 0; i--) {
				for(int j = nzc[1] - 1; j >= 0; j--) {
					for(int k = nzc[2] - 1; k >= 0; k--) {
						if(overlapind[i] != overlapind[j + (TAXON_SIZE + 1)] && overlapind[i] != overlapind[k + (TAXON_SIZE + 1) * 2] && overlapind[j + (TAXON_SIZE + 1)] != overlapind[k + (TAXON_SIZE + 1) * 2])
							weight += F(overlap[i], overlap[j + (TAXON_SIZE + 1)], overlap[k + (TAXON_SIZE + 1) * 2]);
					}
				}
			}
			
			top = top + geneInt + 1;
		}
	}

	weightArray[idx] = weight;
}
