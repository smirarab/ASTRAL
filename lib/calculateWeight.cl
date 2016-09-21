long F(int a, int b, int c) {
	return ((long)(a + b + c - 3))*a*b*c;
}
struct Tripartition {
	__global long* cluster1;
	__global long* cluster2;
	__global long* cluster3;
};
inline uint popcnt(const ulong i) {
        uint n;
        asm("popc.b64 %0, %1;" : "=r"(n) : "l" (i));
        return n;
}
int bitIntersectionSize(__global long input1[SPECIES_WORD_LENGTH], __global long input2[SPECIES_WORD_LENGTH]) {
	int out = 0;
	for (int i = 0; i < SPECIES_WORD_LENGTH; i++) {
		out += popcnt(input1[i]&input2[i]);
	}
	return out;
}
__kernel void calcWeight(
	__global int* geneTreesAsInts,
	int geneTreesAsIntsLength,
	__global long* allArray,
	__global long* tripartitions,
	__global long* weightArray
){
	long weight = 0;
	struct Tripartition trip;
	int idx = get_global_id(0);
	trip.cluster1 = SPECIES_WORD_LENGTH * 3 * idx + tripartitions;
	trip.cluster2 = SPECIES_WORD_LENGTH * 3 * idx + SPECIES_WORD_LENGTH + tripartitions;
	trip.cluster3 = SPECIES_WORD_LENGTH * 3 * idx + SPECIES_WORD_LENGTH * 2 + tripartitions;
	
	int allsides[3];

	int newTree = 1;
	int counter = 0;
	int treeCounter = 0;

	int stack[(STACK_SIZE+2)*3];
	int overlap [(STACK_SIZE+1)*3];
	int overlapind [(STACK_SIZE+1)*3];
	
	int top = 0;
	for (; counter < geneTreesAsIntsLength; counter++) {

		if (newTree) {
			newTree = 0;

			allsides[0] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster1);
			allsides[1] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster2);
			allsides[2] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster3);

			treeCounter++;

		}
		if (geneTreesAsInts[counter] >= 0) {
			stack[top] = ((trip.cluster1[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts[counter] / LONG_BIT_LENGTH])>>(geneTreesAsInts[counter] % LONG_BIT_LENGTH)) & 1;
			stack[top + (STACK_SIZE+2)] = ((trip.cluster2[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts[counter] / LONG_BIT_LENGTH])>>(geneTreesAsInts[counter] % LONG_BIT_LENGTH)) & 1;
			stack[top + (STACK_SIZE+2)*2] = ((trip.cluster3[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts[counter] / LONG_BIT_LENGTH])>>(geneTreesAsInts[counter] % LONG_BIT_LENGTH)) & 1;
			top++;
		}
		else if (geneTreesAsInts[counter] == INT_MIN) {
			top = 0;
			newTree = 1;
		}
		else if (geneTreesAsInts[counter] == -2) {
			top--;
	
			int newSides0 = stack[top] + stack[top - 1];
			int newSides1 = stack[top + (STACK_SIZE+2)] + stack[top - 1 + (STACK_SIZE+2)];
			int newSides2 = stack[top + (STACK_SIZE+2)*2] + stack[top - 1 + (STACK_SIZE+2)*2];
			
			int side3s0 = allsides[0] - newSides0;
			int side3s1 = allsides[1] - newSides1;
			int side3s2 = allsides[2] - newSides2;

			weight += 
				F(stack[top], stack[top - 1 + (STACK_SIZE + 2)], side3s2) +
				F(stack[top], stack[top - 1 + (STACK_SIZE + 2)*2], side3s1) +
				F(stack[top + (STACK_SIZE + 2)], stack[top - 1], side3s2) +
				F(stack[top + (STACK_SIZE + 2)], stack[top - 1 + (STACK_SIZE + 2)*2], side3s0) +
				F(stack[top + (STACK_SIZE + 2)*2], stack[top - 1], side3s1) +
				F(stack[top + (STACK_SIZE + 2)*2], stack[top - 1 + (STACK_SIZE + 2)], side3s0);
				
			stack[top - 1] = newSides0;
			stack[top - 1 + (STACK_SIZE + 2)] = newSides1;
			stack[top - 1 + (STACK_SIZE + 2)*2] = newSides2;
			
		}
		else { //for polytomies
		
			int nzc[3];
			nzc[0] = nzc[1] = nzc[2] = 0;
			int newSides[3];
			newSides[0] = newSides[1] = newSides[2] = 0;
			
			for(int side = 0; side < 3; side++) {
				for(int i = top - 1; i >= top + geneTreesAsInts[counter]; i--) {
					if(stack[i + side * (STACK_SIZE + 2)] > 0) {
						newSides[side] += stack[i + side * (STACK_SIZE + 2)];
						overlap[nzc[side]+ side * (STACK_SIZE + 1)] = stack[i + side * (STACK_SIZE + 2)];
						overlapind[nzc[side] + side * (STACK_SIZE + 1)] = i;
						nzc[side]++;
					}
				}
				
				stack[top + side * (STACK_SIZE + 2)] = allsides[side] - newSides[side];
				
				if(stack[top + side * (STACK_SIZE + 2)] > 0) {
					overlap[nzc[side] + side * (STACK_SIZE + 1)] = stack[top + side * (STACK_SIZE + 2)];
					overlapind[nzc[side] + side * (STACK_SIZE + 1)] = top;
					nzc[side]++;					
				}
				stack[top + geneTreesAsInts[counter] + side * (STACK_SIZE + 2)] = newSides[side];
			}
			
			for(int i = nzc[0] - 1; i >= 0; i--) {
				for(int j = nzc[1] - 1; j >= 0; j--) {
					for(int k = nzc[2] - 1; k >= 0; k--) {
						if(overlapind[i] != overlapind[j + (STACK_SIZE + 1)] && overlapind[i] != overlapind[k + (STACK_SIZE + 1) * 2] && overlapind[j + (STACK_SIZE + 1)] != overlapind[k + (STACK_SIZE + 1) * 2])
							weight += F(overlap[i], overlap[j + (STACK_SIZE + 1)], overlap[k + (STACK_SIZE + 1) * 2]);
					}
				}
			}
			
			top = top + geneTreesAsInts[counter] + 1;
			
		}
	}
	weightArray[idx] = weight;
}
