__constant sampler_t sampler =
      CLK_NORMALIZED_COORDS_FALSE
    | CLK_ADDRESS_CLAMP_TO_EDGE
    | CLK_FILTER_NEAREST;
inline ulong clock() {
	ulong xa;
	asm("mov.u64 %0, %%globaltimer;" : "=l"(xa));
	return xa;
}
inline long F(ushort a, ushort b, ushort c) {
	return ((long)(a + b + c - 3))*a*b*c;
}
inline long FF(ushort a1, ushort a2, ushort a3, ushort b1, ushort b2, ushort b3, ushort c1, ushort c2, ushort c3) {
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
int bitIntersectionSize(__global const long input1[SPECIES_WORD_LENGTH], __global const long input2[SPECIES_WORD_LENGTH]) {
	int out = 0;
	for (int i = 0; i < SPECIES_WORD_LENGTH; i++) {
		out += popcnt(input1[i]&input2[i]);
	}
	return out;
}
__kernel void calcWeight(
	__read_only image2d_t geneTreesAsInts,
	int geneTreesAsIntsLength,
	__global const long* allArray,
	__global const long* tripartitions1glob,
	__global const long* tripartitions2glob,
	__global const long* tripartitions3glob,
	__global long* weightArray,
	__global ushort* stack,
	__global ulong* profile
){
	long weight = 0;
	struct cl_Tripartition trip;
	int idx = get_global_id(0);
	int globallen = get_global_size(0);
	int globallen2 = globallen * 2;
	int globallen3 = globallen * 3;
	
	trip.cluster1 = tripartitions1glob + idx * SPECIES_WORD_LENGTH;
	trip.cluster2 = tripartitions2glob + idx * SPECIES_WORD_LENGTH;
	trip.cluster3 = tripartitions3glob + idx * SPECIES_WORD_LENGTH;

	int allsides[3];

	int newTree = 1;
	int counter = 0;
	int treeCounter = 0;

	ushort overlap [(TAXON_SIZE + 1) * 3];
	ushort overlapind [(TAXON_SIZE + 1) * 3];

	int top = 0;
	int4 geneTreesAsInts4Ints;
	ulong starttime = clock();
	while(1){
		geneTreesAsInts4Ints = read_imagei(geneTreesAsInts, sampler, (int2)((counter / 4) % IMAGE_WIDTH, counter / 4 / IMAGE_WIDTH));
		/*
		if(0 && idx == 0){
			if(geneTreesAsInts2[counter] != geneTreesAsInts4Ints.x)
				printf("x %d %d %d %d\n", counter, geneTreesAsIntsLength, geneTreesAsInts2[counter], geneTreesAsInts4Ints.x);
			if(geneTreesAsInts2[counter + 1] != geneTreesAsInts4Ints.y)
				printf("y %d %d %d %d\n", counter + 1, geneTreesAsIntsLength, geneTreesAsInts2[counter + 1], geneTreesAsInts4Ints.y);
			if(geneTreesAsInts2[counter + 2] != geneTreesAsInts4Ints.z)
				printf("z %d %d %d %d\n", counter + 2, geneTreesAsIntsLength, geneTreesAsInts2[counter + 2], geneTreesAsInts4Ints.z);
			if(geneTreesAsInts2[counter + 3] != geneTreesAsInts4Ints.w)
				printf("w %d %d %d %d\n", counter + 3, geneTreesAsIntsLength, geneTreesAsInts2[counter + 3], geneTreesAsInts4Ints.w);
		}
		*/
		/*if(top * globallen3 + globallen2 + idx >= globallen3 * (STACK_SIZE + 2))
			printf("SOMETHINGWRONG %d %d %d %d\n", top * globallen3 + globallen2 + idx, globallen3 * (STACK_SIZE + 2), top, globallen);
		*/
		if(counter < geneTreesAsIntsLength) {
			counter++;
			if (newTree) {
				if(idx == 0){
					starttime = clock();
				}
				newTree = 0;

				allsides[0] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster1);
				allsides[1] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster2);
				allsides[2] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster3);

				treeCounter++;
				if(idx == 0){
					profile[0] += clock() - starttime;
				}
			}
			if (geneTreesAsInts4Ints.x >= 0) {
				if(idx == 0)
					starttime = clock();
				stack[top * globallen3 + idx] = ((trip.cluster1[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.x / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.x % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen + idx] = ((trip.cluster2[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.x / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.x % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen2 + idx] = ((trip.cluster3[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.x / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.x % LONG_BIT_LENGTH)) & 1;
				top++;
				if(idx == 0)
					profile[1] += clock() - starttime;
			}
			else if (geneTreesAsInts4Ints.x == INT_MIN) {
				top = 0;
				newTree = 1;
			}
			else if (geneTreesAsInts4Ints.x == -2) {
				if(idx == 0)
					starttime = clock();
				top--;
				int topminus1 = top - 1;
				ushort newSides0 = stack[top * globallen3 + idx] + stack[topminus1 * globallen3 + idx];
				ushort newSides1 = stack[top * globallen3 + globallen + idx] + stack[topminus1 * globallen3 + globallen + idx];
				ushort newSides2 = stack[top * globallen3 + globallen2 + idx] + stack[topminus1 * globallen3 + globallen2 + idx];
				
				ushort side3s0 = allsides[0] - newSides0;
				ushort side3s1 = allsides[1] - newSides1;
				ushort side3s2 = allsides[2] - newSides2;
				/*	
				weight += 
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s2) +
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s1) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + idx], side3s2) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], side3s1) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s0);
				*/

				weight += FF(stack[top * globallen3 + idx], stack[top * globallen3 + globallen + idx], stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0, side3s1, side3s2);	
				stack[topminus1 * globallen3 + idx] = newSides0;
				stack[topminus1 * globallen3 + globallen + idx] = newSides1;
				stack[topminus1 * globallen3 + globallen2 + idx] = newSides2;

				if(idx == 0)
					profile[2] += clock() - starttime;
			}
			else { //for polytomies
			
				if(idx == 0)
					starttime = clock();
				int nzc[3];
				nzc[0] = nzc[1] = nzc[2] = 0;
				int newSides[3];
				newSides[0] = newSides[1] = newSides[2] = 0;
				
				for(int side = 0; side < 3; side++) {
					for(int i = top - 1; i >= top + geneTreesAsInts4Ints.x; i--) {
						if(stack[i * globallen3 + side * globallen + idx] > 0) {
							newSides[side] += stack[i * globallen3 + side * globallen + idx];
							overlap[nzc[side]+ side * (TAXON_SIZE + 1)] = stack[i * globallen3 + side * globallen + idx];
							overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = i;
							nzc[side]++;
						}
					}
					
					stack[top * globallen3 + side * globallen + idx] = allsides[side] - newSides[side];
					
					if(stack[top * globallen3 + side * globallen + idx] > 0) {
						overlap[nzc[side] + side * (TAXON_SIZE + 1)] = stack[top * globallen3 + side * globallen + idx];
						overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = top;
						nzc[side]++;					
					}
					stack[(top + geneTreesAsInts4Ints.x) * globallen3 + side * globallen + idx] = newSides[side];
				}
				
				for(int i = nzc[0] - 1; i >= 0; i--) {
					for(int j = nzc[1] - 1; j >= 0; j--) {
						for(int k = nzc[2] - 1; k >= 0; k--) {
							if(overlapind[i] != overlapind[j + (TAXON_SIZE + 1)] && overlapind[i] != overlapind[k + (TAXON_SIZE + 1) * 2] && overlapind[j + (TAXON_SIZE + 1)] != overlapind[k + (TAXON_SIZE + 1) * 2])
								weight += F(overlap[i], overlap[j + (TAXON_SIZE + 1)], overlap[k + (TAXON_SIZE + 1) * 2]);
						}
					}
				}
				
				top = top + geneTreesAsInts4Ints.x + 1;
				if(idx == 0)
					profile[3] += clock() - starttime;		
			}
		}
		else {
			break;
		}
		if(counter < geneTreesAsIntsLength) {
			counter++;
			if (newTree) {
				newTree = 0;

				allsides[0] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster1);
				allsides[1] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster2);
				allsides[2] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster3);

				treeCounter++;

			}
			if (geneTreesAsInts4Ints.y >= 0) {
				stack[top * globallen3 + idx] = ((trip.cluster1[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.y / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.y % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen + idx] = ((trip.cluster2[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.y / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.y % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen2 + idx] = ((trip.cluster3[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.y / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.y % LONG_BIT_LENGTH)) & 1;
				top++;
			}
			else if (geneTreesAsInts4Ints.y == INT_MIN) {
				top = 0;
				newTree = 1;
			}
			else if (geneTreesAsInts4Ints.y == -2) {
				top--;
				int topminus1 = top - 1;
				int newSides0 = stack[top * globallen3 + idx] + stack[topminus1 * globallen3 + idx];
				int newSides1 = stack[top * globallen3 + globallen + idx] + stack[topminus1 * globallen3 + globallen + idx];
				int newSides2 = stack[top * globallen3 + globallen2 + idx] + stack[topminus1 * globallen3 + globallen2 + idx];
				
				int side3s0 = allsides[0] - newSides0;
				int side3s1 = allsides[1] - newSides1;
				int side3s2 = allsides[2] - newSides2;
				/*
				weight += 
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s2) +
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s1) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + idx], side3s2) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], side3s1) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s0);
				*/	
				weight += FF(stack[top * globallen3 + idx], stack[top * globallen3 + globallen + idx], stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0, side3s1, side3s2);	
				stack[topminus1 * globallen3 + idx] = newSides0;
				stack[topminus1 * globallen3 + globallen + idx] = newSides1;
				stack[topminus1 * globallen3 + globallen2 + idx] = newSides2;
				
			}
			else { //for polytomies
			
				int nzc[3];
				nzc[0] = nzc[1] = nzc[2] = 0;
				int newSides[3];
				newSides[0] = newSides[1] = newSides[2] = 0;
				
				for(int side = 0; side < 3; side++) {
					for(int i = top - 1; i >= top + geneTreesAsInts4Ints.y; i--) {
						if(stack[i * globallen3 + side * globallen + idx] > 0) {
							newSides[side] += stack[i * globallen3 + side * globallen + idx];
							overlap[nzc[side]+ side * (TAXON_SIZE + 1)] = stack[i * globallen3 + side * globallen + idx];
							overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = i;
							nzc[side]++;
						}
					}
					
					stack[top * globallen3 + side * globallen + idx] = allsides[side] - newSides[side];
					
					if(stack[top * globallen3 + side * globallen + idx] > 0) {
						overlap[nzc[side] + side * (TAXON_SIZE + 1)] = stack[top * globallen3 + side * globallen + idx];
						overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = top;
						nzc[side]++;					
					}
					stack[(top + geneTreesAsInts4Ints.y) * globallen3 + side * globallen + idx] = newSides[side];
				}
				
				for(int i = nzc[0] - 1; i >= 0; i--) {
					for(int j = nzc[1] - 1; j >= 0; j--) {
						for(int k = nzc[2] - 1; k >= 0; k--) {
							if(overlapind[i] != overlapind[j + (TAXON_SIZE + 1)] && overlapind[i] != overlapind[k + (TAXON_SIZE + 1) * 2] && overlapind[j + (TAXON_SIZE + 1)] != overlapind[k + (TAXON_SIZE + 1) * 2])
								weight += F(overlap[i], overlap[j + (TAXON_SIZE + 1)], overlap[k + (TAXON_SIZE + 1) * 2]);
						}
					}
				}
				
				top = top + geneTreesAsInts4Ints.y + 1;
				
			}
		}
		else {
			break;
		}
		
		if(counter < geneTreesAsIntsLength) {
			counter++;
			if (newTree) {
				newTree = 0;

				allsides[0] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster1);
				allsides[1] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster2);
				allsides[2] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster3);

				treeCounter++;

			}
			if (geneTreesAsInts4Ints.z >= 0) {
				stack[top * globallen3 + idx] = ((trip.cluster1[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.z / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.z % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen + idx] = ((trip.cluster2[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.z / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.z % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen2 + idx] = ((trip.cluster3[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.z / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.z % LONG_BIT_LENGTH)) & 1;
				top++;
			}
			else if (geneTreesAsInts4Ints.z == INT_MIN) {
				top = 0;
				newTree = 1;
			}
			else if (geneTreesAsInts4Ints.z == -2) {
				top--;
				int topminus1 = top - 1;
				int newSides0 = stack[top * globallen3 + idx] + stack[topminus1 * globallen3 + idx];
				int newSides1 = stack[top * globallen3 + globallen + idx] + stack[topminus1 * globallen3 + globallen + idx];
				int newSides2 = stack[top * globallen3 + globallen2 + idx] + stack[topminus1 * globallen3 + globallen2 + idx];
				
				int side3s0 = allsides[0] - newSides0;
				int side3s1 = allsides[1] - newSides1;
				int side3s2 = allsides[2] - newSides2;
				/*
				weight += 
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s2) +
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s1) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + idx], side3s2) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], side3s1) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s0);
				*/	
				weight += FF(stack[top * globallen3 + idx], stack[top * globallen3 + globallen + idx], stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0, side3s1, side3s2);	
				stack[topminus1 * globallen3 + idx] = newSides0;
				stack[topminus1 * globallen3 + globallen + idx] = newSides1;
				stack[topminus1 * globallen3 + globallen2 + idx] = newSides2;
				
			}
			else { //for polytomies
			
				int nzc[3];
				nzc[0] = nzc[1] = nzc[2] = 0;
				int newSides[3];
				newSides[0] = newSides[1] = newSides[2] = 0;
				
				for(int side = 0; side < 3; side++) {
					for(int i = top - 1; i >= top + geneTreesAsInts4Ints.z; i--) {
						if(stack[i * globallen3 + side * globallen + idx] > 0) {
							newSides[side] += stack[i * globallen3 + side * globallen + idx];
							overlap[nzc[side]+ side * (TAXON_SIZE + 1)] = stack[i * globallen3 + side * globallen + idx];
							overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = i;
							nzc[side]++;
						}
					}
					
					stack[top * globallen3 + side * globallen + idx] = allsides[side] - newSides[side];
					
					if(stack[top * globallen3 + side * globallen + idx] > 0) {
						overlap[nzc[side] + side * (TAXON_SIZE + 1)] = stack[top * globallen3 + side * globallen + idx];
						overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = top;
						nzc[side]++;					
					}
					stack[(top + geneTreesAsInts4Ints.z) * globallen3 + side * globallen + idx] = newSides[side];
				}
				
				for(int i = nzc[0] - 1; i >= 0; i--) {
					for(int j = nzc[1] - 1; j >= 0; j--) {
						for(int k = nzc[2] - 1; k >= 0; k--) {
							if(overlapind[i] != overlapind[j + (TAXON_SIZE + 1)] && overlapind[i] != overlapind[k + (TAXON_SIZE + 1) * 2] && overlapind[j + (TAXON_SIZE + 1)] != overlapind[k + (TAXON_SIZE + 1) * 2])
								weight += F(overlap[i], overlap[j + (TAXON_SIZE + 1)], overlap[k + (TAXON_SIZE + 1) * 2]);
						}
					}
				}
				
				top = top + geneTreesAsInts4Ints.z + 1;
				
			}
		}
		else {
			break;
		}
		
		if(counter < geneTreesAsIntsLength) {
			counter++;
			if (newTree) {
				newTree = 0;

				allsides[0] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster1);
				allsides[1] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster2);
				allsides[2] = bitIntersectionSize(&allArray[treeCounter * SPECIES_WORD_LENGTH], trip.cluster3);

				treeCounter++;

			}
			if (geneTreesAsInts4Ints.w >= 0) {
				stack[top * globallen3 + idx] = ((trip.cluster1[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.w / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.w % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen + idx] = ((trip.cluster2[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.w / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.w % LONG_BIT_LENGTH)) & 1;
				stack[top * globallen3 + globallen2 + idx] = ((trip.cluster3[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts4Ints.w / LONG_BIT_LENGTH])>>(geneTreesAsInts4Ints.w % LONG_BIT_LENGTH)) & 1;
				top++;
			}
			else if (geneTreesAsInts4Ints.w == INT_MIN) {
				top = 0;
				newTree = 1;
			}
			else if (geneTreesAsInts4Ints.w == -2) {
				top--;
				int topminus1 = top - 1;
				int newSides0 = stack[top * globallen3 + idx] + stack[topminus1 * globallen3 + idx];
				int newSides1 = stack[top * globallen3 + globallen + idx] + stack[topminus1 * globallen3 + globallen + idx];
				int newSides2 = stack[top * globallen3 + globallen2 + idx] + stack[topminus1 * globallen3 + globallen2 + idx];
				
				int side3s0 = allsides[0] - newSides0;
				int side3s1 = allsides[1] - newSides1;
				int side3s2 = allsides[2] - newSides2;
				/*
				weight += 
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s2) +
					F(stack[top * globallen3 + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s1) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + idx], side3s2) +
					F(stack[top * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], side3s1) +
					F(stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + globallen + idx], side3s0);
				*/	
				weight += FF(stack[top * globallen3 + idx], stack[top * globallen3 + globallen + idx], stack[top * globallen3 + globallen2 + idx], stack[topminus1 * globallen3 + idx], stack[topminus1 * globallen3 + globallen + idx], stack[topminus1 * globallen3 + globallen2 + idx], side3s0, side3s1, side3s2);	
				stack[topminus1 * globallen3 + idx] = newSides0;
				stack[topminus1 * globallen3 + globallen + idx] = newSides1;
				stack[topminus1 * globallen3 + globallen2 + idx] = newSides2;
				
			}
			else { //for polytomies
			
				int nzc[3];
				nzc[0] = nzc[1] = nzc[2] = 0;
				int newSides[3];
				newSides[0] = newSides[1] = newSides[2] = 0;
				
				for(int side = 0; side < 3; side++) {
					for(int i = top - 1; i >= top + geneTreesAsInts4Ints.w; i--) {
						if(stack[i * globallen3 + side * globallen + idx] > 0) {
							newSides[side] += stack[i * globallen3 + side * globallen + idx];
							overlap[nzc[side]+ side * (TAXON_SIZE + 1)] = stack[i * globallen3 + side * globallen + idx];
							overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = i;
							nzc[side]++;
						}
					}
					
					stack[top * globallen3 + side * globallen + idx] = allsides[side] - newSides[side];
					
					if(stack[top * globallen3 + side * globallen + idx] > 0) {
						overlap[nzc[side] + side * (TAXON_SIZE + 1)] = stack[top * globallen3 + side * globallen + idx];
						overlapind[nzc[side] + side * (TAXON_SIZE + 1)] = top;
						nzc[side]++;					
					}
					stack[(top + geneTreesAsInts4Ints.w) * globallen3 + side * globallen + idx] = newSides[side];
				}
				
				for(int i = nzc[0] - 1; i >= 0; i--) {
					for(int j = nzc[1] - 1; j >= 0; j--) {
						for(int k = nzc[2] - 1; k >= 0; k--) {
							if(overlapind[i] != overlapind[j + (TAXON_SIZE + 1)] && overlapind[i] != overlapind[k + (TAXON_SIZE + 1) * 2] && overlapind[j + (TAXON_SIZE + 1)] != overlapind[k + (TAXON_SIZE + 1) * 2])
								weight += F(overlap[i], overlap[j + (TAXON_SIZE + 1)], overlap[k + (TAXON_SIZE + 1) * 2]);
						}
					}
				}
				
				top = top + geneTreesAsInts4Ints.w + 1;
				
			}
		}
		else {
			break;
		}
	}

	weightArray[idx] = weight;
}
