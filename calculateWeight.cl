long F(int a, int b, int c) {
	return ((long)(a + b + c - 3))*a*b*c;
}
struct Tripartition {
	__global long* cluster1;
	__global long* cluster2;
	__global long* cluster3;
};
int bitIntersectionSize(__global long input1[SPECIES_WORD_LENGTH], __global long input2[SPECIES_WORD_LENGTH]) {
	int out = 0;
	for (int i = 0; i < SPECIES_WORD_LENGTH; i++) {
		out += popcount(input1[i]&input2[i]);
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

	int children[(STACK_SIZE+2)*3];
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
			if (((trip.cluster1[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts[counter] / LONG_BIT_LENGTH])>>(geneTreesAsInts[counter] % LONG_BIT_LENGTH)) & 1) {
				children[top] = 1;
				children[top + (STACK_SIZE+2)] = 0;
				children[top + (STACK_SIZE+2)*2] = 0;
			}
			else if (((trip.cluster2[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts[counter] / LONG_BIT_LENGTH])>>(geneTreesAsInts[counter] % LONG_BIT_LENGTH)) & 1) {
				children[top] = 0;
				children[top + (STACK_SIZE+2)] = 1;
				children[top + (STACK_SIZE+2)*2] = 0;
			}
			else if (((trip.cluster3[SPECIES_WORD_LENGTH - 1 - geneTreesAsInts[counter] / LONG_BIT_LENGTH])>>(geneTreesAsInts[counter] % LONG_BIT_LENGTH)) & 1){
				children[top] = 0;
				children[top + (STACK_SIZE+2)] = 0;
				children[top + (STACK_SIZE+2)*2] = 1;
			}
			else {
				children[top] = 0;
				children[top + (STACK_SIZE+2)] = 0;
				children[top + (STACK_SIZE+2)*2] = 0;
			}
			top++;
		}
		else if (geneTreesAsInts[counter] == INT_MIN) {
			top = 0;
			newTree = 1;
		}
		else if (geneTreesAsInts[counter] == -2) {
			top--;
	
			int newSides0 = children[top] + children[top - 1];
			int newSides1 = children[top + (STACK_SIZE+2)] + children[top - 1 + (STACK_SIZE+2)];
			int newSides2 = children[top + (STACK_SIZE+2)*2] + children[top - 1 + (STACK_SIZE+2)*2];
			
			int side3s0 = allsides[0] - newSides0;
			int side3s1 = allsides[1] - newSides1;
			int side3s2 = allsides[2] - newSides2;

			weight += 
				F(children[top], children[top - 1 + (STACK_SIZE + 2)], side3s2) +
				F(children[top], children[top - 1 + (STACK_SIZE + 2)*2], side3s1) +
				F(children[top + (STACK_SIZE + 2)], children[top - 1], side3s2) +
				F(children[top + (STACK_SIZE + 2)], children[top - 1 + (STACK_SIZE + 2)*2], side3s0) +
				F(children[top + (STACK_SIZE + 2)*2], children[top - 1], side3s1) +
				F(children[top + (STACK_SIZE + 2)*2], children[top - 1 + (STACK_SIZE + 2)], side3s0);
				
			children[top - 1] = newSides0;
			children[top - 1 + (STACK_SIZE + 2)] = newSides1;
			children[top - 1 + (STACK_SIZE + 2)*2] = newSides2;
			
			/*
			weight += 
				((long)(children[top] + children[top - 1 + (STACK_SIZE + 2)] + side3.s2 - 3))*children[top]*children[top - 1 + (STACK_SIZE + 2)]*side3.s2 +
				((long)(children[top] + children[top - 1 + (STACK_SIZE + 2)*2] + side3.s1 - 3))*children[top]*children[top - 1 + (STACK_SIZE + 2)*2]*side3.s1 +
				((long)(children[top + (STACK_SIZE + 2)] + children[top - 1] + side3.s2 - 3))*children[top + (STACK_SIZE + 2)]*children[top - 1]*side3.s2 +
				((long)(children[top + (STACK_SIZE + 2)] + children[top - 1 + (STACK_SIZE + 2)*2] + side3.s0 - 3))*children[top + (STACK_SIZE + 2)]*children[top - 1 + (STACK_SIZE + 2)*2]*side3.s0 +
				((long)(children[top + (STACK_SIZE + 2)*2] + children[top - 1] + side3.s1 - 3))*children[top + (STACK_SIZE + 2)*2]*children[top - 1]*side3.s1 +
				((long)(children[top + (STACK_SIZE + 2)*2] + children[top - 1 + (STACK_SIZE + 2)] + side3.s0 - 3))*children[top + (STACK_SIZE + 2)*2]*children[top - 1 + (STACK_SIZE + 2)]*side3.s0;
				*/
			/*if(idx == 32800) {
				printf("|%d %d %d %d %d %d %d %d %d %d |", weight, children[top], children[top + (STACK_SIZE + 2)], children[top + (STACK_SIZE + 2)*2], children[top - 1], children[top - 1 + (STACK_SIZE + 2)], children[top - 1 + (STACK_SIZE + 2)*2], side3.s0, side3.s1, side3.s2);
			}*/
		}
		else { //for polytomies
		
			int len = -geneTreesAsInts[counter] + 1;

			int newSides0 = 0;
			int newSides1 = 0;
			int newSides2 = 0;
			
			for (int i = top - 1; i >= top + geneTreesAsInts[counter]; i--) {
				newSides0 = children[i];
				newSides1 = children[i + (STACK_SIZE+2)];
				newSides2 = children[i + (STACK_SIZE+2)*2];

			}
			
			children[top] = allsides[0] - newSides0;
			children[top + (STACK_SIZE+2)] = allsides[1] - newSides1;
			children[top + (STACK_SIZE+2)*2] = allsides[2] - newSides2;
			for (int i = top; i >= top + geneTreesAsInts[counter]; i--) {
				if(children[i] == 0)
					continue;
				for (int j = top; j >= top + geneTreesAsInts[counter]; j--) {
					if(children[j+(STACK_SIZE+2)] == 0 || i == j)
						continue;
					for (int k = top; k >= top + geneTreesAsInts[counter]; k--) {
						if(children[k+(STACK_SIZE+2)*2] == 0 || i == k || j == k)
							continue;
						weight += F(children[i], children[j+(STACK_SIZE+2)], children[k+(STACK_SIZE+2)*2]);
					}
				}
			}
			
			top = top + geneTreesAsInts[counter] + 1;
			
			children[top - 1] = newSides0;
			children[top - 1 + (STACK_SIZE + 2)] = newSides1;
			children[top - 1 + (STACK_SIZE + 2)*2] = newSides2;
		}
	}
	weightArray[idx] = weight;
}
