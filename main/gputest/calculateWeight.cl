struct Intersects {
	int s0;
	int s1;
	int s2;
};
void addIntersects(struct Intersects * augend, struct Intersects * addend, struct Intersects * result) {
	result->s0 = augend->s0 + addend->s0;
	result->s1 = augend->s1 + addend->s1;
	result->s2 = augend->s2 + addend->s2;
}
void subtractIntersects(struct Intersects * minuend, struct Intersects * subtrahend, struct Intersects * result) {
	result->s0 = minuend->s0 - subtrahend->s0;
	result->s1 = minuend->s1 - subtrahend->s1;
	result->s2 = minuend->s2 - subtrahend->s2;
}
struct IntersectsStack {
	struct Intersects array [SPECIES_LENGTH];
	int currentIndex;
};
void push(struct IntersectsStack * stack, struct Intersects * item) {
	stack->array[++(stack->currentIndex)].s0 = item->s0;
	stack->array[(stack->currentIndex)].s1 = item->s1;
	stack->array[(stack->currentIndex)].s2 = item->s2;
}
void pop(struct IntersectsStack * stack, struct Intersects * item) {
	item->s0 = stack->array[stack->currentIndex].s0;
	item->s1 = stack->array[stack->currentIndex].s1;
	item->s2 = stack->array[stack->currentIndex--].s2;
}
clear(struct IntersectsStack * stack) {
	stack->currentIndex = -1;
}
long F(int a, int b, int c) {
	return ((long)(a + b + c - 3))*a*b*c;
}
struct Tripartition {
	__global int* cluster1;
	__global int* cluster2;
	__global int* cluster3;
};
struct Intersects * getSide(int in, struct Intersects * side, struct Tripartition * trip) {
	if (trip->cluster1[in]) {
		side->s0 = 1;
		side->s1 = 0;
		side->s2 = 0;
	}
	else if (trip->cluster2[in]) {
		side->s0 = 0;
		side->s1 = 1;
		side->s2 = 0;
	}
	else {
		side->s0 = 0;
		side->s1 = 0;
		side->s2 = 1;
	}

	return side;
}
int bitIntersectionSize(__global int input1[SPECIES_LENGTH], __global int input2[SPECIES_LENGTH]) {
	int out = 0;
	for (int i = 0; i < SPECIES_LENGTH; i++) {
		if (input1[i] == input2[i] && input1[i] != 0) {
			out++;
		}
	}
	return out;
}
__kernel void calcWeight(
	__global int* geneTreesAsInts,
	int geneTreesAsIntsLength,
	__global int* allArray,
	__global int* tripartitions,
	__global long* weightArray
){
	long weight = 0;
	struct Tripartition trip;
	int idx = get_global_id(0);
	trip.cluster1 = SPECIES_LENGTH * 3 * idx + tripartitions;
	trip.cluster2 = SPECIES_LENGTH * 3 * idx + SPECIES_LENGTH + tripartitions;
	trip.cluster3 = SPECIES_LENGTH * 3 * idx + SPECIES_LENGTH * 2 + tripartitions;
	
	struct Intersects allsides;
	allsides.s0 = 0;
	allsides.s1 = 0;
	allsides.s2 = 0;

	struct IntersectsStack stack;
	stack.currentIndex = -1;

	int newTree = 1;
	int counter = 0;
	int treeCounter = 0;
	for (; counter < geneTreesAsIntsLength; counter++) {

		if (newTree) {
			newTree = 0;

			allsides.s0 = bitIntersectionSize(&allArray[treeCounter * SPECIES_LENGTH], trip.cluster1);
			allsides.s1 = bitIntersectionSize(&allArray[treeCounter * SPECIES_LENGTH], trip.cluster2);
			allsides.s2 = bitIntersectionSize(&allArray[treeCounter * SPECIES_LENGTH], trip.cluster3);

			treeCounter++;

		}
		if (geneTreesAsInts[counter] >= 0) {
			struct Intersects side;
			push(&stack, getSide(geneTreesAsInts[counter], &side, &trip));
		}
		else if (geneTreesAsInts[counter] == INT_MIN) {
			clear(&stack);
			newTree = 1;
		}
		else if (geneTreesAsInts[counter] == -2) {
			struct Intersects side1;
			struct Intersects side2;
			struct Intersects newSide;
			struct Intersects side3;

			pop(&stack, &side1);
			pop(&stack, &side2);
			
			addIntersects(&side1, &side2, &newSide);

			push(&stack, &newSide);

			subtractIntersects(&allsides, &newSide, &side3);

			weight += 
				F(side1.s0, side2.s1, side3.s2) +
				F(side1.s0, side2.s2, side3.s1) +
				F(side1.s1, side2.s0, side3.s2) +
				F(side1.s1, side2.s2, side3.s0) +
				F(side1.s2, side2.s0, side3.s1) +
				F(side1.s2, side2.s1, side3.s0);
			/*
			weight += 
				((long)(side1.s0 + side2.s1 + side3.s2 - 3))*side1.s0*side2.s1*side3.s2 +
				((long)(side1.s0 + side2.s2 + side3.s1 - 3))*side1.s0*side2.s2*side3.s1 +
				((long)(side1.s1 + side2.s0 + side3.s2 - 3))*side1.s1*side2.s0*side3.s2 +
				((long)(side1.s1 + side2.s2 + side3.s0 - 3))*side1.s1*side2.s2*side3.s0 +
				((long)(side1.s2 + side2.s0 + side3.s1 - 3))*side1.s2*side2.s0*side3.s1 +
				((long)(side1.s2 + side2.s1 + side3.s0 - 3))*side1.s2*side2.s1*side3.s0;
				*/
			/*if(idx == 32800) {
				printf("|%d %d %d %d %d %d %d %d %d %d |", weight, side1.s0, side1.s1, side1.s2, side2.s0, side2.s1, side2.s2, side3.s0, side3.s1, side3.s2);
			}*/
		}
	}
	weightArray[idx] = weight;
}