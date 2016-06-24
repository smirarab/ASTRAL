#include "calculateWeightByTraversal.h"
#include <limits.h>
#include <stdio.h>

//testing purposes only. getside will probably use a read only memory

struct Intersects * getSide(int in, struct Intersects * side, struct Tripartition * trip) {
	/*int index = 0;
	while (index < CLUSTER_SIZE && in != trip->cluster1[index])
		index++;
	if (index < CLUSTER_SIZE) {
		index = 1;
	}
	else {
		index = 0;
		while (index < CLUSTER_SIZE && in != trip->cluster2[index])
			index++;
		if (index < CLUSTER_SIZE) {
			index = 2;
		}
		else
			index = 3;
	}*/
	
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

//for intersections
addIntersects(struct Intersects * augend, struct Intersects * addend, struct Intersects * result) {
	result->s0 = augend->s0 + addend->s0;
	result->s1 = augend->s1 + addend->s1;
	result->s2 = augend->s2 + addend->s2;
}
subtractIntersects(struct Intersects * minuend, struct Intersects * subtrahend, struct Intersects * result) {
	result->s0 = minuend->s0 - subtrahend->s0;
	result->s1 = minuend->s1 - subtrahend->s1;
	result->s2 = minuend->s2 - subtrahend->s2;
}

//for intersection stack
inline push(struct IntersectsStack * stack, struct Intersects * item) {
	stack->array[++(stack->currentIndex)].s0 = item->s0;
	stack->array[(stack->currentIndex)].s1 = item->s1;
	stack->array[(stack->currentIndex)].s2 = item->s2;
}

inline pop(struct IntersectsStack * stack, struct Intersects * item) {
	item->s0 = stack->array[stack->currentIndex].s0;
	item->s1 = stack->array[stack->currentIndex].s1;
	item->s2 = stack->array[stack->currentIndex--].s2;
}

poll(struct IntersectsStack * stack, struct Intersects * item) {
	item->s0 = stack->array[stack->currentIndex].s0;
	item->s1 = stack->array[stack->currentIndex].s1;
	item->s2 = stack->array[stack->currentIndex].s2;
}

clear(struct IntersectsStack * stack) {
	stack->currentIndex = -1;
}

//F function
inline long F(int a, int b, int c) {
	return ((long)(a + b + c - 3))*a*b*c;
}

//there is a popcount method in opencl. refer to opencl reference guide. will do it with naive method here. Assumed sorted.
int bitIntersectionSize(long input1[CLUSTER_SIZE], long input2[CLUSTER_SIZE]) {
	int out = 0;
	/*int input2Counter = 0;
	int i = 0;
	while (i < CLUSTER_SIZE && input1[i] != -1 && input2[input2Counter] != -1) {
		if (input1[i] > input2[input2Counter]) {
			input2Counter++;
		}
		else if (input1[i] < input2[input2Counter]) {
			i++;
		}
		else {
			out++;
			input2Counter++;
		}
	}*/
	for (int i = 0; i < CLUSTER_SIZE; i++) {
		if (input1[i] == input2[i] && input1[i] != 0) {
			out++;
		}
	}
	return out;
}

long calculateWeightByTraversal(struct Tripartition * trip, int all[][CLUSTER_SIZE], int allLength, int geneTreesAsInts[], int geneTreeAsIntsLength) {

	//actual program
	long weight = 0;

	struct Intersects allsides;

	//should be set to 0 and be actually set in the if newtree statement. set to some values for testing purposes.
	allsides.s0 = 0;
	allsides.s1 = 0;
	allsides.s2 = 0;

	struct IntersectsStack stack;
	stack.currentIndex = -1;

	int newTree = 1;
	int counter = 0;
	int treeCounter = 0;
	for (; counter < geneTreeAsIntsLength; counter++) {

		if (newTree) {
			newTree = 0;

			allsides.s0 = bitIntersectionSize(all[treeCounter], trip->cluster1);
			allsides.s1 = bitIntersectionSize(all[treeCounter], trip->cluster2);
			allsides.s2 = bitIntersectionSize(all[treeCounter], trip->cluster3);

			treeCounter++;

		}
		if (geneTreesAsInts[counter] >= 0) {
			struct Intersects side;
			push(&stack, getSide(geneTreesAsInts[counter], &side, trip));
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
			F(1, 2, 3);
		}
	}
	return weight;
}