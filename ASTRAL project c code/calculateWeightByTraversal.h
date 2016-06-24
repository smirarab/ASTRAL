#pragma once
#define STACK_SIZE 200
#define CLUSTER_SIZE 200
#define BITSET_SIZE 1
long calculateWeightByTraversal(struct Tripartition * trip, int all[][CLUSTER_SIZE], int allLength, int geneTreesAsInts[], int geneTreeAsIntsLength);

//intersects
struct Intersects {
	int s0;
	int s1;
	int s2;
};

addIntersects(struct Intersects * augend, struct Intersects * addend, struct Intersects * result);
subtractIntersects(struct Intersects * minuend, struct Intersects * subtrahend, struct Intersects * result);

//intersects stack
struct IntersectsStack {
	struct Intersects array [STACK_SIZE];
	int currentIndex;
};

inline push(struct IntersectsStack * stack, struct Intersects * item);
inline pop(struct IntersectsStack * stack, struct Intersects * item);
poll(struct IntersectsStack * stack, struct Intersects * item);
clear(struct IntersectsStack * stack);

//tripartition/bitset
struct Tripartition {
	long cluster1[CLUSTER_SIZE];
	long cluster2[CLUSTER_SIZE];
	long cluster3[CLUSTER_SIZE];
};

int bitIntersectionSize(int input1[], int input2[]);

inline long F(int a, int b, int c);