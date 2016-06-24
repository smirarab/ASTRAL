#include "main.h"
#include <stdio.h>
int main()
{

	FILE * geneStream = fopen("geneTreesAsInts.txt", "r");
	FILE * tripStream = fopen("Tripartitions.txt", "r");

	if (geneStream == NULL) {
		printf("geneStream null");
		return;
	}
	if (tripStream == NULL) {
		printf("tripStream null");
		return;
	}

	//getting the geneTreesAsInts
	int geneTreesAsInts[99999];
	int geneTreesAsIntsCounter = 0;
	char integerString[BUFSIZ];
	int integerStringCounter = 0;
	int characterObtained;
	while ((characterObtained = fgetc(geneStream)) != -1) {
		if (characterObtained == ' ') {
			integerString[integerStringCounter] = 0;
			geneTreesAsInts[geneTreesAsIntsCounter++] = atoi(integerString);
			integerStringCounter = 0;
		}
		integerString[integerStringCounter++] = characterObtained;
	}
	if (integerStringCounter != 0) {
		integerString[integerStringCounter] = 0;
		geneTreesAsInts[geneTreesAsIntsCounter++] = atoi(integerString);
		integerStringCounter = 0;
	}

	//getting the all array
	int all[BUFSIZ][CLUSTER_SIZE] = { 0 };
	int allCounter = 0;
	int integerObtained = 0;
	while ((characterObtained = fgetc(tripStream)) != '&') {
		if (characterObtained == '}') {
			if (integerStringCounter != 0) {
				integerString[integerStringCounter] = 0;
				integerObtained = atoi(integerString);

				all[allCounter][integerObtained] = 1;
				integerStringCounter = 0;
			}
			allCounter++;
		}
		if (characterObtained <= '9' && characterObtained >= '0') {
			integerString[integerStringCounter++] = characterObtained;
		}
		if (characterObtained == ',') {
			integerString[integerStringCounter] = 0;
			integerObtained = atoi(integerString);

			all[allCounter][integerObtained] = 1;
			integerStringCounter = 0;
		}
	}

	integerStringCounter = 0;
	//cycling through the tripartitions
	struct Tripartition trip;
	for (int i = 0; i < CLUSTER_SIZE; i++) {
		trip.cluster1[i] = 0;
		trip.cluster2[i] = 0;
		trip.cluster3[i] = 0;
	}

	int tripCounter = 1;
	while ((characterObtained = fgetc(tripStream)) != -1) {
		if (characterObtained == '}') {
			if (integerStringCounter != 0) {
				integerString[integerStringCounter] = 0;
				integerObtained = atoi(integerString);
				if (tripCounter == 1) {
					trip.cluster1[integerObtained] = 1;
				}
				if (tripCounter == 2) {
					trip.cluster2[integerObtained] = 1;
				}
				if (tripCounter == 3) {
					trip.cluster3[integerObtained] = 1;
					characterObtained = fgetc(tripStream);
					while (!(characterObtained <= '9' && characterObtained >= '0')) {
						characterObtained = fgetc(tripStream);
					}
					integerStringCounter = 0;
					while (characterObtained <= '9' && characterObtained >= '0') {
						integerString[integerStringCounter++] = characterObtained;
						characterObtained = fgetc(tripStream);
					}
					integerString[integerStringCounter] = 0;

					long weight = atoi(integerString);
					long calculated = calculateWeightByTraversal(&trip, &all, allCounter, geneTreesAsInts, geneTreesAsIntsCounter);
					if (weight != calculated) {
						return -1;
					}
					tripCounter = 0;
					for (int i = 0; i < CLUSTER_SIZE; i++) {
						trip.cluster1[i] = 0;
						trip.cluster2[i] = 0;
						trip.cluster3[i] = 0;
					}
				}
				integerStringCounter = 0;
			}
			tripCounter++;
			
		}
		if (characterObtained <= '9' && characterObtained >= '0') {
			integerString[integerStringCounter++] = characterObtained;
		}
		if (characterObtained == ',') {
			integerString[integerStringCounter] = 0;
			integerObtained = atoi(integerString);

			if (tripCounter == 1) {
				trip.cluster1[integerObtained] = 1;
			}
			if (tripCounter == 2) {
				trip.cluster2[integerObtained] = 1;
			}
			if (tripCounter == 3) {
				trip.cluster3[integerObtained] = 1;
			}
			integerStringCounter = 0;
		}
	}
	fclose(geneStream);
	fclose(tripStream);

	return 0;
}
