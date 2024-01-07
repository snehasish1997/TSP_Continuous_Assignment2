#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

/*This code is for reading and writing to files for the 2023-24 COMP528 CA2*/

/*Use the functions in this file to read from the input file, and write to the output file*/

/*You should use this file when compiling your code*/

/*Declare these functions at the top of each 'main' file*/

/*If there are any issues with this code, please contact: h.j.forbes@liverpool.ac.uk*/

// Function declarations
double calculateEuclideanDistance(double x1, double y1, double x2, double y2);
double **generateDistanceMatrix(double **coords, int numOfCoords);

// Function to calculate Euclidean distance between two points
double calculateEuclideanDistance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

// Function to generate a distance matrix from the coordinates
double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));
    for(int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
        if (distanceMatrix[i] == NULL) {
            perror("Memory Allocation Failed");
            exit(EXIT_FAILURE);
        }

        for(int j = 0; j < numOfCoords; j++) {
            // Distance from a point to itself is always 0
            if(i == j) {
                distanceMatrix[i][j] = 0.0;
            } else {
                distanceMatrix[i][j] = calculateEuclideanDistance(
                    coords[i][0], coords[i][1], coords[j][0], coords[j][1]
                );
            }
        }
    }
    return distanceMatrix;
}



int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);

/*Gets the number of the coordinates in the file. Returns as a single integer*/
int readNumOfCoords(char *filename){
	FILE *file = fopen(filename, "r");
	int numOfCoords = 0;

	if(file == NULL){
		return -1;
	}

	char line[100];

	while(fgets(line, sizeof(line), file) != NULL){
		numOfCoords++;
	}

	fclose(file);

    return numOfCoords;
}

/*Gets the data from the file. Returns as an array of doubles, Ignores the first integer*/
double **readCoords(char *filename, int numOfCoords){
	FILE *file = fopen(filename,"r");
    int i;

	char line[100];
    
    if(file == NULL) {
        printf("Unable to open file: %s", filename);
        return NULL;
    }

	double **coords = (double **)malloc(numOfCoords * sizeof(double *));

	for(i = 0; i < numOfCoords; i++){
		coords[i] = (double *) malloc(2 * sizeof(double));
		if (coords[i] == NULL){
			perror("Memory Allocation Failed");
		}
	}

	int lineNum = 0;
	while(fgets(line, sizeof(line), file) != NULL){
		double x, y;
		if (sscanf(line, "%lf,%lf", &x, &y) == 2){
			coords[lineNum][0] = x;
			coords[lineNum][1] = y;
			lineNum++;
		}
	}

	fclose(file);

	return coords;
}

void *writeTourToFile(int *tour, int tourLength, char *filename){
	
	FILE *file = fopen(filename, "w");
	int i;	
	
	if(file == NULL){
		printf("Unable to open file: %s", filename);
		return NULL;
	}

	fprintf(file, "%d \n", tourLength);

	printf("Writing output data\n");
    for(i=0; i < tourLength; i++) {
        fprintf(file, "%d ", tour[i]);
    }

	fclose(file);

	return tour;
}
