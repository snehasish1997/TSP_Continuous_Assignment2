#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

double calculateEuclideanDistance(double x1, double y1, double x2, double y2);
double **generateDistanceMatrix(double **coords, int numOfCoords);
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);

// Function to calculate the total cost of a tour
double calculateTourCost(double **distanceMatrix, int *tour, int numOfCoords) {
    double cost = 0.0;
    for (int i = 0; i < numOfCoords - 1; ++i) {
        cost += distanceMatrix[tour[i]][tour[i + 1]];
    }
    cost += distanceMatrix[tour[numOfCoords - 1]][tour[0]]; // Closing the tour loop
    return cost;
}

void cheapestInsertion(double **distanceMatrix, int numOfCoords, int *tour, int startVertex) {
    bool *visited = (bool *)calloc(numOfCoords, sizeof(bool));
    int tourSize = 1;
    tour[0] = startVertex;
    visited[startVertex] = true;

    while (tourSize < numOfCoords) {
        double minCost = -1;
        int minVertex = -1;
        int insertPos = -1;

        for (int i = 0; i < numOfCoords; ++i) {
            if (!visited[i]) {
                for (int j = 0; j < tourSize; ++j) {
                    int current = tour[j];
                    int next = tour[(j + 1) % tourSize];
                    double cost = distanceMatrix[current][i] + distanceMatrix[i][next] - distanceMatrix[current][next];
                    
                    if (minCost < 0 || cost < minCost) {
                        minCost = cost;
                        minVertex = i;
                        insertPos = j;
                    }
                }
            }
        }

        for (int k = tourSize; k > insertPos; --k) {
            tour[k] = tour[k - 1];
        }
        tour[insertPos + 1] = minVertex;
        visited[minVertex] = true;
        tourSize++;
    }

    free(visited);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <coord_file> <out_cheapest> <out_farthest> <out_nearest>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *coordFileName = argv[1];
    int numOfCoords = readNumOfCoords(coordFileName);
    double **coords = readCoords(coordFileName, numOfCoords);
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    int *tour = (int *)malloc(numOfCoords * sizeof(int));
    int *bestTour = (int *)malloc(numOfCoords * sizeof(int));
    if (!tour || !bestTour) {
        perror("Memory allocation for tour failed");
        exit(EXIT_FAILURE);
    }

    double bestTourCost = -1;
    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        cheapestInsertion(distanceMatrix, numOfCoords, tour, startVertex);
        
        double currentTourCost = calculateTourCost(distanceMatrix, tour, numOfCoords);
        if (bestTourCost < 0 || currentTourCost < bestTourCost) {
            bestTourCost = currentTourCost;
            memcpy(bestTour, tour, numOfCoords * sizeof(int));
        }
    }

    char *outCheapestFileName = argv[2];
    writeTourToFile(bestTour, numOfCoords, outCheapestFileName);

    free(tour);
    free(bestTour);
    for (int i = 0; i < numOfCoords; ++i) {
        free(distanceMatrix[i]);
        free(coords[i]);
    }
    free(distanceMatrix);
    free(coords);

    return 0;
}
