#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

double calculateEuclideanDistance(double x1, double y1, double x2, double y2);
double **generateDistanceMatrix(double **coords, int numOfCoords);
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);

// Function to calculate the total cost of a tour (as in previous implementations)
double calculateTourCost(double **distanceMatrix, int *tour, int numOfCoords) {
    // Implementation remains the same as in previous algorithms
}

void nearestAddition(double **distanceMatrix, int numOfCoords, int *tour) {
    bool *visited = (bool *)calloc(numOfCoords, sizeof(bool));
    int tourSize = 1;

    // Start with vertex 0 (or any other vertex of your choice)
    tour[0] = 0;
    visited[0] = true;

    while (tourSize < numOfCoords) {
        double minDistance = -1;
        int nearestVertex = -1;
        int insertPos = -1;

        for (int i = 0; i < numOfCoords; ++i) {
            if (!visited[i]) {
                for (int j = 0; j < tourSize; ++j) {
                    double distance = distanceMatrix[tour[j]][i];
                    if (minDistance < 0 || distance < minDistance) {
                        minDistance = distance;
                        nearestVertex = i;
                        insertPos = j;
                    }
                }
            }
        }

        // Insert nearestVertex after insertPos
        for (int k = tourSize; k > insertPos; --k) {
            tour[k] = tour[k - 1];
        }
        tour[insertPos + 1] = nearestVertex;
        visited[nearestVertex] = true;
        tourSize++;
    }

    free(visited);
}

int main(int argc, char *argv[]) {
    // Similar setup as in previous implementations
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
        nearestAddition(distanceMatrix, numOfCoords, tour);
        
        double currentTourCost = calculateTourCost(distanceMatrix, tour, numOfCoords);
        if (bestTourCost < 0 || currentTourCost < bestTourCost) {
            bestTourCost = currentTourCost;
            memcpy(bestTour, tour, numOfCoords * sizeof(int));
        }
    }

    char *outNearestFileName = argv[4];
    writeTourToFile(bestTour, numOfCoords, outNearestFileName);

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
