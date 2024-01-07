#include <stdio.h>
#include <stdlib.h>

double calculateEuclideanDistance(double x1, double y1, double x2, double y2);
double **generateDistanceMatrix(double **coords, int numOfCoords);
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void writeTourToFile(int *tour, int tourLength, char *filename);

void cheapestInsertion(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
void farthestInsertion(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
void nearestAddition(double **distanceMatrix, int numOfCoords, int *tour);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <coord_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *coordFileName = argv[1];
    int numOfCoords = readNumOfCoords(coordFileName);
    double **coords = readCoords(coordFileName, numOfCoords);
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    int *tour = (int *)malloc(numOfCoords * sizeof(int));
    if (!tour) {
        perror("Memory allocation for tour failed");
        exit(EXIT_FAILURE);
    }

    // Iterate over all vertices as starting points for each algorithm
    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        printf("Starting Vertex: %d\n", startVertex);

        // Cheapest Insertion
        cheapestInsertion(distanceMatrix, numOfCoords, tour, startVertex);
        printf("Cheapest Insertion Tour: ");
        for (int i = 0; i < numOfCoords; i++) {
            printf("%d ", tour[i]);
        }
        printf("\n");

        // Farthest Insertion
        farthestInsertion(distanceMatrix, numOfCoords, tour, startVertex);
        printf("Farthest Insertion Tour: ");
        for (int i = 0; i < numOfCoords; i++) {
            printf("%d ", tour[i]);
        }
        printf("\n");

        // Nearest Addition
        nearestAddition(distanceMatrix, numOfCoords, tour);
        printf("Nearest Addition Tour: ");
        for (int i = 0; i < numOfCoords; i++) {
            printf("%d ", tour[i]);
        }
        printf("\n\n");
    }

    // Cleanup
    free(tour);
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
        free(coords[i]);
    }
    free(distanceMatrix);
    free(coords);

    return 0;
}
