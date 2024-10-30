#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
    int rank, size, tosses_per_process, total_tosses;
    long long number_in_circle = 0, global_number_in_circle = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        printf("Ingrese el número total de lanzamientos (tosses): ");
        scanf("%d", &total_tosses);
    }

    MPI_Bcast(&total_tosses, 1, MPI_INT, 0, MPI_COMM_WORLD);
    tosses_per_process = total_tosses / size;
    srand(rank + 1);

    for (int toss = 0; toss < tosses_per_process; toss++) {
        double x = (double)rand() / RAND_MAX * 2 - 1;
        double y = (double)rand() / RAND_MAX * 2 - 1;
        if (x * x + y * y <= 1) {
            number_in_circle++;
        }
    }

    MPI_Reduce(&number_in_circle, &global_number_in_circle, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi_estimate = 4.0 * global_number_in_circle / total_tosses;
        printf("Estimación de π: %f\n", pi_estimate);
    }

    MPI_Finalize();
    return 0;
}