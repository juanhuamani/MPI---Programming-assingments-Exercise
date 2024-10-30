#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    int my_rank, comm_sz;
    int local_value, global_sum;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    local_value = my_rank;

    int step = 1;
    while (step < comm_sz) {
        if (my_rank % (2 * step) == 0) {
            int received_value;
            MPI_Recv(&received_value, 1, MPI_INT, my_rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_value += received_value;
        } else if (my_rank % step == 0) {
            MPI_Send(&local_value, 1, MPI_INT, my_rank - step, 0, MPI_COMM_WORLD);
            break;
        }
        step *= 2;
    }

    if (my_rank == 0) {
        global_sum = local_value;
        printf("Global sum is %d\n", global_sum);
    }

    MPI_Finalize();
    return 0;
}