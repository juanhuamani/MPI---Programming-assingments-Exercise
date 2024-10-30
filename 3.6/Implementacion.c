#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void read_matrix(int n, int *A, int rank, int comm_sz, int submatrix_size);
void distribute_vector(int *vector, int *local_vector, int n, int submatrix_size, int rank, int sqrt_comm_sz);
void gather_result(int *local_result, int n, int rank, int comm_sz);

int main(int argc, char *argv[]) {
    int n = 4; 
    int comm_sz, my_rank;
    int sqrt_comm_sz, submatrix_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    sqrt_comm_sz = (int) sqrt(comm_sz);
    submatrix_size = n / sqrt_comm_sz;

    int *local_A = (int *)malloc(submatrix_size * submatrix_size * sizeof(int));
    int *local_vector = (int *)malloc(submatrix_size * sizeof(int));
    int *local_result = (int *)malloc(submatrix_size * sizeof(int));

    if (my_rank == 0) {
        int *A = (int *)malloc(n * n * sizeof(int));
        int *vector = (int *)malloc(n * sizeof(int));
        
        read_matrix(n, A, my_rank, comm_sz, submatrix_size);
        for (int i = 0; i < n; i++) vector[i] = 1;  // Vector ejemplo

        for (int i = 0; i < sqrt_comm_sz; i++) {
            for (int j = 0; j < sqrt_comm_sz; j++) {
                int dest = i * sqrt_comm_sz + j;
                if (dest == 0) {
                    // Copia local para el proceso 0
                    for (int k = 0; k < submatrix_size; k++)
                        for (int l = 0; l < submatrix_size; l++)
                            local_A[k * submatrix_size + l] = A[k * n + l];
                } else {
                    MPI_Send(&A[(i * n + j) * submatrix_size], submatrix_size * submatrix_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
                }
            }
        }
        
        distribute_vector(vector, local_vector, n, submatrix_size, my_rank, sqrt_comm_sz);

        free(A);
        free(vector);
    } else {
        MPI_Recv(local_A, submatrix_size * submatrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        distribute_vector(NULL, local_vector, n, submatrix_size, my_rank, sqrt_comm_sz);
    }

    for (int i = 0; i < submatrix_size; i++) {
        local_result[i] = 0;
        for (int j = 0; j < submatrix_size; j++) {
            local_result[i] += local_A[i * submatrix_size + j] * local_vector[j];
        }
    }

    gather_result(local_result, n, my_rank, comm_sz);

    free(local_A);
    free(local_vector);
    free(local_result);

    MPI_Finalize();
    return 0;
}

void read_matrix(int n, int *A, int rank, int comm_sz, int submatrix_size) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i * n + j] = 1; 
}

void distribute_vector(int *vector, int *local_vector, int n, int submatrix_size, int rank, int sqrt_comm_sz) {
    int diag_process = rank / sqrt_comm_sz == rank % sqrt_comm_sz;
    if (rank == 0) {
        for (int i = 0; i < sqrt_comm_sz; i++) {
            int dest = i * sqrt_comm_sz + i;
            if (dest == 0) {
                for (int j = 0; j < submatrix_size; j++)
                    local_vector[j] = vector[j];
            } else {
                MPI_Send(&vector[i * submatrix_size], submatrix_size, MPI_INT, dest, 0, MPI_COMM_WORLD);
            }
        }
    } else if (diag_process) {
        MPI_Recv(local_vector, submatrix_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void gather_result(int *local_result, int n, int rank, int comm_sz) {
    if (rank == 0) {
        int *result = (int *)malloc(n * sizeof(int));
        for (int i = 0; i < n / comm_sz; i++) result[i] = local_result[i];
        
        for (int i = 1; i < comm_sz; i++) {
            MPI_Recv(&result[i * (n / comm_sz)], n / comm_sz, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        printf("Resultado de la multiplicación matriz-vector:\n");
        for (int i = 0; i < n; i++) printf("%d ", result[i]);
        printf("\n");
        free(result);
    } else {
        MPI_Send(local_result, n / comm_sz, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}