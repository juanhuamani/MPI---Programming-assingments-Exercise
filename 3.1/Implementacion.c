#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void Get_input(int *bin_count_p, float *min_meas_p, float *max_meas_p,
               int *data_count_p, int *local_data_count_p, int my_rank,
               int comm_sz, MPI_Comm comm);

void Gen_data(float local_data[], int local_data_count, int data_count,
              float min_meas, float max_meas, int my_rank, MPI_Comm comm);

void Set_bins(float bin_maxes[], int loc_bin_cts[], float min_meas,
              float max_meas, int bin_count, int my_rank, MPI_Comm comm);

void Find_bins(int bin_counts[], float local_data[], int loc_bin_cts[],
               int local_data_count, float bin_maxes[], int bin_count,
               float min_meas, MPI_Comm comm);

int Which_bin(float data, float bin_maxes[], int bin_count, float min_meas);

void Print_histo(float bin_maxes[], int bin_counts[], int bin_count,
                 float min_meas);

void e(int error);

int main(int argc, char *argv[])
{
    int bin_count, data_count, local_data_count;
    float min_meas, max_meas;
    float *bin_maxes, *data, *local_data;
    int *bin_counts, *loc_bin_cts;
    int my_rank, comm_sz;
    MPI_Comm comm;

    e(MPI_Init(&argc, &argv));
    comm = MPI_COMM_WORLD;
    e(MPI_Comm_size(comm, &comm_sz));
    e(MPI_Comm_rank(comm, &my_rank));

    Get_input(&bin_count, &min_meas, &max_meas, &data_count,
              &local_data_count, my_rank, comm_sz, comm);

    bin_maxes = malloc(bin_count * sizeof(float));
    bin_counts = malloc(bin_count * sizeof(int));
    loc_bin_cts = malloc(bin_count * sizeof(int));
    data = malloc(data_count * sizeof(float));
    local_data = malloc(local_data_count * sizeof(float));

    Set_bins(bin_maxes, loc_bin_cts, min_meas, max_meas, bin_count, my_rank, comm);
    Gen_data(local_data, local_data_count, data_count, min_meas, max_meas, my_rank, comm);
    Find_bins(bin_counts, local_data, loc_bin_cts, local_data_count,
              bin_maxes, bin_count, min_meas, comm);
    e(MPI_Reduce(loc_bin_cts, bin_counts, bin_count, MPI_INT, MPI_SUM, 0, comm));
    if (my_rank == 0)
        Print_histo(bin_maxes, bin_counts, bin_count, min_meas);
    free(bin_maxes);
    free(bin_counts);
    free(loc_bin_cts);
    free(data);
    free(local_data);

    MPI_Finalize();
    return 0;
}
void e(int error){
    if (error != MPI_SUCCESS)
    {
        fprintf(stderr, "Error al iniciar el programa MPI. Terminando.\n");
        MPI_Abort(MPI_COMM_WORLD, error);
        MPI_Finalize();
        exit(1);
    }
}
void Print_histo(float bin_maxes[], int bin_counts[], int bin_count, float min_meas){
    int width = 40, max = 0, row_width, i, j;

    // Obtener el conteo máximo
    for (i = 0; i < bin_count; i++){
        if (bin_counts[i] > max)
            max = bin_counts[i];
    }
    for (i = 0; i < bin_count; i++){
        printf("%10.3f |", bin_maxes[i]);
        row_width = (float)bin_counts[i] / max * width;
        for (j = 0; j < row_width; j++)
        {
            printf("#");
        }
        printf("  %d\n", bin_counts[i]);
    }
}

void Find_bins(int bin_counts[], float local_data[], int loc_bin_cts[],
               int local_data_count, float bin_maxes[], int bin_count,
               float min_meas, MPI_Comm comm){
    int i, bin;
    for (i = 0; i < local_data_count; i++)
    {
        bin = Which_bin(local_data[i], bin_maxes, bin_count, min_meas);
        loc_bin_cts[bin]++;
    }
}
int Which_bin(float data, float bin_maxes[], int bin_count, float min_meas){
    int i;
    for (i = 0; i < bin_count - 1; i++)
    {
        if (data <= bin_maxes[i])
            break;
    }
    return i;
}
void Set_bins(float bin_maxes[], int loc_bin_cts[], float min_meas,
              float max_meas, int bin_count, int my_rank, MPI_Comm comm){
    float range = max_meas - min_meas;
    float interval = range / bin_count;
    for (int i = 0; i < bin_count; i++)
    {
        bin_maxes[i] = interval * (i + 1) + min_meas;
        loc_bin_cts[i] = 0;
    }
}
void Gen_data(float local_data[], int local_data_count, int data_count,
              float min_meas, float max_meas, int my_rank, MPI_Comm comm){
    float *data;
    if (my_rank == 0)
    {
        float range = max_meas - min_meas;
        data = malloc(data_count * sizeof(float));
        for (int i = 0; i < data_count; i++)
        {
            data[i] = (float)rand() / RAND_MAX * range + min_meas;
        }
    }
    e(MPI_Scatter(data, local_data_count, MPI_FLOAT, local_data, local_data_count,
                  MPI_FLOAT, 0, comm));
    if (my_rank == 0)
        free(data);
}
void Get_input(int *bin_count_p, float *min_meas_p, float *max_meas_p,
               int *data_count_p, int *local_data_count_p, int my_rank,
               int comm_sz, MPI_Comm comm){
    if (my_rank == 0)
    {
        printf("Número de contenedores (int): ");
        scanf("%d", bin_count_p);
        printf("Valor mínimo (float): ");
        scanf("%f", min_meas_p);
        printf("Valor máximo (float): ");
        scanf("%f", max_meas_p);
        if (*max_meas_p < *min_meas_p)
        {
            float *temp = max_meas_p;
            max_meas_p = min_meas_p;
            min_meas_p = temp;
        }

        printf("Número de valores (int): ");
        scanf("%d", data_count_p);
        *local_data_count_p = *data_count_p / comm_sz;
        *data_count_p = *local_data_count_p * comm_sz;
        printf("\n");
    }

    e(MPI_Bcast(bin_count_p, 1, MPI_INT, 0, comm));
    e(MPI_Bcast(min_meas_p, 1, MPI_FLOAT, 0, comm));
    e(MPI_Bcast(max_meas_p, 1, MPI_FLOAT, 0, comm));
    e(MPI_Bcast(data_count_p, 1, MPI_INT, 0, comm));
    e(MPI_Bcast(local_data_count_p, 1, MPI_INT, 0, comm));
}