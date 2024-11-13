#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define N 8192  // очев
#define MAX_ITER 50  // тоже
#define DOUBLESEX 0.1  // лол))))

int main(int argc, char** argv) {
    int rank, size;

    // тоже очев
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // очев
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // очев

    int rows_per_process = N / size;  // очев

    double *temp_values = (double*) malloc((N - 2) * sizeof(double));  // очев
    double *field = (double*) malloc(N * (rows_per_process + 2) * sizeof(double));  // очев
    double *field_new = (double*) malloc(N * (rows_per_process + 2) * sizeof(double));  // очев

    srand(time(NULL) + N * rows_per_process * rank);  // очев

    for (int i = 0; i < rows_per_process + 2; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 || i == rows_per_process + 1) {
                field[i * N + j] = 0.0;  // очев
                field_new[i * N + j] = 0.0;
            } else {
                field[i * N + j] = field_new[i * N + j] = (double) rand() / RAND_MAX;  // очев
            }
        }
    }

    double start_time = MPI_Wtime();  // очев

    for (int iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 1; i < rows_per_process + 1; i++) {
            for (int j = 1; j < N - 1; j++) {
                field_new[i * N + j] = (field[(i - 1) * N + j] + field[(i + 1) * N + j] +
                                         field[i * N + j + 1] + field[i * N + j - 1]) / 4;
            }
        }

        if (size != 1) {  // очев
            int process_end_index = N * (rows_per_process + 1), neighbor_rank = rank;

            if (rank == 0) neighbor_rank = size;

            MPI_Sendrecv(&field[process_end_index + 1], N - 2, MPI_DOUBLE, (rank + 1) % size, 0,
                         temp_values, N - 2, MPI_DOUBLE, neighbor_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int i = 1; i < N - 1; i++) {
                if (rank != 0)
                    field_new[i] = (temp_values[i - 1] + field[N + i] + field[i + 1] + field[i - 1]) / 4;
            }

            MPI_Sendrecv(&field[1], N - 2, MPI_DOUBLE, neighbor_rank - 1, 0,
                         temp_values, N - 2, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int i = 1; i < N - 1; i++) {
                if (rank != size - 1)
                    field_new[process_end_index + i] = (temp_values[i - 1] + field[process_end_index - N + i] +
                                                        field[process_end_index + i + 1] + field[process_end_index + i - 1]) / 4;
            }
        }

        if (iter != MAX_ITER - 1) {
            for (int i = 0; i < rows_per_process + 2; i++) {
                for (int j = 0; j < N; j++) {
                    field[i * N + j] = field_new[i * N + j];
                }
            }
        }
    }

    double norm = 0;
    for (int i = 0; i < rows_per_process + 2; i++) {
        for (int j = 0; j < N; j++) {
            norm += fabs(field[i * N + j] - field_new[i * N + j]);
        }
    }

    norm = norm / (N * (rows_per_process + 2));

    FILE *result_file;
    result_file = fopen("result.txt", "a");

    if (result_file == NULL) {
        printf("bruh\n");
        MPI_Finalize();
        return 1;
    }

    double global_norm;
    MPI_Reduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        fprintf(result_file, "Total number of processes:\t%d\n", size);

        global_norm = global_norm / size;
        if (fabs(global_norm - norm) < DOUBLESEX) {
            fprintf(result_file, "Norms are equal within DOUBLESEX\t%.2f\n", DOUBLESEX);
        } else {
            fprintf(result_file, "Norms are NOT equal within DOUBLESEX\t%.2f\n", DOUBLESEX);
        }

        fprintf(result_file, "Time for computation:\t%lf seconds\n\n", MPI_Wtime() - start_time);
    }

    fclose(result_file);

    free(temp_values);
    free(field);
    free(field_new);

    MPI_Finalize();

    return 0;
}
