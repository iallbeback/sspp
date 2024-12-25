#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#define MATRIX_SIZE 1536 // Размер матрицы

int main(int argc, char** argv) {
    int process_rank, total_processes;
    int block_size;
    MPI_Comm cartesian_comm, row_comm, col_comm;
    int grid_dimensions[2], periodicity[2], allow_reorder;

    // Инициализация MPI
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
	
	int grid_size = sqrt(total_processes);

    // Определение размера подматрицы для каждого процесса
    int local_matrix_size = MATRIX_SIZE / grid_size;

    // Главный процесс считывает размер блока из аргументов
    if (process_rank == 0) sscanf(argv[1], "%d", &block_size);

    // Распространение значения block_size среди всех процессов
    MPI_Bcast(&block_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int num_blocks = local_matrix_size / block_size;

    // Выделение памяти для локальных матриц
    float *matrix_A = (float*) malloc(local_matrix_size * local_matrix_size * sizeof(float));
    float *matrix_B = (float*) malloc(local_matrix_size * local_matrix_size * sizeof(float));
    float *matrix_C = (float*) malloc(local_matrix_size * local_matrix_size * sizeof(float));

    // Инициализация локальных матриц случайными числами
    srand(time(NULL) + local_matrix_size * local_matrix_size * 2 * process_rank);
    for (int i = 0; i < local_matrix_size * local_matrix_size; i++) {
        matrix_A[i] = (float) rand() / RAND_MAX;
        matrix_B[i] = (float) rand() / RAND_MAX;
        matrix_C[i] = 0.0;
    }

    // Создание декартовой топологии
    grid_dimensions[0] = grid_dimensions[1] = grid_size;
    periodicity[0] = periodicity[1] = 0; // Без замыкания
    allow_reorder = 0;

    MPI_Cart_create(MPI_COMM_WORLD, 2, grid_dimensions, periodicity, allow_reorder, &cartesian_comm);

    int cartesian_rank, coordinates[2];
    MPI_Comm_rank(cartesian_comm, &cartesian_rank);
    MPI_Cart_coords(cartesian_comm, cartesian_rank, 2, coordinates);

    // Создание подкоммуникаторов для строк и столбцов
    int remain_dims[2];
    remain_dims[0] = 0; remain_dims[1] = 1;
    MPI_Cart_sub(cartesian_comm, remain_dims, &row_comm);

    remain_dims[0] = 1; remain_dims[1] = 0;
    MPI_Cart_sub(cartesian_comm, remain_dims, &col_comm);

    int row_rank, col_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);

    // Начало замера времени
    double start_time = MPI_Wtime();

    // Алгоритм SUMMA
    for (int step = 0; step < grid_size; step++) {
        int active_row = step;
        int active_col = step;

        for (int block = 0; block < num_blocks; block++) {
            float *submatrix_A = (float*) malloc(local_matrix_size * block_size * sizeof(float));
            float *submatrix_B = (float*) malloc(local_matrix_size * block_size * sizeof(float));

            if (row_rank == active_row) {
                for (int i = 0; i < local_matrix_size; i++) {
                    for (int j = 0; j < block_size; j++) {
                        submatrix_A[i * block_size + j] = matrix_A[i * local_matrix_size + j + block * block_size];
                    }
                }
            }
            if (col_rank == active_col) {
                for (int i = 0; i < block_size; i++) {
                    for (int j = 0; j < local_matrix_size; j++) {
                        submatrix_B[i * local_matrix_size + j] = matrix_B[(i + block * block_size) * local_matrix_size + j];
                    }
                }
            }

            MPI_Bcast(submatrix_A, local_matrix_size * block_size, MPI_FLOAT, active_row, row_comm);
            MPI_Bcast(submatrix_B, local_matrix_size * block_size, MPI_FLOAT, active_col, col_comm);

            // Обновление матрицы C
            for (int i = 0; i < local_matrix_size; i++) {
                for (int j = 0; j < local_matrix_size; j++) {
                    for (int k = 0; k < block_size; k++) {
                        matrix_C[i * local_matrix_size + j] += submatrix_A[i * block_size + k] * submatrix_B[k * local_matrix_size + j];
                    }
                }
            }

            free(submatrix_A);
            free(submatrix_B);
        }
    }

    // Конец замера времени
    double end_time = MPI_Wtime();
	FILE* no_boobs;
	no_boobs = fopen("result.txt", "a");
    if (process_rank == 0){
		printf("Время выполнения SUMMA: %lf секунд\n", end_time - start_time);
		fprintf(no_boobs, "Время выполнения SUMMA: %lf секунд\n", end_time - start_time);
	}

    // Освобождение памяти
    free(matrix_A);
    free(matrix_B);
    free(matrix_C);

    MPI_Finalize();
    return 0;
}
