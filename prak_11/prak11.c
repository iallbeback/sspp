#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define TOTAL_ELEMENTS 36864 // Размер общего массива (например, 4096*9)

int main(int argc, char** argv) {
    int process_rank, total_processes;

    // Инициализация MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank); // Определение ранга текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &total_processes); // Определение общего числа процессов

    int grid_dimensions[2] = {0, 0}; // Размеры процесса в 2D-топологии

    // Автоматическое определение размеров 2D-топологии
    MPI_Dims_create(total_processes, 2, grid_dimensions);

    if (process_rank == 0) 
        printf("Topology dimensions: %d x %d\n", grid_dimensions[0], grid_dimensions[1]);

    // Определение координат процесса в топологии
    int process_coordinates[2] = {process_rank % grid_dimensions[0], process_rank / grid_dimensions[0]};
    int rows_per_process = TOTAL_ELEMENTS / grid_dimensions[0]; // Число строк на процесс
    int cols_per_process = TOTAL_ELEMENTS / grid_dimensions[1]; // Число столбцов на процесс

    // Выделение памяти для локального подмассива
    double* local_matrix = (double*)malloc(rows_per_process * cols_per_process * sizeof(double));
    srand(time(NULL) + rows_per_process * cols_per_process * process_rank);
    for (int i = 0; i < rows_per_process * cols_per_process; i++) {
        local_matrix[i] = (double)rand() / RAND_MAX; // Инициализация случайными числами
    }

    // Выделение памяти для глобальных массивов
    double* global_vector = (double*)malloc(TOTAL_ELEMENTS * sizeof(double)); // Глобальный вектор
    double* result_vector = (double*)malloc(TOTAL_ELEMENTS * sizeof(double)); // Вектор результата

    // Инициализация вектора результата
    for (int i = 0; i < TOTAL_ELEMENTS; i++) {
        result_vector[i] = 0.0;
    }

    // Инициализация глобального вектора на процессе с рангом 0
    if (process_rank == 0) {
        srand(time(NULL) + total_processes * TOTAL_ELEMENTS * TOTAL_ELEMENTS);
        for (int i = 0; i < TOTAL_ELEMENTS; i++) {
            global_vector[i] = (double)rand() / RAND_MAX;
        }
    }

    // Создание окон MPI для работы с RMA
    MPI_Win global_vector_window, result_vector_window;
    MPI_Win_create(global_vector, TOTAL_ELEMENTS * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &global_vector_window);
    MPI_Win_create(result_vector, TOTAL_ELEMENTS * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &result_vector_window);

    // Начало RMA-операции для получения глобального вектора
    MPI_Win_fence(0, global_vector_window);

    // Каждый процесс получает часть глобального вектора
    if (process_rank != 0) {
        MPI_Get(&global_vector[process_coordinates[0] * rows_per_process], rows_per_process, MPI_DOUBLE, 0, process_coordinates[0] * rows_per_process, rows_per_process, MPI_DOUBLE, global_vector_window);
    }

    MPI_Win_fence(0, global_vector_window);

    // Начало расчётов
    double start_time = MPI_Wtime();

    // Умножение локальной матрицы на часть глобального вектора
    for (int i = 0; i < cols_per_process; i++) {
        for (int j = 0; j < rows_per_process; j++) {
            result_vector[i + process_coordinates[1] * cols_per_process] += local_matrix[i * rows_per_process + j] * global_vector[j + process_coordinates[0] * rows_per_process];
        }
    }

    // Синхронизация результатов с использованием RMA
    MPI_Win_fence(0, result_vector_window);

    if (process_rank != 0) {
        MPI_Accumulate(&result_vector[process_coordinates[1] * cols_per_process], cols_per_process, MPI_DOUBLE, 0, process_coordinates[1] * cols_per_process, cols_per_process, MPI_DOUBLE, MPI_SUM, result_vector_window);
    }

    MPI_Win_fence(0, result_vector_window);

    double end_time = MPI_Wtime();
	
	FILE* z_end = fopen("result.txt", "a");

    // Вывод времени выполнения на процессе с рангом 0
    if (process_rank == 0) {
        printf("Time for RMA operations: %lf seconds\n", end_time - start_time);
		fprintf(z_end, "Time for RMA operations: %lf seconds\n", end_time - start_time);
    }

	// Очистка памяти
    free(local_matrix);
    free(global_vector);
    free(result_vector);
    MPI_Win_free(&global_vector_window);
    MPI_Win_free(&result_vector_window);

    // Завершение работы MPI
    MPI_Finalize();

    return 0;
}
