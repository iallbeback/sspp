#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define N 192
#define MAX_ITER 50  // Максимальное количество итераций

int main(int argc, char** argv) {
    int rank, size;
    MPI_Comm comm3d;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[3] = {0, 0, 0}, period[3] = {0, 0, 0}, reorder = 0;

    MPI_Dims_create(size, 3, dims);

    // Размеры локальной области в каждой из трех осей с учетом добавочных слоев
    int localNx = N / dims[0] + 2, localNy = N / dims[1] + 2, localNz = N / dims[2] + 2;
    int localDims[3] = {localNx, localNy, localNz};
    int localSize = localNx * localNy * localNz;

    // Массивы для данных
    double* data = (double*) malloc(localSize * sizeof(double));
    double* newData = (double*) malloc(localSize * sizeof(double));

    srand(time(NULL) + localSize * rank);

    // Инициализация массива случайными значениями
    for (int i = 0; i < localSize; i++) {
        newData[i] = (double) rand() / RAND_MAX;
    }

    // Создание трехмерной топологии процесса
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, period, reorder, &comm3d);

    int localRank, coords[3];
    MPI_Comm_rank(comm3d, &localRank);
    MPI_Cart_coords(comm3d, localRank, 3, coords);

    // Определение типов данных для слоев
    MPI_Datatype layerType[3];
    MPI_Type_vector(localNz * localNy, 1, localNx, MPI_DOUBLE, &layerType[0]);
    MPI_Type_vector(localNz, localNx, localNx * localNy, MPI_DOUBLE, &layerType[1]);
    MPI_Type_vector(1, localNx * localNy, 0, MPI_DOUBLE, &layerType[2]);
    MPI_Type_commit(&layerType[0]);
    MPI_Type_commit(&layerType[1]);
    MPI_Type_commit(&layerType[2]);

    MPI_Request requests[12];
	for(int i = 0; i<12; i++) requests[i] = MPI_REQUEST_NULL;

    int neighborRanks[6], neighborCoords[3];
    for (int i = 0; i < 3; i++) neighborCoords[i] = coords[i];

    // Расстояния для отправки и получения данных между соседними процессами
    int displacements[12] = {1, localNx - 2, localNx, localNx * (localNy - 2), localNx * localNy, localNx * localNy * (localNz - 2),
                              0, localNx - 1, 0, localNx * (localNy - 1), 0, localNx * localNy * (localNz - 1)};
	//намёк понял

    double startTime = MPI_Wtime();

    // Основной цикл вычислений
    for (int n = 0; n < MAX_ITER; n++) {
        // Копирование новых данных в массив старых данных
        for (int i = 0; i < localSize; i++) {
            data[i] = newData[i];
        }

        // Обмен данными с соседями
        for (int i = 0; i < 3; i++) {
            if (coords[i] > 0) {
                neighborCoords[i] -= 1;
                MPI_Cart_rank(comm3d, neighborCoords, &neighborRanks[i * 2]);
                MPI_Isend(&data[displacements[i * 2]], 1, layerType[i], neighborRanks[i * 2], 0, comm3d, &requests[i * 2]);
                MPI_Irecv(&data[displacements[i * 2 + 6]], 1, layerType[i], neighborRanks[i * 2], 0, comm3d, &requests[6 + i * 2]);
                neighborCoords[i] += 1;
            }
            if (coords[i] < dims[i] - 1) {
                neighborCoords[i] += 1;
                MPI_Cart_rank(comm3d, neighborCoords, &neighborRanks[i * 2 + 1]);
                MPI_Isend(&data[displacements[i * 2 + 1]], 1, layerType[i], neighborRanks[i * 2 + 1], 0, comm3d, &requests[i * 2 + 1]);
                MPI_Irecv(&data[displacements[i * 2 + 7]], 1, layerType[i], neighborRanks[i * 2 + 1], 0, comm3d, &requests[7 + i * 2]);
                neighborCoords[i] -= 1;
            }
        }

        // Вычисление нового значения для каждой ячейки
        for (int i = 2; i < localNz - 2; i++) {
            for (int j = 2; j < localNy - 2; j++) {
                for (int k = 2; k < localNx - 2; k++) {
                    double xSum, ySum, zSum;
                    int id = i * localNx * localNy + j * localNx + k;
                    xSum = data[id + 1] + data[id - 1];
                    ySum = data[id + localNx] + data[id - localNx];
                    zSum = data[id + localNx * localNy] + data[id - localNx * localNy];
                    newData[id] = (xSum + ySum + zSum) / 6;
                }
            }
        }

        MPI_Waitall(12, requests, MPI_STATUS_IGNORE); // кабанчик подскакал

        // Повторение вычислений для крайних граней
        // Левые и правые грани
        for (int i = 1; i < localNz - 1; i++) {
            for (int j = 1; j < localNy - 1; j++) {
                for (int k = 1; k < localNx - 1; k = localNx - 1) {
                    double xSum, ySum, zSum;
                    int id = i * localNx * localNy + j * localNx + k;
                    xSum = data[id + 1] + data[id - 1];
                    ySum = data[id + localNx] + data[id - localNx];
                    zSum = data[id + localNx * localNy] + data[id - localNx * localNy];
                    newData[id] = (xSum + ySum + zSum) / 6;
                }
            }
        }

        // Верхняя и нижняя грани
        for (int i = 1; i < localNz - 1; i = localNz - 1) {
            for (int j = 1; j < localNy - 1; j++) {
                for (int k = 2; k < localNx - 2; k++) {
                    double xSum, ySum, zSum;
                    int id = i * localNx * localNy + j * localNx + k;
                    xSum = data[id + 1] + data[id - 1];
                    ySum = data[id + localNx] + data[id - localNx];
                    zSum = data[id + localNx * localNy] + data[id - localNx * localNy];
                    newData[id] = (xSum + ySum + zSum) / 6;
                }
            }
        }

        // Передняя и задняя грани
        for (int i = 2; i < localNz - 2; i++) {
            for (int j = 1; j < localNy - 1; j = localNy - 1) {
                for (int k = 2; k < localNx - 2; k++) {
                    double xSum, ySum, zSum;
                    int id = i * localNx * localNy + j * localNx + k;
                    xSum = data[id + 1] + data[id - 1];
                    ySum = data[id + localNx] + data[id - localNx];
                    zSum = data[id + localNx * localNy] + data[id - localNx * localNy];
                    newData[id] = (xSum + ySum + zSum) / 6;
                }
            }
        }
    }

    // Вычисление нормы
    double norm = 0;
    for (int i = 0; i < localSize; i++) {
        norm += fabs(data[i] - newData[i]);
    }

    norm = norm / localSize; //нормально так

    printf("Norm rank %d, coord %d-%d-%d: %lf\n", localRank, coords[0], coords[1], coords[2], norm);

    free(data);
    free(newData);

    double endTime = MPI_Wtime();
	
	MPI_Type_free(&layerType[0]);
    MPI_Type_free(&layerType[1]);
    MPI_Type_free(&layerType[2]);

    MPI_Finalize();

	FILE* yes_ass; // эту переменную не я назвал
	yes_ass = fopen("result.txt", "a");

    if (rank == 0) {
		printf("Time Ycobi: %lf\n", endTime - startTime);
		fprintf(yes_ass, "Time Ycobi: %lf\n", endTime - startTime);
	}

    return 0;
}
