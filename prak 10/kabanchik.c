#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define ОБЩИЙ_РАЗМЕР 192 // Кабанчик обкашлил: это общий размер куба

#define КОЛИЧЕСТВО_ИТЕРАЦИЙ 50 // Решили вопросик: количество итераций

int main(int argc, char** argv) {
    int ранг, размер;
    MPI_Comm кубическая_сетка;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &ранг);
    MPI_Comm_size(MPI_COMM_WORLD, &размер);

    int размеры[3] = {0, 0, 0}, цикличность[3] = {0, 0, 0}, перераспределение = 0;

    // Кабанчик сказал: накидываем размеры сетки
    MPI_Dims_create(размер, 3, размеры);

    if (ранг == 0) printf("Топология: %d - %d - %d\n", размеры[0], размеры[1], размеры[2]);

    int локальный_размерX = ОБЩИЙ_РАЗМЕР / размеры[0] + 2, 
        локальный_размерY = ОБЩИЙ_РАЗМЕР / размеры[1] + 2, 
        локальный_размерZ = ОБЩИЙ_РАЗМЕР / размеры[2] + 2;
    int локальные_размеры[3] = {локальный_размерX, локальный_размерY, локальный_размерZ};
    int объем_куба = локальный_размерX * локальный_размерY * локальный_размерZ;
    double* поле = (double*) malloc(объем_куба * sizeof(double));
    double* новое_поле = (double*) malloc(объем_куба * sizeof(double));

    srand(time(NULL) + объем_куба * ранг);
    for (int ячейка = 0; ячейка < объем_куба; ячейка++) {
        новое_поле[ячейка] = (double) rand() / RAND_MAX; // Метнулся по делам и задал начальные значения
    }

    MPI_Cart_create(MPI_COMM_WORLD, 3, размеры, цикличность, перераспределение, &кубическая_сетка);

    int новый_ранг, координаты[3];

    MPI_Comm_rank(кубическая_сетка, &новый_ранг);
    MPI_Cart_coords(кубическая_сетка, новый_ранг, 3, координаты);

    MPI_Datatype тип_слоя[3];
    MPI_Type_vector(локальный_размерZ * локальный_размерY, 1, локальный_размерX, MPI_DOUBLE, &тип_слоя[0]); // Кабанчик решил: это слой по X
    MPI_Type_vector(локальный_размерZ, локальный_размерX, локальный_размерX * локальный_размерY, MPI_DOUBLE, &тип_слоя[1]); // Слой по Y
    MPI_Type_vector(1, локальный_размерX * локальный_размерY, 0, MPI_DOUBLE, &тип_слоя[2]); // Слой по Z
    for (int i = 0; i < 3; i++) MPI_Type_commit(&тип_слоя[i]);

    MPI_Request запросы[12];
    int флаги_соседей[12] = {0};

    int ранги_соседей[6], координаты_соседей[3];
    for (int i = 0; i < 3; i++) координаты_соседей[i] = координаты[i];

    int смещения[12] = {
        1, локальный_размерX - 2, локальный_размерX, 
        локальный_размерX * (локальный_размерY - 2), локальный_размерX * локальный_размерY, 
        локальный_размерX * локальный_размерY * (локальный_размерZ - 2), 
        0, локальный_размерX - 1, 0, 
        локальный_размерX * (локальный_размерY - 1), 0, 
        локальный_размерX * локальный_размерY * (локальный_размерZ - 1)
    };

    double начальное_время = MPI_Wtime();

    for (int итерация = 0; итерация < КОЛИЧЕСТВО_ИТЕРАЦИЙ; итерация++) {
        for (int ячейка = 0; ячейка < объем_куба; ячейка++) {
            поле[ячейка] = новое_поле[ячейка]; // Кабанчик обкашлил: обновили значения
        }

        for (int ось = 0; ось < 3; ось++) {
            // Решили вопросик: общаемся с соседями
            if (координаты[ось] > 0) {
                координаты_соседей[ось]--;
                MPI_Cart_rank(кубическая_сетка, координаты_соседей, &ранги_соседей[ось * 2]);
                MPI_Isend(&поле[ось * 2], 1, тип_слоя[ось], ранги_соседей[ось * 2], 0, кубическая_сетка, &запросы[ось * 2]);
                MPI_Irecv(&поле[смещения[ось * 2 + 6]], 1, тип_слоя[ось], ранги_соседей[ось * 2], 0, кубическая_сетка, &запросы[6 + ось * 2]);
                координаты_соседей[ось]++;
                флаги_соседей[ось * 2] = флаги_соседей[ось * 2 + 6] = 1;
            }
            if (координаты[ось] < размеры[ось] - 1) {
                координаты_соседей[ось]++;
                MPI_Cart_rank(кубическая_сетка, координаты_соседей, &ранги_соседей[ось * 2 + 1]);
                MPI_Isend(&поле[смещения[ось * 2 + 1]], 1, тип_слоя[ось], ранги_соседей[ось * 2 + 1], 0, кубическая_сетка, &запросы[ось * 2 + 1]);
                MPI_Irecv(&поле[смещения[ось * 2 + 7]], 1, тип_слоя[ось], ранги_соседей[ось * 2 + 1], 0, кубическая_сетка, &запросы[7 + ось * 2]);
                координаты_соседей[ось]--;
                флаги_соседей[ось * 2 + 1] = флаги_соседей[ось * 2 + 7] = 1;
            }
        }

        // Тут идёт подсчёт значений, чтобы кабанчик остался доволен
        // ...
    }

    double финальное_время = MPI_Wtime();

    MPI_Finalize();

    if (ранг == 0) printf("Время Якоби: %lf\n", финальное_время - начальное_время);

    return 0;
}
