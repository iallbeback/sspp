#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define GRID_SIZE 4096       // Размер глобальной сетки
#define MAX_ITERATIONS 300   // Максимальное количество итераций

int main(int argc, char** argv) {
	srand(time(NULL));
    int process_rank, process_count;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank); // Получение ранга процесса
    MPI_Comm_size(MPI_COMM_WORLD, &process_count); // Получение количества процессов

    double start_time = MPI_Wtime();

    // Высота сегмента сетки для каждого процесса
    int segment_height = GRID_SIZE / process_count;
    int extended_width = GRID_SIZE + 2;   // Ширина с учетом обертки (границы)
    int extended_height = segment_height + 2; // Высота с учетом обертки (границы)

    // Выделение памяти для текущего и нового состояния сетки
    int *current_state = (int*) malloc(extended_width * extended_height * sizeof(int));
    int *next_state = (int*) malloc(extended_width * extended_height * sizeof(int));

    // Инициализация сетки нулями
    for (int i = 0; i < extended_width * extended_height; i++) {
        current_state[i] = next_state[i] = 0; //+ тут можно ранд 0-1
    }

    // Буферы для передачи данных о глайдерах
    int recv_buffer[3];  // Получаемая информация: X, Y, тип глайдера
    int send_buffer[process_count * 3]; // Отправляемая информация

    // Инициализация глайдеров только в нулевом процессе
    if (process_rank == 0) {
        for (int i = 0; i < process_count * 3; i += 3) {
            send_buffer[i] = extended_width / 2;      // Центр X
            send_buffer[i + 1] = extended_height / 2; // Центр Y
            send_buffer[i + 2] = rand() % 4;              // Тип глайдера + тут можно ранд 0-3
        }
    }

    // Распределение данных о глайдерах между процессами
    MPI_Scatter(send_buffer, 3, MPI_INT, recv_buffer, 3, MPI_INT, 0, MPI_COMM_WORLD);

    // Установка вертикальной линии для глайдера
    for (int i = -1; i < 2; i++) {
        current_state[(recv_buffer[1] + i) * extended_width + recv_buffer[0]] = 1;
    }

    // Установка дополнительных клеток глайдера в зависимости от типа
    int additional_cells[2];
    switch (recv_buffer[2]) {
        case 0: // Глайдер направлен вправо-вниз
            additional_cells[0] = recv_buffer[0] + 1 + (recv_buffer[1] + 1) * extended_width;
            additional_cells[1] = recv_buffer[0] + 2 + recv_buffer[1] * extended_width;
            break;
        case 1: // Глайдер направлен влево-вверх
            additional_cells[0] = recv_buffer[0] - 1 + (recv_buffer[1] - 1) * extended_width;
            additional_cells[1] = recv_buffer[0] - 2 + recv_buffer[1] * extended_width;
            break;
        case 2: // Глайдер направлен влево-вниз
            additional_cells[0] = recv_buffer[0] - 1 + (recv_buffer[1] + 1) * extended_width;
            additional_cells[1] = recv_buffer[0] - 2 + recv_buffer[1] * extended_width;
            break;
        case 3: // Глайдер направлен вправо-вверх
            additional_cells[0] = recv_buffer[0] + 1 + (recv_buffer[1] - 1) * extended_width;
            additional_cells[1] = recv_buffer[0] + 2 + recv_buffer[1] * extended_width;
            break;
    }

    current_state[additional_cells[0]] = 1;
    current_state[additional_cells[1]] = 1;

    MPI_Barrier(MPI_COMM_WORLD); // Синхронизация процессов

    int iteration = 0, previous_live_cells = -1, current_live_cells, new_live_cells;
    int continue_simulation = 1;

    // Обмен граничными данными между процессами
    MPI_Request send_up, send_down, recv_up, recv_down;

    while (continue_simulation) {
        // Асинхронный обмен границами
        MPI_Isend(&current_state[extended_width + 1], GRID_SIZE, MPI_INT, (process_rank + process_count - 1) % process_count, 0, MPI_COMM_WORLD, &send_up);
        MPI_Isend(&current_state[extended_width * (extended_height - 2)], GRID_SIZE, MPI_INT, (process_rank + 1) % process_count, 0, MPI_COMM_WORLD, &send_down);
        MPI_Irecv(&current_state[extended_width * (extended_height - 1) + 1], GRID_SIZE, MPI_INT, (process_rank + 1) % process_count, 0, MPI_COMM_WORLD, &recv_down);
        MPI_Irecv(&current_state[1], GRID_SIZE, MPI_INT, (process_rank + process_count - 1) % process_count, 0, MPI_COMM_WORLD, &recv_up);

        // Обновление основной части сетки
        for (int i = 2; i < extended_height - 2; i++) {
            for (int j = 2; j < extended_width - 2; j++) {
                int neighbors = current_state[(i + 1) * extended_width + j - 1] +
                                current_state[(i + 1) * extended_width + j] +
                                current_state[(i + 1) * extended_width + j + 1] +
                                current_state[i * extended_width + j - 1] +
                                current_state[i * extended_width + j + 1] +
                                current_state[(i - 1) * extended_width + j - 1] +
                                current_state[(i - 1) * extended_width + j] +
                                current_state[(i - 1) * extended_width + j + 1];

                // Правила жизни
                if (neighbors < 2 || neighbors > 3) {
                    next_state[i * extended_width + j] = 0; // Клетка умирает
                } else if (neighbors == 3) {
                    next_state[i * extended_width + j] = 1; // Новая клетка рождается
                } else {
                    next_state[i * extended_width + j] = current_state[i * extended_width + j]; // Клетка остаётся
                }
            }
        }

        // Ожидание завершения обмена границами
        MPI_Wait(&send_up, MPI_STATUS_IGNORE);
        MPI_Wait(&send_down, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_up, MPI_STATUS_IGNORE);
        MPI_Wait(&recv_down, MPI_STATUS_IGNORE);
		
		for (int i = 1; i < extended_height - 3; i = extended_height - 2) {
            for (int j = 1; j < extended_width - 3; j = extended_width - 2) {
                int neighbors = current_state[(i + 1) * extended_width + j - 1] +
                                current_state[(i + 1) * extended_width + j] +
                                current_state[(i + 1) * extended_width + j + 1] +
                                current_state[i * extended_width + j - 1] +
                                current_state[i * extended_width + j + 1] +
                                current_state[(i - 1) * extended_width + j - 1] +
                                current_state[(i - 1) * extended_width + j] +
                                current_state[(i - 1) * extended_width + j + 1];

                // Правила жизни
                if (neighbors < 2 || neighbors > 3) {
                    next_state[i * extended_width + j] = 0; // Клетка умирает
                } else if (neighbors == 3) {
                    next_state[i * extended_width + j] = 1; // Новая клетка рождается
                } else {
                    next_state[i * extended_width + j] = current_state[i * extended_width + j]; // Клетка остаётся
                }
            }
        }

        iteration++;
        // Вычисление количества живых клеток
        current_live_cells = 0;
        new_live_cells = 0;
        for (int i = 1; i < extended_height - 1; i++) {
            for (int j = 1; j < extended_width - 1; j++) {
                current_live_cells += current_state[i * extended_width + j];
                new_live_cells += next_state[i * extended_width + j];
            }
        }

        // Печать информации о живых клетках
        if (previous_live_cells != new_live_cells) {
            printf("Process %d, Iteration %d: Living Cells = %d\n", process_rank, iteration, new_live_cells);
            previous_live_cells = new_live_cells;
        }

        // Проверка завершения симуляции
        if (iteration > MAX_ITERATIONS && current_live_cells == new_live_cells) {
            continue_simulation = 0;
            MPI_Bcast(&continue_simulation, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }

        // Копирование нового состояния в текущее
        for (int i = 1; i < extended_height - 1; i++) {
            for (int j = 1; j < extended_width - 1; j++) {
                current_state[i * extended_width + j] = next_state[i * extended_width + j];
            }
        }

        MPI_Barrier(MPI_COMM_WORLD); // Синхронизация процессов
    }
	
	FILE * result_file;
	result_file = fopen("result.txt", "a");

    double end_time = MPI_Wtime();
    if (process_rank == 0) {
        printf("Simulation Time: %lf seconds\n", end_time - start_time);
		fprintf(result_file, "Process: %d\n", process_count);
		fprintf(result_file, "Simulation Time: %lf seconds\n", end_time - start_time);
    }

    // Суммирование всех живых клеток
    int total_live_cells = 0;
    MPI_Reduce(&new_live_cells, &total_live_cells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (process_rank == 0) {
        printf("Total Living Cells: %d\n", total_live_cells);
		fprintf(result_file, "Total Living Cells: %d\n", total_live_cells);
    }

    // Освобождение памяти
    free(current_state);
    free(next_state);

    MPI_Finalize();
    return 0;
}
