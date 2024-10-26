#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <stdbool.h>

int compare(const void *a, const void *b) {
    int int_a = *(int *)a;
    int int_b = *(int *)b;
    return (int_a > int_b) - (int_a < int_b);
}

void merge(int *A, int *B, int left, int mid, int right) {
    int i = left, j = mid, k = left;
    while (i < mid && j < right) {
        if (A[i] <= A[j]) B[k++] = A[i++];
        else B[k++] = A[j++];
    }
    while (i < mid) B[k++] = A[i++];
    while (j < right) B[k++] = A[j++];
    for (i = left; i < right; i++) A[i] = B[i];
}

void parallel_merge_sort(int *A, int *B, int left, int right, int depth) {
    if (right - left <= 1) return;
    
    if (depth > 0) {
        int mid = (left + right) / 2;
        
        #pragma omp task shared(A, B) if (depth > 0)
        parallel_merge_sort(A, B, left, mid, depth - 1);
        
        #pragma omp task shared(A, B) if (depth > 0)
        parallel_merge_sort(A, B, mid, right, depth - 1);
        
        #pragma omp taskwait
        merge(A, B, left, mid, right);
    } else {
        qsort(A + left, right - left, sizeof(int), compare);
    }
}

bool arrays_are_equal(int *arr1, int *arr2, int N) {
    for (int i = 0; i < N; i++) {
        if (arr1[i] != arr2[i]) {
            return false;
        }
    }
    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <array_size> <num_threads>\n", argv[0]);
        return -1;
    }

    int N = atoi(argv[1]);
    int p = atoi(argv[2]);

    int *A = malloc(N * sizeof(int));
    int *B = malloc(N * sizeof(int));
    int *C = malloc(N * sizeof(int));
    if (!A || !B || !C) {
        perror("Failed to allocate memory");
        return -1;
    }

    // Инициализация массива случайными значениями с сидом по времени
    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        A[i] = rand();
    }

    // Копирование массива для qsort
    for (int i = 0; i < N; i++) {
        C[i] = A[i];
    }

    // Открытие файла для записи результатов
    FILE *file = fopen("result.txt", "w");
    if (file == NULL) {
        perror("Failed to open result.txt");
        free(A);
        free(B);
        free(C);
        return -1;
    }

    // Время выполнения qsort
    double start_time = omp_get_wtime();
    qsort(C, N, sizeof(int), compare);
    double end_time = omp_get_wtime();
    double qsort_time = end_time - start_time;
    
    printf("qsort time:\t%f seconds\n\n", qsort_time);
    fprintf(file, "qsort time:\t%f seconds\n\n", qsort_time);

    // Копируем исходный массив A, чтобы сортировать заново при каждом p
    for (int j = 0; j < N; j++) {
        B[j] = A[j];
    }

    // Время параллельной сортировки с введенным значением p
    start_time = omp_get_wtime();
    
    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        parallel_merge_sort(B, C, 0, N, p);
    }

    end_time = omp_get_wtime();
    double parallel_time = end_time - start_time;

    // Проверка идентичности массивов
    bool is_equal = arrays_are_equal(B, C, N);
    double percentage = (parallel_time / qsort_time) * 100;

    // Вывод результатов на экран и в файл
    printf("Threads:\t%d\n", p);
    printf("Parallel:\t%f seconds\n", parallel_time);
    printf("Par/qsort:\t%.2f%%\n", percentage);
    printf("Arr %s\n\n", is_equal ? "equal" : "diff");

    fprintf(file, "Threads:\t%d\n", p);
    fprintf(file, "Parallel:\t%f seconds\n", parallel_time);
    fprintf(file, "Par/qsort:\t%.2f%%\n", percentage);
    fprintf(file, "Arr %s\n\n", is_equal ? "equal" : "diff");

    // Дополнительные тесты для p = 1, 2, 4, 8, 16
    int thread_counts[] = {1, 2, 4, 8, 16};
    
    // Массивы для хранения результатов
    double times[5];
    double speeds[5];
    double efficiencies[5];

    for (int i = 0; i < sizeof(thread_counts) / sizeof(thread_counts[0]); i++) {
        p = thread_counts[i];

        // Копируем исходный массив A, чтобы сортировать заново при каждом p
        for (int j = 0; j < N; j++) {
            B[j] = A[j];
        }

        // Время параллельной сортировки
        start_time = omp_get_wtime();
        
        #pragma omp parallel num_threads(p)
        {
            #pragma omp single
            parallel_merge_sort(B, C, 0, N, p);
        }

        end_time = omp_get_wtime();
        double parallel_time = end_time - start_time;

        // Проверка идентичности массивов
        is_equal = arrays_are_equal(B, C, N);
        percentage = (parallel_time / qsort_time) * 100;

        // Вычисление ускорения и эффективности
        double S = qsort_time / parallel_time;  // Ускорение
        double E = S / p;  // Эффективность

        // Сохранение результатов для таблицы
        times[i] = parallel_time;
        speeds[i] = S;
        efficiencies[i] = E;

        // Вывод результатов на экран и в файл
        printf("Threads:\t%d\n", p);
        printf("Parallel:\t%f seconds\n", parallel_time);
        printf("Par/qsort:\t%.2f%%\n", percentage);
        printf("Arr %s\n\n", is_equal ? "equal" : "diff");

        fprintf(file, "Threads:\t%d\n", p);
        fprintf(file, "Parallel:\t%f seconds\n", parallel_time);
        fprintf(file, "Par/qsort:\t%.2f%%\n", percentage);
        fprintf(file, "Arr %s\n\n", is_equal ? "equal" : "diff");
    }

    // Вывод таблицы результатов
    printf("\nResults Table:\n");
    printf("Threads\tTime (s)\tSpeedup\tEfficiency\n");
    fprintf(file, "\nResults Table:\n");
    fprintf(file, "Threads\tTime (s)\tSpeedup\tEfficiency\n");

    for (int i = 0; i < sizeof(thread_counts) / sizeof(thread_counts[0]); i++) {
        printf("%d\t%f\t%f\t%f\n", thread_counts[i], times[i], speeds[i], efficiencies[i]);
        fprintf(file, "%d\t%f\t%f\t%f\n", thread_counts[i], times[i], speeds[i], efficiencies[i]);
    }

    // Закрытие файла и освобождение памяти
    fclose(file);
    free(A);
    free(B);
    free(C);
    
    return 0;
}
