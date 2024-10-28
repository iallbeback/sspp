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

void heapify(int *arr, int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest])
        largest = left;

    if (right < n && arr[right] > arr[largest])
        largest = right;

    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        heapify(arr, n, largest);
    }
}

void heap_sort(int *arr, int n) {
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    for (int i = n - 1; i > 0; i--) {
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;
        heapify(arr, i, 0);
    }
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
        heap_sort(A + left, right - left);
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

    srand(time(NULL));
    for (int i = 0; i < N; i++) {
        A[i] = rand();
    }

    for (int i = 0; i < N; i++) {
        C[i] = A[i];
    }

    FILE *file = fopen("result.txt", "w");
    if (file == NULL) {
        perror("Failed to open result.txt");
        free(A);
        free(B);
        free(C);
        return -1;
    }

    double start_time = omp_get_wtime();
    qsort(C, N, sizeof(int), compare);
    double end_time = omp_get_wtime();
    double qsort_time = end_time - start_time;
    
    printf("qsort time:\t%f seconds\n\n", qsort_time);
    fprintf(file, "qsort time:\t%f seconds\n\n", qsort_time);

    for (int j = 0; j < N; j++) {
        B[j] = A[j];
    }

    start_time = omp_get_wtime();
    
    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        parallel_merge_sort(B, C, 0, N, p);
    }

    end_time = omp_get_wtime();
    double parallel_time = end_time - start_time;

    bool is_equal = arrays_are_equal(B, C, N);
    double percentage = (parallel_time / qsort_time) * 100;

    printf("Threads:\t%d\n", p);
    printf("Parallel:\t%f seconds\n", parallel_time);
    printf("Par/qsort:\t%.2f%%\n", percentage);
    printf("Arr %s\n\n", is_equal ? "equal" : "diff");

    fprintf(file, "Threads:\t%d\n", p);
    fprintf(file, "Parallel:\t%f seconds\n", parallel_time);
    fprintf(file, "Par/qsort:\t%.2f%%\n", percentage);
    fprintf(file, "Arr %s\n\n", is_equal ? "equal" : "diff");

    int thread_counts[] = {1, 2, 4, 8, 16};
    
    double times[5];
    double speeds[5];
    double efficiencies[5];

    for (int i = 0; i < sizeof(thread_counts) / sizeof(thread_counts[0]); i++) {
        p = thread_counts[i];

        for (int j = 0; j < N; j++) {
            B[j] = A[j];
        }

        start_time = omp_get_wtime();
        
        #pragma omp parallel num_threads(p)
        {
            #pragma omp single
            parallel_merge_sort(B, C, 0, N, p);
        }

        end_time = omp_get_wtime();
        double parallel_time = end_time - start_time;

        is_equal = arrays_are_equal(B, C, N);
        percentage = (parallel_time / qsort_time) * 100;

        double S = qsort_time / parallel_time;
        double E = S / p;

        times[i] = parallel_time;
        speeds[i] = S;
        efficiencies[i] = E;

        printf("Threads:\t%d\n", p);
        printf("Parallel:\t%f seconds\n", parallel_time);
        printf("Par/qsort:\t%.2f%%\n", percentage);
        printf("Arr %s\n\n", is_equal ? "equal" : "diff");

        fprintf(file, "Threads:\t%d\n", p);
        fprintf(file, "Parallel:\t%f seconds\n", parallel_time);
        fprintf(file, "Par/qsort:\t%.2f%%\n", percentage);
        fprintf(file, "Arr %s\n\n", is_equal ? "equal" : "diff");
    }

    printf("\nResults Table:\n");
    printf("Threads\tTime (s)\tSpeedup\tEfficiency\n");
    fprintf(file, "\nResults Table:\n");
    fprintf(file, "Threads\tTime (s)\tSpeedup\tEfficiency\n");

    for (int i = 0; i < sizeof(thread_counts) / sizeof(thread_counts[0]); i++) {
        printf("%d\t%f\t%f\t%f\n", thread_counts[i], times[i], speeds[i], efficiencies[i]);
        fprintf(file, "%d\t%f\t%f\t%f\n", thread_counts[i], times[i], speeds[i], efficiencies[i]);
    }

    fclose(file);
    free(A);
    free(B);
    free(C);
    
    return 0;
}
