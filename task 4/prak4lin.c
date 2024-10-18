#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>

void init_matrix(double* mat, int N) {
    for (int i = 0; i < N * N; i++) {
        mat[i] = (double)rand() / RAND_MAX; 
    }
}

void matmul_sequential(double* A, double* B, double* C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

void matmul_avx(double* A, double* B, double* C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            __m256d c_vec = _mm256_setzero_pd();

            for (int k = 0; k < N; k += 4) {
                __m256d a_vec = _mm256_loadu_pd(&A[i * N + k]);
                __m256d b_vec = _mm256_loadu_pd(&B[k * N + j]);
                c_vec = _mm256_add_pd(c_vec, _mm256_mul_pd(a_vec, b_vec));
            }

            _mm256_storeu_pd(&C[i * N + j], c_vec);
        }
    }
}

void benchmark(int N) {
    double *A, *B, *C_seq, *C_avx;

    posix_memalign((void**)&A, 32, N * N * sizeof(double));
    posix_memalign((void**)&B, 32, N * N * sizeof(double));
    posix_memalign((void**)&C_seq, 32, N * N * sizeof(double));
    posix_memalign((void**)&C_avx, 32, N * N * sizeof(double));

    init_matrix(A, N);
    init_matrix(B, N);

    clock_t start_seq = clock();
    matmul_sequential(A, B, C_seq, N);
    clock_t end_seq = clock();
    double time_taken_seq = (double)(end_seq - start_seq) / CLOCKS_PER_SEC;
    printf("N = %d, seq mull:\t %f secs\n", N, time_taken_seq);

    clock_t start_avx = clock();
    matmul_avx(A, B, C_avx, N);
    clock_t end_avx = clock();
    double time_taken_avx = (double)(end_avx - start_avx) / CLOCKS_PER_SEC;
    printf("N = %d, AVX mull:\t %f secs\n", N, time_taken_avx);

    free(A);
    free(B);
    free(C_seq);
    free(C_avx);
}

int main() {
    srand(time(NULL));

    int sizes[] = {512, 1024, 2048}; 
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    for (int i = 0; i < num_sizes; i++) {
        int N = sizes[i];
        benchmark(N);
    }

    return 0;
}
