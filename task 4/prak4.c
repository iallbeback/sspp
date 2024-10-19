#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <immintrin.h>

void init_matrix(double* mat, int N) {
    for (int i = 0; i < N * N; i++) {
        mat[i] = (double)rand() / RAND_MAX; 
    }
}

void print_matrix(double* mat, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f \t", mat[i * N + j]);
        }
        printf("\n");
    }
}

void print_matrices_and_diff_to_file(double* mat1, double* mat2, int N, const char* filename) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file %s for writing!\n", filename);
        return;
    }

    fprintf(file, "Matrix 1 (Sequential result):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", mat1[i * N + j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "\nMatrix 2 (AVX result):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", mat2[i * N + j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "\nMatrix Difference (Matrix 1 - Matrix 2):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double diff = mat1[i * N + j] - mat2[i * N + j];
            fprintf(file, "%f ", diff);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

void print_matrices_and_diff_to_file_d(double* mat1, double* mat2, int N, const char* filename) {
    FILE* file = fopen(filename, "a");
    if (file == NULL) {
        printf("Error opening file %s for writing!\n", filename);
        return;
    }

    fprintf(file, "Matrix 1 (Sequential result):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", mat1[i * N + j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "\nMatrix 2 (AVX result):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", mat2[i * N + j]);
        }
        fprintf(file, "\n");
    }

    fprintf(file, "\nMatrix Difference (Matrix 1 - Matrix 2):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double diff = mat1[i * N + j] - mat2[i * N + j];
            fprintf(file, "%f ", diff);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

int are_matrices_equal(double* mat1, double* mat2, int N) {
    for (int i = 0; i < N * N; i++) {
		double dif = mat1[i] - mat2[i];
		if (dif < 0) dif *= -1;
        if (dif > 0.00001) {
            return 0;
        }
    }
    return 1;
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

void matmul_diagonal(double* A, double* B, double* C, int N) {
    for (int d = 0; d < 2 * N - 1; d++) {
        int i_start = (d < N) ? 0 : d - N + 1;
        int i_end = (d < N) ? d : N - 1;

        for (int i = i_start; i <= i_end; i++) {
            int j = d - i;

            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

void transpose(double* B, double* B_T, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B_T[j * N + i] = B[i * N + j];
        }
    }
}

void matmul_avx(double* A, double* B, double* B_T, double* C, int N) {
	transpose(B, B_T, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            __m256d c_vec = _mm256_setzero_pd();
            for (int k = 0; k < N; k += 4) {
                __m256d a_vec = _mm256_loadu_pd(&A[i * N + k]);
                __m256d b_vec = _mm256_loadu_pd(&B_T[j * N + k]);
                c_vec = _mm256_add_pd(c_vec, _mm256_mul_pd(a_vec, b_vec));
            }

            //double c[4];
            //_mm256_storeu_pd(c, c_vec);
            //C[i * N + j] = c[0] + c[1] + c[2] + c[3];
			__m256d hsum = _mm256_hadd_pd(c_vec, c_vec);
            double result = ((double*)&hsum)[0] + ((double*)&hsum)[2];
            C[i * N + j] = result;
        }
    }
}

void matmul_avx_d(double* A, double* B, double* B_T, double* C, int N) {
	transpose(B, B_T, N);
    for (int d = 0; d < 2 * N - 1; d++) {
        int i_start = (d < N) ? 0 : d - N + 1;
        int i_end = (d < N) ? d : N - 1;

        for (int i = i_start; i <= i_end; i++) {
            int j = d - i;  

            __m256d c_vec = _mm256_setzero_pd();

            for (int k = 0; k < N; k += 4) {  
                __m256d a_vec = _mm256_loadu_pd(&A[i * N + k]);
                __m256d b_vec = _mm256_loadu_pd(&B_T[j * N + k]);

                c_vec = _mm256_add_pd(c_vec, _mm256_mul_pd(a_vec, b_vec));
            }

            //double temp[4];
            //_mm256_storeu_pd(temp, c_vec);
            //C[i * N + j] = temp[0] + temp[1] + temp[2] + temp[3];
			__m256d hsum = _mm256_hadd_pd(c_vec, c_vec);
            double result = ((double*)&hsum)[0] + ((double*)&hsum)[2];
            C[i * N + j] = result;
        }
    }
}

void benchmark(int N) {
    double *A, *B, *B_T, *B_T_d, *C_seq, *C_avx, *C_avx_d, *C_diag;

    int r = posix_memalign((void**)&A, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&B, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&B_T, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&B_T_d, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&C_seq, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&C_avx, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&C_avx_d, 32, N * N * sizeof(double));
    r = posix_memalign((void**)&C_diag, 32, N * N * sizeof(double));

    init_matrix(A, N);
    init_matrix(B, N);

    double start_seq = omp_get_wtime();
    matmul_sequential(A, B, C_seq, N);
    double end_seq = omp_get_wtime();
    double time_taken_seq = end_seq - start_seq;
    printf("N = %d, seq mull:\t %f secs\n", N, time_taken_seq);
	
	double start_diag = omp_get_wtime();
    matmul_diagonal(A, B, C_diag, N);
    double end_diag = omp_get_wtime();
    double time_taken_diag = end_diag - start_diag;
    printf("N = %d, diag mull:\t %f secs\n", N, time_taken_diag);

    double start_avx = omp_get_wtime();
    matmul_avx(A, B, B_T, C_avx, N);
    double end_avx = omp_get_wtime();
    double time_taken_avx = end_avx - start_avx;
    printf("N = %d, AVX mull:\t %f secs\n", N, time_taken_avx);

    double start_avx_d = omp_get_wtime();
    matmul_avx_d(A, B, B_T_d, C_avx_d, N);
    double end_avx_d = omp_get_wtime();
    double time_taken_avx_d = end_avx_d - start_avx_d;
    printf("N = %d, AVX diag mull:\t %f secs\n", N, time_taken_avx_d);

    if (N == 8) {
		print_matrices_and_diff_to_file(C_seq, C_diag, N, "matrices_and_diff.txt");
        print_matrices_and_diff_to_file_d(C_seq, C_avx, N, "matrices_and_diff.txt");
        print_matrices_and_diff_to_file_d(C_seq, C_avx_d, N, "matrices_and_diff.txt");
    }

    if (are_matrices_equal(C_seq, C_avx, N) &&
        are_matrices_equal(C_seq, C_avx_d, N) &&
        are_matrices_equal(C_seq, C_diag, N)) {
        printf("equal\n");
    } else {
        printf("not equal!!!\n");
    }

    free(A);
    free(B);
    free(B_T);
    free(B_T_d);
    free(C_seq);
    free(C_avx);
    free(C_avx_d);
    free(C_diag);
}


int main() {
    srand(time(NULL));

    int sizes[] = {8, 512, 1024, 2048}; 
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    for (int i = 0; i < num_sizes; i++) {
        int N = sizes[i];
        benchmark(N);
    }

    return 0;
}
