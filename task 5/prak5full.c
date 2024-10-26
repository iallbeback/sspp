#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void run_simulation(int a, int b, int x, int N, double p, int P, double *prob_b, double *avg_time, double *elapsed_time) {
    int hits_b = 0;
    long total_steps = 0;
    double total_time = 0.0;

    double start_time = omp_get_wtime();

    #pragma omp parallel num_threads(P)
    {
        struct drand48_data state; 
        long seed = time(NULL) ^ omp_get_thread_num(); 
        srand48_r(seed, &state);  

        #pragma omp for reduction(+:hits_b, total_steps, total_time)
        for (int i = 0; i < N; i++) {
            int x_pos = x;
            double particle_start_time = omp_get_wtime();
            int steps = 0;

            while (x_pos > a && x_pos < b) {
                double r;
                drand48_r(&state, &r);
                if (r < p) {
                    x_pos++;
                } else {
                    x_pos--;
                }
                steps++;
            }

            double particle_end_time = omp_get_wtime();
            total_steps += steps;
            total_time += (particle_end_time - particle_start_time);

            if (x_pos == b) {
                hits_b++;
            }
        }
    }

    double end_time = omp_get_wtime();
    *elapsed_time = end_time - start_time;
    *prob_b = (double)hits_b / N;
    *avg_time = total_time / N;
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s a b x N p P\n", argv[0]);
        return 1;
    }

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int x = atoi(argv[3]);
    int N = atoi(argv[4]);
    double p = atof(argv[5]);
    int P = atoi(argv[6]);

    int test_N[] = {100, 1000, 10000, 100000, 1000000};
    int test_P[] = {1, 2, 4, 8, 16};
    int num_tests_N = sizeof(test_N) / sizeof(test_N[0]);
    int num_tests_P = sizeof(test_P) / sizeof(test_P[0]);

    FILE *file = fopen("result.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    fprintf(file, "Testing with varying N, fixed P = %d\n", P);
    fprintf(file, "%-10s %-20s %-25s %-25s\n", "N", "Prob(b)", "Avg Lifetime (s)", "Elapsed Time (s)");
    for (int i = 0; i < num_tests_N; i++) {
        double prob_b, avg_time, elapsed_time;
        run_simulation(a, b, x, test_N[i], p, P, &prob_b, &avg_time, &elapsed_time);
        fprintf(file, "%-10d %-20.6f %-25.6f %-25.6f\n", test_N[i], prob_b, avg_time, elapsed_time);
        printf("%-10d %-20.6f %-25.6f %-25.6f\n", test_N[i], prob_b, avg_time, elapsed_time);
    }

    fprintf(file, "\nTesting with varying P, fixed N = %d\n", N);
    fprintf(file, "%-10s %-20s %-25s %-25s\n", "P", "Prob(b)", "Avg Lifetime (s)", "Elapsed Time (s)");
    for (int i = 0; i < num_tests_P; i++) {
        double prob_b, avg_time, elapsed_time;
        run_simulation(a, b, x, N, p, test_P[i], &prob_b, &avg_time, &elapsed_time);
        fprintf(file, "%-10d %-20.6f %-25.6f %-25.6f\n", test_P[i], prob_b, avg_time, elapsed_time);
        printf("%-10d %-20.6f %-25.6f %-25.6f\n", test_P[i], prob_b, avg_time, elapsed_time);
    }

    fclose(file);
    printf("\nResults have been saved to result.txt\n");

    return 0;
}
