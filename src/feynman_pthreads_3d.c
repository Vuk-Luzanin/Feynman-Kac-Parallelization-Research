#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <omp.h> // used only for time measurement and passing number of threads -> to uniform the run.py script
#include "util.h"

#define NI 16
#define NJ 11
#define NK 6

int num_threads = 8;
static int ni = NI;
static int nj = NJ;
static int nk = NK;

static double wt[NI+1][NJ+1][NK+1] = {{{0}}};
static double w_exact[NI+1][NJ+1][NK+1] = {{{0}}};

static double a = 3.0;
static double b = 2.0;
static double c = 1.0;
static double h = 0.001;

static double stepsz;

typedef struct 
{
    int i, j, k;
    double x0, y0, z0;
    int start_trial;
    int end_trial;
} trial_arg_t;

static pthread_mutex_t wt_mutex;

double potential(double a, double b, double c, double x, double y, double z)
{
    return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) +
           1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

double r8_uniform_01(int *seed)
{
    int k = *seed / 127773;
    *seed = 16807 * (*seed - k * 127773) - k * 2836;
    if (*seed < 0)
        *seed += 2147483647;
    return (double)(*seed) * 4.656612875E-10;
}

void* trial_worker(void *varg)
{
    trial_arg_t *arg = (trial_arg_t*) varg;

    for (int trial_id = arg->start_trial; trial_id < arg->end_trial; trial_id++) {
        int seed = 123456789u + trial_id;
        double x1 = arg->x0;
        double x2 = arg->y0;
        double x3 = arg->z0;

        double w = 1.0;
        double chk = 0.0;

        while (chk < 1.0) {
            double ut = r8_uniform_01(&seed);
            double us, dx = 0, dy = 0, dz = 0;

            if (ut < 1.0 / 3.0) 
            {
                us = r8_uniform_01(&seed) - 0.5;
                dx = (us < 0.0) ? -stepsz : stepsz;
            }

            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0) 
            {
                us = r8_uniform_01(&seed) - 0.5;
                dy = (us < 0.0) ? -stepsz : stepsz;
            }

            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0) 
            {
                us = r8_uniform_01(&seed) - 0.5;
                dz = (us < 0.0) ? -stepsz : stepsz;
            }

            double vs = potential(a, b, c, x1, x2, x3);
            x1 += dx; x2 += dy; x3 += dz;
            double vh = potential(a, b, c, x1, x2, x3);
            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);
            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
        }

        pthread_mutex_lock(&wt_mutex);
        wt[arg->i][arg->j][arg->k] += w;
        pthread_mutex_unlock(&wt_mutex);
    }

    return NULL;
}

double feynman_pthreads(double a, double b, double c, int ni, int nj, int nk, int N)
{
    int n_inside = 0;

    for (int i = 1; i <= ni; i++) 
    {
        for (int j = 1; j <= nj; j++) 
        {
            for (int k = 1; k <= nk; k++) 
            {
                double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
                double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
                double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);
                double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);
                w_exact[i][j][k] = 0.0;
                wt[i][j][k] = 0.0;

                if (1.0 < chk)
                    continue;

                n_inside++;
                w_exact[i][j][k] = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);

                pthread_t threads[num_threads];
                trial_arg_t args[num_threads];
                int trials_per_thread = N / num_threads;
                int remainder = N % num_threads;
                int current = 0;

                for (int t = 0; t < num_threads; t++) 
                {
                    int start = current;
                    int count = trials_per_thread + (t < remainder ? 1 : 0);
                    int end = start + count;
                    current = end;

                    args[t].i = i;
                    args[t].j = j;
                    args[t].k = k;
                    args[t].x0 = x;
                    args[t].y0 = y;
                    args[t].z0 = z;
                    args[t].start_trial = start;
                    args[t].end_trial = end;

                    pthread_create(&threads[t], NULL, trial_worker, &args[t]);
                }

                for (int t = 0; t < num_threads; t++) 
                {
                    pthread_join(threads[t], NULL);
                }
            }
        }
    }

    double err = 0.0;
    for (int i = 0; i <= NI; ++i)
        for (int j = 0; j <= NJ; ++j)
            for (int k = 0; k <= NK; ++k)
                if (w_exact[i][j][k] != 0.0)
                    err += pow(w_exact[i][j][k] - (wt[i][j][k] / (double)(N)), 2);

    return sqrt(err / (double)(n_inside));
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Invalid number of arguments passed.\n");
        return 1;
    }

    const int N = atoi(argv[1]);
    num_threads = get_num_threads();

    stepsz = sqrt(3 * h);
    pthread_mutex_init(&wt_mutex, NULL);

    printf("TEST: N=%d, num_threads=%d\n", N, num_threads);
    double wtime = omp_get_wtime();
    double err = feynman_pthreads(a, b, c, ni, nj, nk, N);
    wtime = omp_get_wtime() - wtime;
    printf("%d    %lf    %lf\n", N, err, wtime);
    printf("TEST END\n");

    pthread_mutex_destroy(&wt_mutex);

    return 0;
}
