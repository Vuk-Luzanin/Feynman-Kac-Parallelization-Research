#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <omp.h> // used only for time measurement and passing number of threads -> to uniform the run.py script
#include "util.h"

#define NUM_LOCKS   256
#define DIMENSIONS  3
#define NI          16
#define NJ          11
#define NK          6

int num_threads = 8;

static double wt[NI+1][NJ+1][NK+1] = {{{0}}};
static double w_exact[NI+1][NJ+1][NK+1] = {{{0}}};

static pthread_mutex_t wt_mutexes[NUM_LOCKS];

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


void init_locks() {
    for (int i = 0; i < NUM_LOCKS; i++) {
        pthread_mutex_init(&wt_mutexes[i], NULL);
    }
}

void destroy_locks() {
    for (int i = 0; i < NUM_LOCKS; i++) {
        pthread_mutex_destroy(&wt_mutexes[i]);
    }   
}


// something like hash function that maps indexes (i, j and k) into index of lock that is used for that group of elements
// treba obratiti paznju na to sto se nece sve brave podjednako koristiti (brave za tacke van elipsoida nece biti koriscene)
unsigned int get_lock_index(int i, int j, int k) 
{
  unsigned int hash = (unsigned int)(
      i * 73856093 ^ 
      j * 19349663 ^ 
      k * 83492791
  );
  return hash % NUM_LOCKS;
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

        while (chk < 1.0) 
        {
#ifdef SMALL_STEP
            double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
            double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
            double dz = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
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
#endif
            double vs = potential(a, b, c, x1, x2, x3);
            x1 += dx; x2 += dy; x3 += dz;
            double vh = potential(a, b, c, x1, x2, x3);
            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);
            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
        }

        int lock_id = get_lock_index(arg->i, arg->j, arg->k);
        pthread_mutex_lock(&wt_mutexes[lock_id]);
        wt[arg->i][arg->j][arg->k] += w;
        pthread_mutex_unlock(&wt_mutexes[lock_id]);
    }
    return NULL;
}

double feynman_pthreads_3d(double a, double b, double c, int N)
{
    int n_inside = 0;

    for (int i = 1; i <= NI; i++) 
    {
        for (int j = 1; j <= NJ; j++) 
        {
            for (int k = 1; k <= NK; k++) 
            {
                double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
                double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
                double z = ((double)(NK - k) * (-c) + (double)(k - 1) * c) / (double)(NK - 1);
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

    stepsz = sqrt(DIMENSIONS * h);

    init_locks();

    printf("TEST: N=%d, num_threads=%d\n", N, num_threads);
    double wtime = omp_get_wtime();
    double err = feynman_pthreads_3d(a, b, c, N);
    wtime = omp_get_wtime() - wtime;
    printf("%d    %lf    %lf\n", N, err, wtime);
    printf("TEST END\n");

    destroy_locks();

    return 0;
}
