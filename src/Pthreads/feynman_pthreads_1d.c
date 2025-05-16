#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <omp.h> // used only for time measurement and passing number of threads -> to uniform the run.py script
#include "util.h"

#define NUM_LOCKS   256
#define DIMENSIONS  1
#define NI          11

static double a = 2.0;
static double h = 0.0001;

static double stepsz;

int num_threads = 8;

static double wt[NI+1] = {0};
static double w_exact[NI+1] = {0};

static pthread_mutex_t wt_mutexes[NUM_LOCKS];

typedef struct 
{
    int i;
    double x0;
    int start_trial;
    int end_trial;
} trial_arg_t;


double potential ( double a, double x )
{
  double value;
  value = 2.0 * pow ( x / a / a, 2 ) + 1.0 / a / a;
  return value;
}

double r8_uniform_01(int *seed)
{
  int k;
  double r;

  k = *seed / 127773;
  *seed = 16807 * (*seed - k * 127773) - k * 2836;

  if (*seed < 0)
  {
    *seed = *seed + 2147483647;
  }
  r = (double)(*seed) * 4.656612875E-10;

  return r;
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
unsigned int get_lock_index(int i) 
{
  unsigned int hash = (unsigned int)(
      i * 73856093 
  );
  return hash % NUM_LOCKS;
}


void* trial_worker(void *varg)
{
    trial_arg_t *arg = (trial_arg_t*) varg;

    for (int trial_id = arg->start_trial; trial_id < arg->end_trial; trial_id++) 
    {
        int seed = 123456789u + trial_id;
        double x1 = arg->x0;

        double w = 1.0;
        double chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0) 
        {
#ifdef SMALL_STEP
            double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
            double us, dx = 0;

            us = r8_uniform_01(&seed) - 0.5;
            if (us < 0.0)
            {
                dx = -stepsz;
            }
            else
            {
                dx = stepsz;
            }
#endif
            // potential before moving
            double vs = potential(a, x1);

            // move
            x1 = x1 + dx;

            // potential after moving
            double vh = potential(a, x1);

            double we = (1.0 - h * vs) * w;           // Euler-ov korak
            w = w - 0.5 * h * (vh * we + vs * w);     // trapezna aproksimacija

            chk = pow(x1 / a, 2);
        }

        int lock_id = get_lock_index(arg->i);
        pthread_mutex_lock(&wt_mutexes[lock_id]);
        wt[arg->i] += w;
        pthread_mutex_unlock(&wt_mutexes[lock_id]);
    }
    return NULL;
}

double feynman_pthreads_1d(const double a, const int N) 
{
    int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

    for (int i = 1; i <= NI; i++)
    {
        // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = NI -> x = a
        double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
        double chk = pow(x / a, 2);
        w_exact[i] = 0.0;
        wt[i] = 0.0;

        if (1.0 < chk)
        {
            continue;
        }

        // tacka je unutar 1-D elipsoida
        n_inside++;

        w_exact[i] = exp(pow(x / a, 2) - 1.0);

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
            args[t].x0 = x;
            args[t].start_trial = start;
            args[t].end_trial = end;

            pthread_create(&threads[t], NULL, trial_worker, &args[t]);
        }

        for (int t = 0; t < num_threads; t++) 
        {
            pthread_join(threads[t], NULL);
        }
    }
    double err = 0.0;
    for (int i = 0; i <= NI; ++i)
        if (w_exact[i] != 0.0)
            err += pow(w_exact[i] - (wt[i] / (double)(N)), 2);

  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}


int main ( int argc, char **argv )
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
    double err = feynman_pthreads_1d(a, N);
    wtime = omp_get_wtime() - wtime;
    printf("%d    %lf    %lf\n", N, err, wtime);
    printf("TEST END\n");

    destroy_locks();

    return 0;
}
