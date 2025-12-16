#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <omp.h>
#include "util.h"

#define DIMENSIONS  1
#define NI          11
#define MAX_INSIDE  NI

static double a = 2.0;
static double h = 0.0001;
static double stepsz = 0.0;
static int num_threads = 8;

static double wt[NI+1] = {0};
static double w_exact[NI+1] = {0};

static int global_trial_index;  // globalni atomic brojac za dinamicku raspodelu posla


typedef struct {
    int i;       // to add on seed
    double x;   // coordinate
} grid_point_t;

typedef struct {
    int thread_id;
    int N;
    int n_inside;
    // pamti sve tacke koje su unutar elipsoida - samo za njih prolazimo algoritam
    double wt_local[NI+1];  // lokalni mutex za svaku nit (svaka nit rezultat akumulira u svoju lokalnu promenljivu)
} thread_arg_t;

static grid_point_t inside_points[MAX_INSIDE];


inline double potential(double a, double x) {
    double a2_inv = 1.0 / (a * a);
    return 2.0 * (x * x) * a2_inv * a2_inv + a2_inv;
}

inline double r8_uniform_01(int *seed) {
    int k = *seed / 127773;
    *seed = 16807 * (*seed - k * 127773) - k * 2836;
    if (*seed < 0) *seed += 2147483647;
    return (double)(*seed) * 4.656612875E-10;
}

// restrict - preko ovog pokazivaca se jedino pristupa ovoj lokaciji
// dynamic work distribution 
void* worker_dynamic(void* restrict arg_) {
    thread_arg_t *arg = (thread_arg_t *)arg_;
    int trials = arg->N;
    int n_inside = arg->n_inside;
    int total_trials = n_inside * trials;   // Ukupno pokušaja = broj_tačaka * N

    // inicijalizujemo lokalni niz gresaka
    for (int i = 0; i <= NI; i++) 
        arg->wt_local[i] = 0.0;

    // lokalno, svaka tacka racuna svoje pokusaje - bolje za neravnomeran broj iteracija - dinamicko rasporedjivanje
    int t;
    while ((t = __atomic_fetch_add(&global_trial_index, 1, __ATOMIC_RELAXED)) < total_trials) {
        int grid_idx = t / trials;  // Određuje tačku (0..n_inside-1)
        int trial_idx = t % trials; // Redni broj pokušaja (0..trials-1)

        int i = inside_points[grid_idx].i;      // Globalni indeks tačke (1..NI)
        double x0 = inside_points[grid_idx].x;  // Koordinata tačke

        int seed = 123456789u + i*t + trial_idx;  // Seed zavisi od tačke i pokušaja
        double x1 = x0;
        double w = 1.0;

        while ((x1 / a) * (x1 / a) < 1.0) {
            double dx = (r8_uniform_01(&seed) < 0.5) ? -stepsz : stepsz;
            double vs = potential(a, x1);
            x1 += dx;
            double vh = potential(a, x1);
            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);
        }

        arg->wt_local[i] += w;  
    }

    return NULL;
}

// Worker with static work distribution via round-robin stepping 
void* worker_static(void* restrict arg_) {
    thread_arg_t *arg = (thread_arg_t *)arg_;        
    int trials      = arg->N;
    int n_inside    = arg->n_inside;
    int tid         = arg->thread_id;
    int total_trials = n_inside * trials;   // Ukupno pokušaja = broj_tačaka * N

    // initialize local accumulator
    for (int i = 0; i <= NI; i++) {
        arg->wt_local[i] = 0.0;
    }

    // static distribution: each thread processes t = tid, tid + n_threads, ...
    for (int t = tid; t < total_trials; t += num_threads) {
        int grid_idx = t / trials;
        int trial_idx= t % trials;

        int i = inside_points[grid_idx].i;
        double x0 = inside_points[grid_idx].x;

        int seed = 123456789u + i * t + trial_idx;
        double x1 = x0;
        double w  = 1.0;

        // walk until exit
        while ((x1 / a) * (x1 / a) < 1.0) {
            double dx = (r8_uniform_01(&seed) < 0.5) ? -stepsz : stepsz;
            double vs = potential(a, x1);
            x1 += dx;
            double vh = potential(a, x1);
            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);
        }

        // accumulate result
        arg->wt_local[i] += w;
    }

    return NULL;
}


double feynman_pthreads(const double a, const int N) {
    int n_inside = 0;

    // pronalazimo tacke unutar elipsoida
    for (int i = 1; i <= NI; ++i) {
        // double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);

        // uzima indekse od krajeva ka pocetku
        int alt_i = (i % 2 == 1) ? (i - 1) / 2 + 1 : NI - (i / 2) + 1;
        double x = ((double)(NI - alt_i) * (-a) + (double)(alt_i - 1) * a) / (double)(NI - 1);

        double r = x / a;
        double chk = r * r;
        w_exact[i] = 0.0;
        wt[i] = 0.0;

        if (chk > 1.0) continue;

        w_exact[i] = exp(r * r - 1.0);

        inside_points[n_inside].i = i;
        inside_points[n_inside].x = x;
        n_inside++;
    }

    pthread_t threads[num_threads];
    thread_arg_t args[num_threads];

    global_trial_index = 0;  // resetujemo pre pokretanja niti

    for (int t = 0; t < num_threads; ++t) {
        args[t].N = N;
        args[t].n_inside = n_inside;
        args[t].thread_id = t;
        
        pthread_create(&threads[t], NULL, worker_dynamic, &args[t]);
    }

    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], NULL);
    }


    for (int i = 0; i <= NI; i++) {
        for (int t = 0; t < num_threads; t++) {
            wt[i] += args[t].wt_local[i];
        }
    }

    double err = 0.0;
    for (int i = 1; i <= NI; ++i)
        if (w_exact[i] != 0.0)
        {
            double diff = w_exact[i] - (wt[i] / (double)(N));
            err += diff * diff;
        }

    return sqrt(err / (double)(n_inside));
}

int main(int argc, char **argv) {
    if (argc < 2) {
        printf("Usage: %s <N>\n", argv[0]);
        return 1;
    }

    const int N = atoi(argv[1]);

    num_threads = get_num_threads();
    stepsz = sqrt(DIMENSIONS * h);

    printf("TEST: N=%d, num_threads=%d\n", N, num_threads);

    double wtime = omp_get_wtime();  // početno vreme
    double err = feynman_pthreads(a, N);
    wtime = omp_get_wtime() - wtime;

    printf("%d    %lf    %lf\n", N, err, wtime);
    printf("TEST END\n");

    return 0;
}
