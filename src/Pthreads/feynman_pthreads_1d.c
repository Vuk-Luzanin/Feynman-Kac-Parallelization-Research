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
static double stepsz;
static int num_threads = 8;

static double wt[NI+1] = {0};
static double w_exact[NI+1] = {0};

typedef struct {
    int i;       // to add on seed
    double x;   // coordinate
} grid_point_t;

typedef struct {
    int thread_id;
    int N;
    int n_inside;
    // pamti sve tacke koje su unutar elipsoida - samo za njih prolazimo algoritam
    grid_point_t inside_points[MAX_INSIDE];
    double wt_local[NI+1];  // lokalni mutex za svaku nit (svaka nit rezultat akumulira u svoju lokalnu promenljivu)
} thread_arg_t;

double potential(double a, double x) {
    return 2.0 * pow(x / a / a, 2) + 1.0 / a / a;
}

double r8_uniform_01(int *seed) {
    int k = *seed / 127773;
    *seed = 16807 * (*seed - k * 127773) - k * 2836;
    if (*seed < 0) *seed += 2147483647;
    return (double)(*seed) * 4.656612875E-10;
}

void* worker(void *arg_) {
    thread_arg_t *arg = (thread_arg_t *)arg_;
    int tid = arg->thread_id;
    int trials = arg->N;
    int n_inside = arg->n_inside;

    // inicijalizujemo lokalni niz gresaka
    for (int i=0; i<=NI; i++)
        arg->wt_local[i] = 0.0;

    for (int i_idx = 0; i_idx < n_inside; ++i_idx) {
        int i = arg->inside_points[i_idx].i;
        double x0 = arg->inside_points[i_idx].x;

        // suma w za sve pokusaje koje radi jedna nit iz jedne tacke
        double local_sum = 0.0;

        // lokalno, svaka tacka racuna svoje pokusaje - pomeren brojac za num_threads pri svakoj iteraciji - bolje za neravnomeran broj iteracija
        // svaka nit radi N / num_threads iteracija
        for (int t = tid; t < trials; t += num_threads) {
            int seed = 123456789u + i * trials + t;
            double x1 = x0;
            double w = 1.0, chk = 0.0;

            while (chk < 1.0) {
                double dx = (r8_uniform_01(&seed) - 0.5 < 0) ? -stepsz : stepsz;
                double vs = potential(a, x1);
                x1 += dx;
                double vh = potential(a, x1);
                double we = (1.0 - h * vs) * w;
                w = w - 0.5 * h * (vh * we + vs * w);
                chk = pow(x1 / a, 2);
            }

            local_sum += w;
        }

        // dodajemo na lokalnu promenljivu
        arg->wt_local[i] += local_sum;
    }

    return NULL;
}

double feynman_pthreads(const double a, const int N) {
    int n_inside = 0;
    grid_point_t inside_points[MAX_INSIDE];     // necemo iskoristiti ceo niz (vec samo koliko ima tacaka unutar elipse)

    // pronalazimo tacke unutar elipsoida
    for (int i = 1; i <= NI; ++i) {
        double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
        double chk = pow(x / a, 2);
        w_exact[i] = 0.0;
        wt[i] = 0.0;

        if (chk >= 1.0) continue;

        w_exact[i] = exp(pow(x / a, 2) - 1.0);
        inside_points[n_inside].i = i;
        inside_points[n_inside].x = x;
        n_inside++;
    }

    pthread_t threads[num_threads];
    thread_arg_t args[num_threads];

    for (int t = 0; t < num_threads; ++t) {
        args[t].thread_id = t;
        args[t].N = N;
        args[t].n_inside = n_inside;
        for (int j = 0; j < n_inside; ++j)
            args[t].inside_points[j] = inside_points[j];    // svaka nit dobije sve tacke

        pthread_create(&threads[t], NULL, worker, &args[t]);
    }

    for (int t = 0; t < num_threads; ++t) {
        pthread_join(threads[t], NULL);
    }

    // redukujemo sve lokalne rezultate
    for (int t=0; t<num_threads; t++)
        for (int i=0; i<=NI; i++)
            wt[i] += args[t].wt_local[i];

    double err = 0.0;
    for (int i = 1; i <= NI; ++i)
        if (w_exact[i] != 0.0)
            err += pow(w_exact[i] - (wt[i] / (double)(N)), 2);

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
    double wtime = omp_get_wtime();
    double err = feynman_pthreads(a, N);
    wtime = omp_get_wtime() - wtime;
    printf("%d    %lf    %lf\n", N, err, wtime);
    printf("TEST END\n");

    return 0;
}
