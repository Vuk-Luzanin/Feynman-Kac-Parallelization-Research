#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "util.h"

/*  PARAMETERS  */
static const double x0    = 0.0;
static const double T     = 1.0;
static const double sigma = 2.0;
static const double K     = 6.0;

/*
Standard normal RNG (Box–Muller)
randn(): generates N(0,1) needed for Brownian motion;
plain rand() gives uniform numbers and is not sufficient 
*/
double randn(void)
{
    static int has_spare = 0;
    static double spare;

    if (has_spare) {
        has_spare = 0;
        return spare;
    }

    has_spare = 1;
    double u, v, s;
    do {
        u = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
        v = (rand() / (double)RAND_MAX) * 2.0 - 1.0;
        s = u*u + v*v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

// thread-safe normal RNG
static inline double randn_r(unsigned int *seed)
{
    double u1 = rand_r(seed) / (double)RAND_MAX;
    double u2 = rand_r(seed) / (double)RAND_MAX;

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;

    return r * cos(theta);
}

/*  DRIFT FUNCTIONS  */
double drift_zero(double x)
{
    (void)x;
    return 0.0;
}

double drift_importance(double x)
{
    double b = K / T;
    return (x < K) ? b : 0.0;
}

/*  INDICATOR FUNCTION  */
int f_indicator(double *x, int N)
{
    double max = x[0];
    for (int i = 1; i <= N; ++i)
        if (x[i] > max) max = x[i];

    return (max >= K);
}

/*  IMPORTANCE SAMPLING  */
void importance_sampling(
    double (*b_fun)(double),
    int M,
    int N,
    double *mean,
    double *var
)
{
    double dt = T / N;
    double sqrt_dt = sqrt(dt);

    double sum = 0.0;
    double sum_sq = 0.0;

#pragma omp parallel default(none) \
    shared(M, N, dt, sqrt_dt, b_fun) \
    reduction(+:sum, sum_sq)
{
    unsigned int seed = 1234 + 1000 * omp_get_thread_num();

#pragma omp for schedule(static)
    for (int m = 0; m < M; ++m) {

        double X[N + 1];
        X[0] = x0;

        /* Simulate path */
        for (int n = 0; n < N; ++n) {
            X[n + 1] = X[n]
                     + b_fun(X[n]) * dt
                     + sigma * sqrt_dt * randn_r(&seed);
        }

        /* Likelihood ratio */
        double logM = 0.0;
        for (int n = 0; n < N; ++n) {
            double b = b_fun(X[n]);
            logM += b * (X[n + 1] - X[n])
                  - 0.5 * b * b * dt;
        }

        double M_weight = exp(-(1.0 / (sigma * sigma)) * logM);
        double value = f_indicator(X, N) * M_weight;

        sum     += value;
        sum_sq += value * value;
    }
} // parallel

    *mean = sum / M;
    *var  = sum_sq / M - (*mean) * (*mean);
}

// void print_confidence_interval(double mean, double var, int M) {
//     const double z = 1.96;  // 95% confidence
//     double std_err = sqrt(var / M);
//     printf("         95%% CI: [%.6f, %.6f]\n",
//            mean - z * std_err,
//            mean + z * std_err);
// }


int main(int argc, char **argv)
{
    srand(42);

    int M = 10000;
    // int Ns[] = {10, 100, 1000};   // steps for time discretization - now sent trough command line

    if (argc < 2)
    {
        printf("Invalid number of arguments passed.\n");
        return 1;
    }

    // numer of steps for time discretization
    const int N = atoi(argv[1]);

    // printf("Exact probability (reflection principle):\n");
    // double exact = 2.0 * (1.0 - 0.5 *
    //     (1.0 + erf((K - x0) / (sigma * sqrt(2.0 * T)))));
    // printf("P = %.6f\n\n", exact);

    printf("TEST: girsanov importance sampling arguments [%d] and sequential\n", N);
    double wtime = omp_get_wtime();
    double mean, var;

    // printf("N = %d\n", N);

    /* Without importance sampling */
    // importance_sampling(drift_zero, M, N, &mean, &var);
    // printf("  No IS : mean = %.6f, var = %.6e\n", mean, var / M);
    // print_confidence_interval(mean, var, M);

    /* With importance sampling */
    importance_sampling(drift_importance, M, N, &mean, &var);
    // printf("  IS    : mean = %.6f, var = %.6e\n\n", mean, var / M);
    // print_confidence_interval(mean, var, M);

    wtime = omp_get_wtime() - wtime;
    printf("%d    %lf    %lf\n", N, var, wtime);
    printf("TEST END\n");

    return 0;
}
