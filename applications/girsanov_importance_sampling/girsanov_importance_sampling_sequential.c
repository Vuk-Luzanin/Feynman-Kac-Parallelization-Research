#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ---------------- PARAMETERS ---------------- */
static const double x0    = 0.0;
static const double T     = 1.0;
static const double sigma = 2.0;
static const double K     = 6.0;

/* ---------------- RANDOM NORMAL ---------------- */
/* Box–Muller: rand() is uniform, not Gaussian */
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
        u = 2.0 * rand() / RAND_MAX - 1.0;
        v = 2.0 * rand() / RAND_MAX - 1.0;
        s = u*u + v*v;
    } while (s >= 1.0 || s == 0.0);

    s = sqrt(-2.0 * log(s) / s);
    spare = v * s;
    return u * s;
}

/* ---------------- DRIFT FUNCTIONS ---------------- */
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

/* ---------------- INDICATOR FUNCTION ---------------- */
int f_indicator(double *x, int N)
{
    double max = x[0];
    for (int i = 1; i <= N; ++i)
        if (x[i] > max) max = x[i];

    return (max >= K);
}

/* ---------------- IMPORTANCE SAMPLING CORE ---------------- */
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

    /* Monte Carlo loop (PARALLELIZABLE) */
    for (int m = 0; m < M; ++m) {

        double X[N + 1];
        X[0] = x0;

        /* Simulate path */
        for (int n = 0; n < N; ++n) {
            X[n + 1] = X[n]
                     + b_fun(X[n]) * dt
                     + sigma * sqrt_dt * randn();
        }

        /* Likelihood ratio */
        double logM = 0.0;
        for (int n = 0; n < N; ++n) {
            double b = b_fun(X[n]);
            logM += b * (X[n + 1] - X[n])
                  - 0.5 * b * b * dt;
        }
        double M_weight = exp(-(1.0 / (sigma * sigma)) * logM);

        /* Indicator */
        double value = f_indicator(X, N) * M_weight;

        sum += value;
        sum_sq += value * value;
    }

    *mean = sum / M;
    *var  = sum_sq / M - (*mean) * (*mean);
}

/* ---------------- MAIN ---------------- */
int main(void)
{
    srand(42);

    int M = 10000;
    int Ns[] = {10, 100, 1000, 10000};

    printf("Exact probability (reflection principle):\n");
    double exact = 2.0 * (1.0 - 0.5 *
        (1.0 + erf((K - x0) / (sigma * sqrt(2.0 * T)))));
    printf("P = %.6f\n\n", exact);

    for (int i = 0; i < 4; ++i) {
        int N = Ns[i];

        double mean, var;

        printf("N = %d\n", N);

        /* Without importance sampling */
        importance_sampling(drift_zero, M, N, &mean, &var);
        printf("  No IS : mean = %.6f, var = %.6e\n", mean, var / M);

        /* With importance sampling */
        importance_sampling(drift_importance, M, N, &mean, &var);
        printf("  IS    : mean = %.6f, var = %.6e\n\n", mean, var / M);
    }

    return 0;
}
