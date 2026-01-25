#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>    // for time measurement


/* 
Indicator function 1_{[-1,1]}
 */
static inline int initial_condition(double x)
{
    return (x >= -1.0 && x <= 1.0);
}

/* 
Standard normal CDF via erf
Φ(z) = 0.5 * (1 + erf(z / sqrt(2)))
 */
static inline double normal_cdf(double z)
{
    return 0.5 * (1.0 + erf(z / sqrt(2.0)));
}

/* 
Exact solution of heat equation
via Feynman–Kac
 */
double exact_solution(double t, double x)
{
    if (t <= 0.0)
        return initial_condition(x);

    double denom = sqrt(2.0 * t);
    double a = (1.0 + x) / denom;
    double b = (-1.0 + x) / denom;

    return normal_cdf(a) - normal_cdf(b);
}

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


int main(void)
{
    /* Parameters */
    const int    N    = 100;     /* time steps, number of steps until the end of movement */
    const int    M    = 1000;    /* Monte Carlo paths from one point */
    const int    n_mc = 20;      /* spatial points, number of start points */
    const double T    = 1.0;
    const double L    = 3.0; 

    const double dt      = T / N;
    const double sqrt_dt = sqrt(dt);
    const double sqrt_2  = sqrt(2.0);

    printf("TEST: heat equation arguments [%d] and sequential\n", N);
    double wtime = omp_get_wtime();


    /* Spatial grid */
    double x_mc[n_mc];
    for (int i = 0; i < n_mc; ++i)
        x_mc[i] = -L + 2.0 * L * i / (n_mc - 1);


    /* MC estimator: mc[i][n] ≈ u(t_n, x_i) */
    double mc_estimator[n_mc][N + 1];
    for (int i = 0; i < n_mc; ++i)
        for (int n = 0; n <= N; ++n)
            mc_estimator[i][n] = 0.0;

    srand(42);

    // MONTE CARLO SIMULATION
    for (int m = 0; m < M; ++m) {

        double W = 0.0;  /* Brownian motion */

        for (int n = 0; n <= N; ++n) {

            if (n > 0)
                W += sqrt_dt * randn();

            for (int i = 0; i < n_mc; ++i) {
                double val = x_mc[i] + sqrt_2 * W;
                mc_estimator[i][n] += initial_condition(val);
            }
        }
    }

    /* Normalize */
    for (int i = 0; i < n_mc; ++i)
        for (int n = 0; n <= N; ++n)
            mc_estimator[i][n] /= M;

    // OUTPUT: exact vs MC at t = T
    // printf("\nComparison at final time t = %.2f\n", T);
    // printf("   x        exact        MC\n");
    // printf("--------------------------------\n");

    double err = 0.0;
    for (int i = 0; i < n_mc; ++i) {
        double exact = exact_solution(T, x_mc[i]);
        double mcval = mc_estimator[i][N];
        err += (exact - mcval) * (exact - mcval);
        // printf("%+7.3f   %.6f   %.6f\n", x_mc[i], exact, mcval);
    }
    wtime = omp_get_wtime() - wtime;
    printf("%d    %lf    %lf\n", N, err, wtime);
    printf("TEST END\n");

    return 0; 
}
