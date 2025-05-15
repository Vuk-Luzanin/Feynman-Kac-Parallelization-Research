# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>
#include <omp.h>
#include "util.h"

#define DIMENSIONS 1

// potencijalna energiju u nekom centralno simetricnom polju unutar elipsoida - ovde je to 1 dimenzija
// energija veca sto je tacka dalja od centra i sto su dimenzije elipsoida manje - ovde je to 1 dimenzija
double potential ( double a, double x )
{
  double value;
  value = 2.0 * pow ( x / a / a, 2 ) + 1.0 / a / a;
  return value;
}
// generator pseudoslučajnih brojeva po uniformnoj raspodeli - svaka nit ima svoj seed, jer kada bi bio shared, uniformnost ne bi bila garantovana
// real 8-byte number in [0,1)
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

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double feynman_1(const double a, const double h, const double stepsz, const int ni, const int N) 
{
  int seed = 123456789;
  double err = 0.0;
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)


#pragma omp parallel default(none) shared(a, h, stepsz, ni, N) \
                                   firstprivate(seed) \
                                   reduction(+ : err) \
                                   reduction(+ : n_inside)
{
  // seed is private variable, so numbers can generate uniformly
  seed += omp_get_thread_num();

#pragma omp for schedule(static, 1)
  for (int i = 1; i <= ni; i++)
  {
    // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
    double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
    double chk = pow(x / a, 2);

    double w_exact = 0;
    double wt = 0;

    if (1.0 < chk)
    {
      // tacka nije unutar 1-D elipsoida
      w_exact = 1.0;
      wt = 1.0;
      continue;
    }
    // tacka je unutar 1-D elipsoida
    n_inside++;

    // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
    w_exact = exp(pow(x / a, 2) - 1.0);
    wt = 0.0;

    // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
    for (int trial = 0; trial < N; trial++)
    {
      double x1 = x;

      double w = 1.0;
      chk = 0.0;

      // kretanje cestice - dok se nalazi unutar elipsoida
      while (chk < 1.0)
      {
        // da li se pomeramo za +stepsz ili -stepsz
        double us;
        double dx;

        us = r8_uniform_01(&seed) - 0.5;
        if (us < 0.0)
        {
          dx = -stepsz;
        }
        else
        {
          dx = stepsz;
        }
        
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
      wt = wt + w;
    }
    // srednja vrenost tezine za N pokusaja
    wt = wt / (double)(N);

    // kvadrat razlike tacne i numericki dobijene vrednosti
    err += pow(w_exact - wt, 2);

#ifdef DEBUG
    printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
            x, y, z, wt, w_exact, fabs(w_exact - wt), steps / N);
#endif
  
  }
} // parallel
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}


double (*FUNCS[])(const double, const double, const double, const int, const int) = {feynman_1};

int main ( int argc, char **argv )
{
  const double a = 2.0;
  const double h = 0.0001;
  const int ni = 11;

  const double stepsz = sqrt(DIMENSIONS * h);

  if (argc < 3)
  {
    printf("Invalid number of arguments passed.\n");
    return 1;
  }

  // index of function
  const int func = atoi(argv[1]);
  
  // numer of walks per point
  const int N = atoi(argv[2]);

  printf("TEST: func=%d, N=%d, num_threads=%ld\n", func, N, get_num_threads());
  double wtime = omp_get_wtime();
  double err = FUNCS[func](a, h, stepsz, ni, N);
  wtime = omp_get_wtime() - wtime;
  printf("%d    %lf    %lf\n", N, err, wtime);
  printf("TEST END\n");

  return 0;
}


// --------------------------------------------------------------------
/*
  Purpose:

    MAIN is the main program for FEYNMAN_KAC_2D.

  Discussion:

    This program is derived from section 2.5, exercise 2.2 of Petersen and Arbenz.

    The problem is to determine the solution U(X,Y) of the following 
    partial differential equation:

      (1/2) Laplacian U - V(X,Y) * U = 0,

    inside the elliptic domain D:
 
      D = { (X,Y) | (X/A)^2+(Y/B)^2 <= 1 }
   
    with the boundary condition U(boundary(D)) = 1.

    The V(X,Y) is the potential function:

      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 ) + 1/A^2 + 1/B^2.

    The analytic solution of this problem is already known:

      U(X,Y) = exp ( (X/A)^2 + (Y/B)^2 - 1 ).

    Our method is via the Feynman-Kac Formula.

    The idea is to start from any (x,y) in D, and
    compute (x+Wx(t),y+Wy(t)) where 2D Brownian motion
    (Wx,Wy) is updated each step by sqrt(h)*(z1,z2),
    each z1,z2 are independent approximately Gaussian 
    random variables with zero mean and variance 1. 

    Each (x1(t),x2(t)) is advanced until (x1,x2) exits 
    the domain D.  

    Upon its first exit from D, the sample path (x1,x2) is stopped and a 
    new sample path at (x,y) is started until N such paths are completed.
 
    The Feynman-Kac formula gives the solution here as

      U(X,Y) = (1/N) sum(1 <= I <= N) Y(tau_i),

    where

      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s)) ds),

    and tau = first exit time for path (x1,x2). 

    The integration procedure is a second order weak accurate method:

      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
                [ x2(t) + sqrt ( h ) * z2 ]

    Here Z1, Z2 are approximately normal univariate Gaussians. 

    An Euler predictor approximates Y at the end of the step

      Y_e     = (1 - h*v(X(t)) * Y(t), 

    A trapezoidal rule completes the step:

      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].

  Licensing:

    This code is distributed under the MIT license. 

  Modified:

    31 May 2012

  Author:

    Original C 3D version by Wesley Petersen.
    C 2D version by John Burkardt.

  Reference:

    Peter Arbenz, Wesley Petersen,
    Introduction to Parallel Computing:
    A Practical Guide with Examples in C,
    Oxford, 2004,
    ISBN: 0-19-851577-4,
    LC: QA76.59.P47.
*/