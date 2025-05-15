#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "util.h"

#define NUM_LOCKS   256
#define DIMENSIONS  3
#define NI          16
#define NJ          11
#define NK          6

static double a = 3.0;
static double b = 2.0;
static double c = 1.0;
static double h = 0.001;

static double stepsz;


// potencijalna energiju u nekom centralno simetricnom polju unutar elipsoida
// energija veca sto je tacka dalje od centra i sto su dimenzije elipsoida manje
double potential(double a, double b, double c, double x, double y, double z)
{
  return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
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


double feynman_1(const double a, const double b, const double c, const double h, const double stepsz, const int N) 
{
  int seed = 123456789;
  double err = 0.0;
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

#pragma omp parallel default(none) shared(a, b, c, h, stepsz, N) \
                                   firstprivate(seed) \
                                   reduction(+ : err) \
                                   reduction(+ : n_inside)
{
  // seed is private variable, so numbers can generate uniformly
  seed += omp_get_thread_num();

#pragma omp for collapse(3)
  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++)
    {
      for (int k = 1; k <= NK; k++)
      {
        // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
        double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
        double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
        double z = ((double)(NK - k) * (-c) + (double)(k - 1) * c) / (double)(NK - 1);
        // check if point is inside elipsoid [a, b, c]
        double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);

        double w_exact = 0.0;
        double wt = 0.0;

        if (chk > 1.0)
        {
          // tacka nije unutar elipsoida
// deo oznacen za DEBUG se kompajlira ako se prilikom prevodjenja navede opcija -DDEBUG
#ifdef DEBUG
          printf("  %7.4f  %7.4f  %7.4f  %10.4e  %10.4e  %10.4e  %8d\n",
                 x, y, z, 1.0, 1.0, 0.0, 0);
#endif

          continue;
        }

        // tacka je unutar elipsoida
        n_inside++;

        // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
        w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);
        wt = 0.0;

#ifdef DEBUG
        int steps = 0;
#endif
        // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
        for (int trial = 0; trial < N; trial++)
        {
          double x1 = x;
          double x2 = y;
          double x3 = z;

          double w = 1.0;
          chk = 0.0;

          // kretanje cestice - dok se nalazi unutar elipsoida
          while (chk < 1.0)
          {
            // da li se pomeramo ili ne (po sve tri ose) - sa verovatnocom 1/3 za svaku osu
            double ut = r8_uniform_01(&seed);

            // da li se pomeramo za +stepsz ili -stepsz
            double us;

            double dx = 0;
            double dy = 0;
            double dz = 0;

            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dx = -stepsz;
              } 
              else
              {
                dx = stepsz;
              }
            }
            else
            {
              dx = 0.0;
            }
            
            ut = r8_uniform_01(&seed);

            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dy = -stepsz;
              }
              else
              {
                dy = stepsz;
              }
            }
            else
            {
              dy = 0.0;
            }
           
            ut = r8_uniform_01(&seed);

            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dz = -stepsz;
              }
              else
              {
                dz = stepsz;
              }
            }
            else
            {
              dz = 0.0;
            }
              
            // potential before moving
            double vs = potential(a, b, c, x1, x2, x3);

            // move
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

#ifdef DEBUG
            ++steps;
#endif

            // potential after moving
            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;           // Euler-ov korak
            w = w - 0.5 * h * (vh * we + vs * w);     // trapezna aproksimacija

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
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
    }
  }
} // parallel
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double feynman_2(const double a, const double b, const double c, const double h, const double stepsz, const int N) 
{
  static int seed = 123456789;  // set to be static -> to be global (shared) by default
  int n_inside = 0;

  // for every point in grid - must be initialized (moguce je da se desi da neki taskovi krenu da pristupaju matrici pre postavljanja na 0 -> citaju neinicijalizovanu memoriju)
  double w_exact[17][12][7] = {{{0}}};
  double wt[17][12][7] = {{{0}}};

#pragma omp threadprivate(seed) // then, it is private for every thread (no race condition)
#pragma omp parallel default(none) shared(a, b, c, h, stepsz, N, n_inside, w_exact, wt)
{
  seed += omp_get_thread_num();

#pragma omp single
{
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
        {
          continue;
        }

        n_inside++;

        w_exact[i][j][k] = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);

        for (int trial = 0; trial < N; trial++)
        {
#pragma omp task shared(wt)
{
          double x1 = x;
          double x2 = y;
          double x3 = z;

          double w = 1.0;
          chk = 0.0;

          while (chk < 1.0)
          {
            double ut = r8_uniform_01(&seed);
            double us;

            double dx;
            double dy;
            double dz;
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dx = -stepsz;
              }
              else
              {
                dx = stepsz;
              }
            }
            else
            {
              dx = 0.0;
            }

              ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dy = -stepsz;
              }
              else
              {
                dy = stepsz;
              }
            }
            else
            {
              dy = 0.0;
            }
            
            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dz = -stepsz;
              }
              else
              {
                dz = stepsz;
              }
            }
            else
            {
              dz = 0.0;
            }

            // potential before moving
            double vs = potential(a, b, c, x1, x2, x3);

            // move
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

            // potential after moving
            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }
// jer se menja u kontekstu vise taskova
#pragma omp atomic
          wt[i][j][k] += w;
} // task
        }
      }
    }
  }
} // single
} // parallel
  // na kraju obracunati gresku po svim osama - u matrici wt
  double err = 0.0;
  for (int i = 0; i <= 16; ++i)
  {
      for (int j = 0; j <= 11; ++j)
      {
          for (int k = 0; k <= 6; ++k)
          {
              if (w_exact[i][j][k] == 0.0)
              {
                  // kada tacka nije unutar elipsoida
                  continue;
              }
              err += pow(w_exact[i][j][k] - (wt[i][j][k] / (double)(N)), 2);
          }
      }
  }
  return sqrt(err / (double)(n_inside));
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 
omp_lock_t locks[NUM_LOCKS];

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

// using one lock for more threads
double feynman_3(const double a, const double b, const double c, const double h, const double stepsz, const int N) 
{
  static int seed = 123456789;  // set to be static -> to be global (shared) by default
  int n_inside = 0;

  // for every point in grid
  double w_exact[17][12][7] = {{{0}}};
  double wt[17][12][7] = {{{0}}};

  for (int i = 0; i < NUM_LOCKS; i++) 
  {
    omp_init_lock(&locks[i]);
  }

#pragma omp threadprivate(seed) // then, it is private for every thread (no race condition)
#pragma omp parallel default(none) shared(a, b, c, h, stepsz, N, n_inside, w_exact, wt, locks)
{
  seed += omp_get_thread_num();

#pragma omp single
{
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
        {
          continue;
        }

        n_inside++;

        w_exact[i][j][k] = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);

        for (int trial = 0; trial < N; trial++)
        {
#pragma omp task shared(wt)
{
          double x1 = x;
          double x2 = y;
          double x3 = z;

          double w = 1.0;
          chk = 0.0;

          while (chk < 1.0)
          {
            double ut = r8_uniform_01(&seed);
            double us;

            double dx;
            double dy;
            double dz;
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dx = -stepsz;
              }
              else
              {
                dx = stepsz;
              }
            }
            else
            {
              dx = 0.0;
            }

              ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dy = -stepsz;
              }
              else
              {
                dy = stepsz;
              }
            }
            else
            {
              dy = 0.0;
            }
            
            ut = r8_uniform_01(&seed);
            if (ut < 1.0 / 3.0)
            {
              us = r8_uniform_01(&seed) - 0.5;
              if (us < 0.0)
              {
                dz = -stepsz;
              }
              else
              {
                dz = stepsz;
              }
            }
            else
            {
              dz = 0.0;
            }

            // potential before moving
            double vs = potential(a, b, c, x1, x2, x3);

            // move
            x1 = x1 + dx;
            x2 = x2 + dy;
            x3 = x3 + dz;

            // potential after moving
            double vh = potential(a, b, c, x1, x2, x3);

            double we = (1.0 - h * vs) * w;
            w = w - 0.5 * h * (vh * we + vs * w);

            chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
          }

          // koriscenje lock-a
          int lock_id = get_lock_index(i, j, k);    // izracunaj index lock-a koji je potreban
          omp_set_lock(&locks[lock_id]);
          wt[i][j][k] += w;
         omp_unset_lock(&locks[lock_id]);
} // task
        }
      }
    }
  }
} // single
} // parallel -> implicitna barijera za taskove -> obrada nad matricom je zavrsena
  // na kraju obracunati gresku po svim osama - u matrici wt
  double err = 0.0;
  for (int i = 0; i <= 16; ++i)
  {
      for (int j = 0; j <= 11; ++j)
      {
          for (int k = 0; k <= 6; ++k)
          {
              if (w_exact[i][j][k] == 0.0)
              {
                  // kada tacka nije unutar elipsoida
                  continue;
              }
              err += pow(w_exact[i][j][k] - (wt[i][j][k] / (double)(N)), 2);
          }
      }
  }

  // oslobađanje brava
  for (int i = 0; i < NUM_LOCKS; i++) {
    omp_destroy_lock(&locks[i]);
  }

  return sqrt(err / (double)(n_inside));
}


double (*FUNCS[])(const double, const double, const double, const double, const double, const int) = {feynman_1, feynman_2, feynman_3};


int main(int argc, char **argv)
{
  if (argc < 3)
  {
    printf("Invalid number of arguments passed.\n");
    return 1;
  }

  // index of function
  const int func = atoi(argv[1]);
  
  // numer of walks per point
  const int N = atoi(argv[2]);

  stepsz = sqrt(DIMENSIONS * h);      // velicina koraka kojim se cestica pomera u svakom koraku

  printf("TEST: func=%d, N=%d, num_threads=%ld\n", func, N, get_num_threads());
  double wtime = omp_get_wtime();
  double err = FUNCS[func](a, b, c, h, stepsz, N);
  wtime = omp_get_wtime() - wtime;
  printf("%d    %lf    %lf\n", N, err, wtime);
  printf("TEST END\n");

  return 0;
}



// commands: 
// make
/*
./feyman 0 1000
./feyman 0 5000
./feyman 0 10000
./feyman 0 20000
 */

 /*
 a, b and represents elipsoid - potential energy is like springs in three directions with different stiffness
 (x * x) / (a * a) + (y * y) / (b * b) + (z * z) / (c * c) <= 1.0

 - cestica se krece od tacke do tacke - po diskretizovanoj mrezi
 - iz svake tacke resetke, cestica krece N puta

 - w - trenutna tezina cestice u simulaciji
 */

 /******************************************************************************/
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


