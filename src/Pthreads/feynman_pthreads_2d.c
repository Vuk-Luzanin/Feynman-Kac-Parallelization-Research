#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <omp.h>  // used only for time measurement and passing number of threads -> to uniform the run.py script
#include "util.h"

#define NUM_LOCKS   512
#define DIMENSIONS  2
#define NI          6
#define NJ          11

int num_threads = 8;

static double wt[NI+1][NJ+1] = {{0}};
static double w_exact[NI+1][NJ+1] = {{0}};

static pthread_mutex_t wt_mutexes[NUM_LOCKS];

static double a = 2.0;
static double b = 1.0;
static double h = 0.001;

static double stepsz;

typedef struct 
{
    int i, j;
    double x0, y0;
    int start_trial;
    int end_trial;
} trial_arg_t;


double potential ( double a, double b, double x, double y )
{
  double value;
  value = 2.0 * ( pow ( x / a / a, 2 ) + pow ( y / b / b, 2 ) ) + 1.0 / a / a + 1.0 / b / b;

  return value;
}

int i4_ceiling ( double x )
{
  int value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1;
  }

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

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// something like hash function that maps indexes (i, j and k) into index of lock that is used for that group of elements
// treba obratiti paznju na to sto se nece sve brave podjednako koristiti (brave za tacke van elipsoida nece biti koriscene)
unsigned int get_lock_index(int i, int j) 
{
  unsigned int hash = (unsigned int)(
    i * 73856093 ^ 
    j * 19349663
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

    double w = 1.0;
    double chk = 0.0;

    while (chk < 1.0) 
    {
#ifdef SMALL_STEP
      double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
      double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
      double ut = r8_uniform_01 ( &seed );

      double us;
      double dx = 0;
      double dy = 0;

      if ( ut < 1.0 / 2.0 )
      {
        us = r8_uniform_01 ( &seed ) - 0.5;
        if ( us < 0.0)
        {
          dx = - stepsz;
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
      if ( ut < 1.0 / 2.0 )
      {
        us = r8_uniform_01(&seed) - 0.5;
        if (us < 0.0)
        { 
          dy = - stepsz;
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
#endif
      // potential before moving
      double vs = potential(a, b, x1, x2);

      // move
      x1 = x1 + dx;
      x2 = x2 + dy;
      
      // potential after moving
      double vh = potential(a, b, x1, x2);

      double we = (1.0 - h * vs) * w;           // Euler-ov korak
      w = w - 0.5 * h * (vh * we + vs * w);     // trapezna aproksimacija

      chk = pow(x1 / a, 2) + pow(x2 / b, 2);
    }

    int lock_id = get_lock_index(arg->i, arg->j);
    pthread_mutex_lock(&wt_mutexes[lock_id]);
    wt[arg->i][arg->j] += w;
    pthread_mutex_unlock(&wt_mutexes[lock_id]);
  }

  return NULL;
}


double feynman_pthreads_2d(const double a, const double b, const int N, const int func) 
{
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++ )
    {
      // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

      w_exact[i][j] = 0.0;
      wt[i][j] = 0.0;

      if ( 1.0 < chk )
      {
        // tacka nije unutar 1-D elipsoida
        continue;
      }

      // tacka je unutar 2-D elipsoida
      n_inside++;
  
      // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
      w_exact[i][j] = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);


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
        args[t].x0 = x;
        args[t].y0 = y;
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

  double err = 0.0;
  for (int i = 0; i <= NI; ++i)
    for (int j = 0; j <= NJ; ++j)
      if (w_exact[i][j] != 0.0)
        err += pow(w_exact[i][j] - (wt[i][j] / (double)(N)), 2);

  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}


int main ( int argc, char **argv )
{
  if (argc < 3)
  {
    printf("Invalid number of arguments passed.\n");
    return 1;
  }

  const int func = atoi(argv[1]);
  const int N = atoi(argv[2]);
  num_threads = get_num_threads();
  
  stepsz = sqrt(DIMENSIONS * h);

  init_locks();

  printf("TEST: func=%d, N=%d, num_threads=%ld\n", func, N, get_num_threads());
  double wtime = omp_get_wtime();
  double err = feynman_pthreads_2d(a, b, N, func);
  wtime = omp_get_wtime() - wtime;
  printf("%d    %lf    %lf\n", N, err, wtime);
  printf("TEST END\n");

  destroy_locks();

  return 0;
}




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