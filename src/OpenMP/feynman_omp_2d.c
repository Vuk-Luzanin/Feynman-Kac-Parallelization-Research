#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "util.h"

#define NUM_LOCKS   512
#define DIMENSIONS  2

static double a = 2.0;
static double b = 1.0;

int NI = 0;
int NJ = 0;

// static double b = 50.0;
static double h = 0.001;

static double stepsz;


inline double potential ( double a, double b, double x, double y )
{
  double value; 
  value = 2.0 * ( pow ( x / a / a, 2 ) + pow ( y / b / b, 2 ) ) + 1.0 / a / a + 1.0 / b / b;

  return value;
}

inline int i4_ceiling ( double x )
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
inline double r8_uniform_01(int *seed)
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

// solution with spiral matrix traversal
double feynman_6(const double a, const double b, const double h, const double stepsz, const int N) 
{
  int mat[NI * NJ][2];

  int left = 1, right = NJ + 1;
  int top = 1, bottom = NI + 1;

  int cnt = 0;
  while (left < right && top < bottom) 
  {
    // left to right
    for (int i = left; i < right; i++)
    {
      mat[cnt][0] = top;
      mat[cnt][1] = i;
      cnt++;
    }
    top++;

    // top to bottom
    for (int i = top; i < bottom; i++)
    {
      mat[cnt][0] = i;
      mat[cnt][1] = right - 1;
      cnt++;
    }
    right--;

    // if we had only one "line" left then we are done
    if (!(left < right && top < bottom))
      break;

    // right to left
    for (int i = right - 1; i >= left; i--)
    {
      mat[cnt][0] = bottom - 1;
      mat[cnt][1] = i;
      cnt++;
    }
    bottom--;

    //bottom to top
    for (int i = bottom - 1; i >= top; i--)
    {
      mat[cnt][0] = i;
      mat[cnt][1] = left;
      cnt++;
    }
    left++;

  }

  int seed = 123456789;   
  double err = 0.0;
  int n_inside = 0; 

#pragma omp parallel for schedule(guided) \
                         reduction(+:err, n_inside)
  for (int i = 0; i < NI*NJ; i++)
  {
    double x = ((double)(NI - mat[i][0]) * (-a) + (double)(mat[i][0] - 1) * a) / (double)(NI - 1);
    double y = ((double)(NJ - mat[i][1]) * (-b) + (double)(mat[i][1] - 1) * b) / (double)(NJ - 1);

    double chk = pow(x / a, 2) + pow(y / b, 2);

    if ( 1.0 < chk )
    {
      // tacka nije unutar 2-D elipsoida
      continue;
    }

    n_inside++;
    double w_exact = exp ( pow ( x / a, 2 )
                + pow ( y / b, 2 ) - 1.0 );
    double wt = 0;

    for (int trial = 0; trial < N; trial++)
    {   
        // seed is private variable, so numbers can generate uniformly
        int localseed = seed + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x;
        double x2 = y;
  
        double w = 1.0;
        chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0)
        {
#ifdef SMALL_STEP
          double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
          double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
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
        wt += w;
    }

    wt = wt / ( double ) ( N ); 
    err = err + pow ( w_exact - wt, 2 );
  }

  return err = sqrt ( err / ( double ) ( n_inside ) );
}

double feynman_0(const double a, const double b, const double h, const double stepsz, const int N) 
{
  int seed = 123456789;
  double err = 0.0;
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

  double* w_exact = calloc((NI+1) * (NJ+1), sizeof(double));
  double* wt = calloc((NI+1) * (NJ+1), sizeof(double));

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, n_inside, w_exact, wt, err, seed, NI, NJ)
{
  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++ )
    {
      // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
      // private by default in whole parallel region
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);


      if ( 1.0 < chk )
      {
        // tacka nije unutar 2-D elipsoida
        continue;
      }

        // tacka je unutar 2-D elipsoida
#pragma omp single nowait
{
        n_inside++;
        w_exact[i * (NJ+1) + j] = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);
} // single

      // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
      // x, y and chk are private for thread by default
      // ovde sme nowait, jer niti ne zavise jedna od druge (sinhronizacija postoji vec zbog reduction, pa je i ona zadovoljena)
#pragma omp for nowait reduction(+:wt[:(NI+1)*(NJ+1)])       // probati sa pomeranjem pokazivaca, p dereferencirati
      for (int trial = 0; trial < N; trial++)
      {
        // seed is private variable, so numbers can generate uniformly
        int localseed = seed + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x;
        double x2 = y;
  
        double w = 1.0;
        chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0)
        {
#ifdef SMALL_STEP
          double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
          double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
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
        wt[i * (NJ+1) + j] += w;
      }
    }
  }

  // barijera samo na kraju trial petlje za svaku nit (u for je dodato nowait, kako se niti medjusobno ne bi cekale)
  // po defaultu, na kraju for petlje se sve niti cekaju
#pragma omp barrier

#pragma omp for collapse(2) reduction(+:err)
  for (int i=0; i<=NI; i++)
  {
    for (int j=0; j<=NJ; j++)
    {
      if (w_exact[i * (NJ+1) + j] == 0.0)
      {
        // kada tacka nije unutar elipsoida
        continue;
      }
      // kvadrat razlike tacne i numericki dobijene vrednosti
      err += pow(w_exact[i * (NJ+1) + j] - wt[i * (NJ+1) + j] / (double)(N), 2);
    }
  }
} // parallel
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}

double feynman_1(const double a, const double b, const double h, const double stepsz, const int N)
{
  double err = 0.0;
  int n_inside = 0;

  double w_exact[NI+1][NJ+1];

  for (int i = 1; i <= NI; i++) {
    for (int j = 1; j <= NJ; j++) {
      // Interpolacija koordinata
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

      w_exact[i][j] = 0;

      if (chk > 1.0) continue; // Van elipsoida

      // Tačka je unutar 2D elipsoida
      w_exact[i][j] = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);
      n_inside++;

      double local_sum = 0.0;

      // Paralelizovana petlja po pokušajima
#pragma omp parallel for reduction(+:local_sum)
      for (int trial = 0; trial < N; trial++) 
      {
        // seed is private variable, so numbers can generate uniformly
        int localseed = 123456789 + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x, x2 = y;
        double w = 1.0;
        double chk_inner = 0.0;

        while (chk_inner < 1.0) 
        {
#ifdef SMALL_STEP
          double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt(DIMENSIONS * h);
          double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt(DIMENSIONS * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
#endif
          double vs = potential(a, b, x1, x2);
          x1 += dx;
          x2 += dy;
          double vh = potential(a, b, x1, x2);

          double we = (1.0 - h * vs) * w;
          w = w - 0.5 * h * (vh * we + vs * w);

          chk_inner = pow(x1 / a, 2) + pow(x2 / b, 2);
        }

        local_sum += w;
      }

      err += pow(w_exact[i][j] - local_sum / (double)N, 2);
    }
  }

  return sqrt(err / (double)n_inside);
}


// solution with task per walk
double feynman_2(const double a, const double b, const double h, const double stepsz, const int N) 
{
  int seed = 123456789;
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

  // for every point in grid - must be initialized (moguce je da se desi da neki taskovi krenu da pristupaju matrici pre postavljanja na 0 -> citaju neinicijalizovanu memoriju)
  double* w_exact = calloc((NI+1) * (NJ+1), sizeof(double));
  double* wt = calloc((NI+1) * (NJ+1), sizeof(double));

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, n_inside, w_exact, wt, seed, NI, NJ) 
{
#pragma omp single
{
  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++ )
    {
      // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

      w_exact[i * (NJ+1) + j] = 0.0;
      wt[i * (NJ+1) + j] = 0.0;

      if ( 1.0 < chk )
      {
        // tacka nije unutar 1-D elipsoida
        continue;
      }

      // tacka je unutar 2-D elipsoida
      n_inside++;
 
      // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
      w_exact[i * (NJ+1) + j] = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);

      // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
      for (int trial = 0; trial < N; trial++)
      {
#pragma omp task shared(wt)
{
        // seed is private variable, so numbers can generate uniformly
        int localseed = seed + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x;
        double x2 = y;
  
        double w = 1.0;
        chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0)
        {
#ifdef SMALL_STEP
            double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
            double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
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
// jer se menja u kontekstu vise taskova
#pragma omp atomic
        wt[i * (NJ+1) + j] += w;
} // task
      }
    }
  }
} // single
} // parallel
  // na kraju izracunati gresku po svim osama - u matrici wt
  double err = 0.0;
  for (int i = 0; i <= NI; ++i)
  {
    for (int j = 0; j <= NJ; ++j)
    { 
      if (w_exact[i * (NJ+1) + j] == 0.0)
      {
        // kada tacka nije unutar elipsoida
        continue;
      }
      err += pow(w_exact[i * (NJ+1) + j] - (wt[i * (NJ+1) + j] / (double)(N)), 2);
    }
  }
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}



omp_lock_t locks[NUM_LOCKS];

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


// solution with task per walk
double feynman_3(const double a, const double b, const double h, const double stepsz, const int N) 
{
  int seed = 123456789;     // set to be static -> to be global (shared) by default
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

  // for every point in grid - must be initialized (moguce je da se desi da neki taskovi krenu da pristupaju matrici pre postavljanja na 0 -> citaju neinicijalizovanu memoriju)
  double* w_exact = calloc((NI+1) * (NJ+1), sizeof(double));
  double* wt = calloc((NI+1) * (NJ+1), sizeof(double));

  for (int i = 0; i < NUM_LOCKS; i++) 
  {
    omp_init_lock(&locks[i]);
  }

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, n_inside, w_exact, wt, locks, seed, NI, NJ) 
{
#pragma omp single
{
  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++ )
    {
      // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

      w_exact[i * (NJ+1) + j] = 0.0;
      wt[i * (NJ+1) + j] = 0.0;

      if ( 1.0 < chk )
      {
        // tacka nije unutar 1-D elipsoida
        continue;
      }

      // tacka je unutar 2-D elipsoida
      n_inside++;
 
      // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
      w_exact[i * (NJ+1) + j] = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);

      // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
      for (int trial = 0; trial < N; trial++)
      {
#pragma omp task shared(wt)
{
        // seed is private variable, so numbers can generate uniformly
        int localseed = seed + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x;
        double x2 = y;
  
        double w = 1.0;
        chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0)
        {
#ifdef SMALL_STEP
            double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
            double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
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

        // koriscenje lock-a
        // sigurno jer svaki task simulira jedno kretanje iz tacke, pa ce se svi taskovi za 1 tacku (isti i, j) sinhronizovati nad istom bravom u nizu locks
        // 
        int lock_id = get_lock_index(i, j);    // izracunaj index lock-a koji je potreban
        omp_set_lock(&locks[lock_id]);
        wt[i * (NJ+1) + j] += w;
        omp_unset_lock(&locks[lock_id]);
} // task
      }
    }
  }
} // single
} // parallel
  // na kraju obracunati gresku po svim osama - u matrici wt
  double err = 0.0;
  for (int i = 0; i <= NI; ++i)
  {
    for (int j = 0; j <= NJ; ++j)
    { 
      if (w_exact[i * (NJ+1) + j] == 0.0)
      {
        // kada tacka nije unutar elipsoida
        continue;
      }
      err += pow(w_exact[i * (NJ+1) + j] - (wt[i * (NJ+1) + j] / (double)(N)), 2);
    }
  }
  // oslobađanje brava
  for (int i = 0; i < NUM_LOCKS; i++) {
    omp_destroy_lock(&locks[i]);
  }

  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}


// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// solution with for collapse of outer loops and reduction of error
double feynman_4(const double a, const double b, const double h, const double stepsz, const int N) 
{
  int seed = 123456789;
  double err = 0.0;
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, seed, NI, NJ) \
                                   reduction(+ : err) \
                                   reduction(+ : n_inside)
{
#pragma omp for collapse(2)
  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++ )
    {
      // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

      double w_exact = 0.0;
      double wt = 0.0;

      if ( 1.0 < chk )
      {
        // tacka nije unutar 1-D elipsoida
        continue;
      }

      // tacka je unutar 2-D elipsoida
      n_inside++;
 
      // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
      w_exact = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);
      wt = 0.0;

      // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
      for (int trial = 0; trial < N; trial++)
      {
        // seed is private variable, so numbers can generate uniformly
        int localseed = seed + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x;
        double x2 = y;
  
        double w = 1.0;
        chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0)
        {
#ifdef SMALL_STEP
            double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
            double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
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
        wt = wt + w;
      }
      // srednja vrenost tezine za N pokusaja
      wt = wt / (double)(N);

      // kvadrat razlike tacne i numericki dobijene vrednosti
      err += pow(w_exact - wt, 2);
    }
  }
} // parallel
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}

// solution with for directive for outer loop and reduction of error
double feynman_5(const double a, const double b, const double h, const double stepsz, const int N) 
{
  int seed = 123456789;
  double err = 0.0;
  int n_inside = 0;   // broj tacaka unutar elipsoida (unutar mreze)

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, seed, NI, NJ) \
                                   reduction(+ : err) \
                                   reduction(+ : n_inside)
{
#pragma omp for schedule(dynamic)
  for (int i = 1; i <= NI; i++)
  {
    for (int j = 1; j <= NJ; j++ )
    {
      // interpolacija koordinata kako bi se dobilo kada je i = 1 -> x = -a, kada je i = ni -> x = a
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

      double w_exact = 0.0;
      double wt = 0.0;

      if ( 1.0 < chk )
      {
        // tacka nije unutar 1-D elipsoida
        continue;
      }

      // tacka je unutar 2-D elipsoida
      n_inside++;
 
      // analitička vrednost funkcije gustine/potencijala u tački unutar elipsoida - referentna vrednost koju poredimo u odnosu na numericku - wt
      w_exact = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);
      wt = 0.0;

      // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
      for (int trial = 0; trial < N; trial++)
      {
        // seed is private variable, so numbers can generate uniformly
        int localseed = seed + omp_get_thread_num() * 997 + trial;      // LEAP-FROG

        double x1 = x;
        double x2 = y;
  
        double w = 1.0;
        chk = 0.0;

        // kretanje cestice - dok se nalazi unutar elipsoida
        while (chk < 1.0)
        {
#ifdef SMALL_STEP
            double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
            double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt((DIMENSIONS*1.0) * h);
#else
          double ut = r8_uniform_01(&localseed);
          double dx = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;

          ut = r8_uniform_01(&localseed);
          double dy = (ut < 0.5) ? ((r8_uniform_01(&localseed) - 0.5) < 0.0 ? -stepsz : stepsz) : 0.0;
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
        wt = wt + w;
      }
      // srednja vrenost tezine za N pokusaja
      wt = wt / (double)(N);

      // kvadrat razlike tacne i numericki dobijene vrednosti
      err += pow(w_exact - wt, 2);
    }
  }
} // parallel
  // root-mean-square (RMS) error
  return sqrt(err / (double)(n_inside));
}


double (*FUNCS[])(const double, const double, const double, const double, const int) = {feynman_0, feynman_1, feynman_2, feynman_3, feynman_4, feynman_5, feynman_6};

int main ( int argc, char **argv )
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

  stepsz = sqrt(DIMENSIONS * h);

  if (a < b)
  {
    NJ = 6;
    NI = 1 + i4_ceiling (b / a) * (NJ - 1);
  }
  else
  {
    NI = 6;
    NJ = 1 + i4_ceiling (a / b) * (NI - 1);
  }

  printf("TEST: func=%d, N=%d, num_threads=%ld\n", func, N, get_num_threads());
  double wtime = omp_get_wtime();
  double err = FUNCS[func](a, b, h, stepsz, N);
  wtime = omp_get_wtime() - wtime;
  printf("%d    %lf    %lf\n", N, err, wtime);
  printf("TEST END\n");

  return 0;
}
