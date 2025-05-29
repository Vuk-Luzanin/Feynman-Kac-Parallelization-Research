#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "util.h"

#define NUM_LOCKS   512
#define DIMENSIONS  2
#define NI          6
#define NJ          11
// #define NI          51
// #define NJ          6

static double a = 2.0;
static double b = 1.0;
// static double b = 20.0;
static double h = 0.001;

static double stepsz;


inline double potential ( double a, double b, double x, double y )
{
  double value; 
  value = 2.0 * ( pow ( x / a / a, 2 ) + pow ( y / b / b, 2 ) ) + 1.0 / a / a + 1.0 / b / b;

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

  double w_exact[NI+1][NJ+1] = {{0}};
  double wt[NI+1][NJ+1] = {{0}};

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, n_inside, w_exact, wt, err, seed)
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
        w_exact[i][j] = exp(pow(x / a, 2) + pow(y / b, 2) - 1.0);
} // single

      // pustamo N tacaka iz izabrane koordinate - visestruki pokusaji kako bi se dobila bolja aproksimacija
      // x, y and chk are private for thread by default
      // ovde sme nowait, jer niti ne zavise jedna od druge (sinhronizacija postoji vec zbog reduction, pa je i ona zadovoljena)
#pragma omp for nowait reduction(+:wt[i][j])
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
        wt[i][j] += w;
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
      if (w_exact[i][j] == 0.0)
      {
        // kada tacka nije unutar elipsoida
        continue;
      }
      // kvadrat razlike tacne i numericki dobijene vrednosti
      err += pow(w_exact[i][j] - wt[i][j] / (double)(N), 2);
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

  double w_exact[NI+1][NJ+1] = {{0}};

  for (int i = 1; i <= NI; i++) {
    for (int j = 1; j <= NJ; j++) {
      // Interpolacija koordinata
      double x = ((double)(NI - i) * (-a) + (double)(i - 1) * a) / (double)(NI - 1);
      double y = ((double)(NJ - j) * (-b) + (double)(j - 1) * b) / (double)(NJ - 1);
      double chk = pow(x / a, 2) + pow(y / b, 2);

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
  double w_exact[NI+1][NJ+1] = {{0}};
  double wt[NI+1][NJ+1] = {{0}};

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, n_inside, w_exact, wt, seed) 
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
        wt[i][j] += w;
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
      if (w_exact[i][j] == 0.0)
      {
        // kada tacka nije unutar elipsoida
        continue;
      }
      err += pow(w_exact[i][j] - (wt[i][j] / (double)(N)), 2);
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
  double w_exact[NI+1][NJ+1] = {{0}};
  double wt[NI+1][NJ+1] = {{0}};

  for (int i = 0; i < NUM_LOCKS; i++) 
  {
    omp_init_lock(&locks[i]);
  }

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, n_inside, w_exact, wt, locks, seed) 
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
        wt[i][j] += w;
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
      if (w_exact[i][j] == 0.0)
      {
        // kada tacka nije unutar elipsoida
        continue;
      }
      err += pow(w_exact[i][j] - (wt[i][j] / (double)(N)), 2);
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

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, seed) \
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

#pragma omp parallel default(none) shared(a, b, h, stepsz, N, seed) \
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


double (*FUNCS[])(const double, const double, const double, const double, const int) = {feynman_0, feynman_1, feynman_2, feynman_3, feynman_4, feynman_5, feynman_6, feynman_7};

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

  printf("TEST: func=%d, N=%d, num_threads=%ld\n", func, N, get_num_threads());
  double wtime = omp_get_wtime();
  double err = FUNCS[func](a, b, h, stepsz, N);
  wtime = omp_get_wtime() - wtime;
  printf("%d    %lf    %lf\n", N, err, wtime);
  printf("TEST END\n");

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