#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "util.h"

double cpu_time(void)
{
  double value;
  // Calls the C standard function clock(), which returns the number of processor clock ticks since the start of the program’s execution.
  // That number is then divided by CLOCKS_PER_SEC (a constant that indicates how many clock ticks occur in one second).
  value = (double)clock() / (double)CLOCKS_PER_SEC;

  return value;
}

// prints current time
void timestamp(void)
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE    //  allows the size of the array to be determined within the function, while avoiding the global use of the macro afterward
}

// returns how many prime numbers are between 2 and n
int prime_number_1(int n)
{
  int total = 0;

#pragma omp parallel default(none) shared(n) \
                                   reduction(+ : total)
{
  int num_threads = omp_get_num_threads();
  int thread_id = omp_get_thread_num();

  for (int i = 2 + thread_id; i <= n; i += num_threads)
  {
    int found = 1;
    for (int j = 2; j < i; j++)
    {
      if ((i % j) == 0)
      {
        found = 0;
        break;
      }
    }
    total = total + found;
  }
} // parallel
  return total;
}

// returns how many prime numbers are between 2 and n
int prime_number_2(int n)
{
  int total = 0;

#pragma omp parallel default(none) shared(n) \
                                   reduction(+ : total)
{

// schedule can be set to 2, because check for even numbers is much faster
#pragma omp for schedule(static, 1)
  for (int i = 2; i <= n; i++)
  {
    int found = 1;
    for (int j = 2; j < i; j++)
    {
      if ((i % j) == 0)
      {
        found = 0;
        break;
      }
    }
    total = total + found;
  }
} // parallel
  return total;
}


void test(int (*func)(int), int n_lo, int n_hi, int n_factor)
{
  //double ctime;

  // printf("\n");
  // printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
  // printf("\n");
  // printf("         N        Pi          Time\n");
  // printf("\n");

  int n = n_lo;

  while (n <= n_hi)
  {
    // ctime = cpu_time();
    double wtime = omp_get_wtime();

    int primes = func(n);

    wtime = omp_get_wtime() - wtime;
    //ctime = cpu_time() - ctime;

    printf("  %8d  %8d  %14f\n", n, primes, wtime);
    n = n * n_factor;
  }

  return;
}

int (*FUNCS[])(int) = {prime_number_1, prime_number_2};

int main(int argc, char *argv[])
{
  // 'func' determines which prime number counting function to use for the test.
  // 'lo' is the starting value of N for the prime number test.
  // 'hi' is the maximum value up to which the prime number test will run.
  // 'factor' is the multiplier used to increase 'N' after each iteration of the test.
  int func = 0;
  int factor;
  int hi;
  int lo;

  // timestamp();
  // printf("\n");
  // printf("PRIME TEST\n");

  if (argc != 5)
  {
    func = 0;
    lo = 1;
    hi = 131072;
    factor = 2;
  }
  else
  {
    func = atoi(argv[1]);
    lo = atoi(argv[2]);
    hi = atoi(argv[3]);
    factor = atoi(argv[4]);
  }

  printf("TEST: func=%d, lo=%d, hi=%d, factor=%d, num_threads=%ld\n", func, lo, hi, factor, get_num_threads());
  test(FUNCS[func], lo, hi, factor);


  // printf("\n");
  // printf("PRIME_TEST\n");
  // printf("  Normal end of execution.\n");
  // timestamp();

  return 0;
}


// commands: 
// make
/*
./prime 0 1 131072 2
./prime 0 5 500000 10
./prime 0 1 65536 4
 */