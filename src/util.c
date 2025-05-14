#include "util.h"
#include <unistd.h>
#include <stdlib.h>

// Function to get the number of threads to be used in parallel execution
long get_num_threads(void)
{
    // Get the number of available logical CPU cores on the system
    long num_threads = sysconf(_SC_NPROCESSORS_ONLN);

    // Check if the environment variable 'OMP_NUM_THREADS' is set (exported)
    char* num_threads_str = getenv("OMP_NUM_THREADS");

    // If the environment variable is set, convert its value from string to long
    if (num_threads_str != NULL)
    {
        num_threads = atol(num_threads_str);
    }

    return num_threads;
}
