#!/usr/bin/env python
from subprocess import Popen, PIPE
import subprocess
import sys
from typing import Any, Dict, List
from os import environ as env
from os.path import dirname, realpath, join
from sys import argv, exit, stderr
from matplotlib import pyplot as plt
import numpy as np

SCRIPT_DIR = dirname(realpath(__file__))
BUILD_DIR = join(SCRIPT_DIR, 'gen')
ACCURACY = 0.01

Result = List[List[str]]

def run_make():
    try:
        # runs make command
        subprocess.run(['make'], check=True)
        print("Make completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Make failed with error: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Make command not found. Please ensure that 'make' is installed.", file=sys.stderr)
        sys.exit(1)

TESTS = {
    'feynman_omp_1d': {
        'type': 'omp',
        'args': [[1000], [5000], [10000], [20000]],
        'funcs': 1,
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY),
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_omp_2d': {
        'type': 'omp',
        'args': [[1000], [5000], [10000], [20000]],
        'funcs': 1,
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY),
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_omp_3d': {
        'type': 'omp',
        'args': [[1000], [5000], [10000], [20000]],
        'funcs': 3,
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY),
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_pthreads_1d': {
        'type': 'pthreads',
        'args': [[1000], [5000], [10000], [20000]],
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY),
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_pthreads_2d': {
        'type': 'pthreads',
        'args': [[1000], [5000], [10000], [20000]],
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY),
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_pthreads_3d': {
        'type': 'pthreads',
        'args': [[1000], [5000], [10000], [20000]],
        'x': lambda result: [int(result[0][0])],
        'y': lambda result, seq_result: [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)],
        'same': lambda result1, result2: (abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY),
        'threads': [1, 2, 4, 8, 16]
    }
}

WIDTH = 1.0

# Defines a function that runs a test executable based on the test type (OpenMP for now)
def run_test(func_num: int, test_type: str, exe_name: str, args: List[int], num_threads: int) -> Result:
    # Make a copy of the current environment so we can modify it locally for this specific test (OMP_NUM_THREADS variable)
    process_env = env.copy()
    
    # Convert all provided arguments to strings, since command-line arguments must be strings
    stringified_args = [str(arg) for arg in args]

    # If the test type is OpenMP
    if test_type == 'omp':
        # Set the number of OpenMP threads in the environment for this process
        process_env['OMP_NUM_THREADS'] = str(num_threads)
        # Prepare the base arguments for executing the program (path to executable + function number)
        process_args = [f'{BUILD_DIR}/{exe_name}', str(func_num)]
        # Build the log file name based on the command arguments
        log_filename = ' '.join(process_args + stringified_args + [str(num_threads)])
    elif test_type == 'pthreads':       # IMPORTANT: PTHREADS reads number of threads from OMP_NUM_THREADS variable
        # Set the number of OpenMP threads in the environment for this process
        process_env['OMP_NUM_THREADS'] = str(num_threads)
        # Prepare the base arguments for executing the program (path to executable + function number)
        process_args = [f'{BUILD_DIR}/{exe_name}']
        # Build the log file name based on the command arguments
        log_filename = ' '.join(process_args + stringified_args + [str(num_threads)])

    # If the test type is unknown, raise an exception
    else:
        raise BaseException('Unknown test type.')

    # Append the stringified arguments to the command to be executed
    process_args += stringified_args
    
    # Launch the process with the specified environment and capture its stdout
    process = Popen(process_args, env=process_env, stdout=PIPE) # env=process_env -> Environment variables that will apply only to this specific process.

    # stdout=PIPE means that the standard output of the process (what the program would normally print to the screen)
    #  is redirected so it can be read from Python code through process.stdout


    # Wait for the process to finish, and if it returns a non-zero status or has no output, return an empty list
    if process.wait() != 0 or not process.stdout:
        return []

    # List to store the parsed results
    results = []

    # Define the log file name and its full path
    log_filename = join(BUILD_DIR, f'{log_filename}.log')

    # Open the log file for writing
    with open(log_filename, 'w', encoding='utf-8') as log_file:
        # Read each line of the process output
        for line in process.stdout:
            # Decode the line from bytes to string
            print(f'Line: "{line}"')
            line = line.decode('utf-8')
            # Write the line to the log file
            log_file.write(line)
            # If the line doesn't start with 'TEST', split it into components and store it in results
            if not line.startswith('TEST'):
                results.append(line.split())

    # Return the collected results
    return results


def run_tests(test_name: str, test_data: Dict[str, Any]):
    # Print the name of the test being run
    print('Running', test_name, 'tests')
    
    # Get the test data corresponding to the test_name from the TESTS dictionary
    test_data = TESTS[test_name]

    # Number of functions to test (default to 1 if 'funcs' not specified)
    num_funcs = test_data['funcs'] if 'funcs' in test_data else 1
    
    # Arguments for the tests (default to empty list if 'args' not specified)
    test_args = test_data['args'] if 'args' in test_data else [[]]

    # Function to get x-axis labels and y-axis data from the test data
    get_x_axis = test_data['x']
    get_y_axis = test_data['y']

    # Function to check if results from different thread counts are the same
    check_same = test_data['same']  # lambda function

    # Type of test (e.g., OpenMP)
    test_type = test_data['type']
    
    # List of thread counts to test
    threads = test_data['threads']

    # Iterate over the functions to test (if more than one function to test)
    for func_num in range(num_funcs):
        
        # Iterate over the argument sets for the function
        for arg_idx, args in enumerate(test_args):
            
            # Initialize variables to store results and x-axis data
            seq_results = []
            x_axis = np.array([])
            x_labels = []
            
            # Set up the plot figure size
            plt.figure(figsize=(15, 6))
            
            # Iterate over the number of threads for the test
            for num_threads in threads:
                print('Running test with function', func_num, 'arguments', args, 'and', num_threads, 'threads')

                # If running with the first number of threads, get sequential results
                if num_threads == threads[0]:
                    seq_results = run_test(func_num, test_type, test_name, args, num_threads)
                    x_labels = get_x_axis(seq_results)  # Get x-axis labels from results
                    x_axis = np.arange(len(x_labels)) * WIDTH  # Generate x-axis values based on number of labels
                
                else:
                    # Run the test with the current number of threads
                    results = run_test(func_num, test_type, test_name, args, num_threads)
                    
                    # If results are empty, print error and exit
                    if len(results) == 0:
                        print('An error occurred while getting results for ', func_num, args, num_threads, file=stderr)
                        exit(1)

                    # If results do not match sequential results, print error and exit
                    if not check_same(seq_results, results):
                        print('Results mismatch for ', func_num, args, num_threads, seq_results, results, file=stderr)
                        print('Test FAILED')
                        exit(2)
                    
                    # Calculate speedups based on sequential and parallel results
                    speedups = get_y_axis(results, seq_results)
                    
                    # Calculate bar width for plotting
                    bar_width = WIDTH / (len(threads) - 1)
                    
                    # Adjust x-axis positions for bars based on thread count
                    x_my = x_axis - (WIDTH / 2) + (threads.index(num_threads) - 1) * bar_width + (bar_width / 2)
                    
                    # Create a bar plot for the current thread count
                    bar = plt.bar(x_my, speedups, label=f'threads={num_threads}', width=bar_width)
                    
                    # Label the bars with speedup values
                    plt.bar_label(bar, [round(speedup, 1) for speedup in speedups])
            
            # Set the title, labels, and ticks for the plot
            plt.title(f'Results for function index {func_num} and arguments {args}')
            plt.xlabel('$N$')
            plt.ylabel('Speedup')
            plt.xticks(x_axis, x_labels)
            plt.legend()

            # Save the plot as an SVG file
            plt.savefig(join(BUILD_DIR, f'results-{test_name}-{func_num}-{arg_idx}.svg'))
            plt.close()
    
    # After all tests are run, print that the test has passed
    print('Test PASSED')


def main():
    if len(sys.argv) > 1:
        test_name = sys.argv[1]
        if test_name not in TESTS:
            print('Invalid test name.')
            exit(3)
        run_tests(test_name, TESTS[test_name])
    else:
        # runs all tests
        for test_name, test_data in TESTS.items():
            run_tests(test_name, test_data)

# To run the script, use 'python run.py' to run all tests or 'python run.py test_name' to run a specific test (e.g., 'feynman_omp_1d'). Results are saved as .svg charts and log files in the 'gen' directory.
if __name__ == "__main__":
    run_make()
    main()


"""

(name, function number (0-worksharing, 1-tasks), number_of_points)
./gen/feyman 0 1000
./gen/feyman 0 5000
./gen/feyman 0 10000
./gen/feyman 0 20000
 
time ./md 

"""