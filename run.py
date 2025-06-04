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
import time

SCRIPT_DIR = dirname(realpath(__file__))
BUILD_DIR = join(SCRIPT_DIR, 'result')
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

def get_x_N(result):
    return [int(result[0][0])]

def get_y_SPEEDUP(result, seq_result):
    return [max(float(seq_result[0][2]), 0.0000001) / max(float(result[0][2]), 0.0000001)]

def get_same(result1, result2):
    return abs(float(result1[0][1]) - float(result2[0][1])) <= ACCURACY


# IMPORTANT: Testove nazivati u formatu feynman_{tehnologija}_{DIMENSION}d
TESTS = {
    'feynman_sequential_1d': {
        'type': 'sequential',
        'args_N': [[1000], [5000], [10000], [20000]],
        'results' : {},              # key (args) : value (results)
        'x_axis' : [],
        'threads': [1]
    },
    'feynman_sequential_2d': {
        'type': 'sequential',
        'args_N': [[1000], [5000], [10000], [20000]],
        'results' : {},              # key (args) : value (results)
        'x_axis' : [],
        'threads': [1]
    },
    'feynman_sequential_3d': {
        'type': 'sequential',
        'args_N': [[1000], [5000], [10000], [20000]],
        'results' : {},              # key (args) : value (results)
        'x_axis' : [],
        'threads': [1]
    },
    'feynman_omp_1d': {
        'type': 'omp',
        'args_N': [[1000], [5000], [10000], [20000]],
        'funcs': 5,
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_omp_2d': {
        'type': 'omp',
        'args_N': [[1000], [5000], [10000], [20000]],
        'funcs': 7,
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_omp_3d': {
        'type': 'omp',
        'args_N': [[1000], [5000], [10000], [20000]],
        'funcs': 6,
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_pthreads_1d': {
        'type': 'pthreads',
        'args_N': [[1000], [5000], [10000], [20000]],
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_pthreads_2d': {
        'type': 'pthreads',
        'args_N': [[1000], [5000], [10000], [20000]],
        'threads': [1, 2, 4, 8, 16]
    },
    'feynman_pthreads_3d': {
        'type': 'pthreads',
        'args_N': [[1000], [5000], [10000], [20000]],
        'threads': [1, 2, 4, 8, 16]
    }
}

thread_colors = {
    2: 'tab:blue',
    4: 'tab:orange',
    8: 'tab:green',
    16: 'tab:red'
}


WIDTH = 2.0
GROUP_WIDTH = 1.2
BAR_WIDTH_INDEX = 0.7

# pokrece test i cuva stdout, i pravi log fajlove
# Defines a function that runs a test executable based on the test type (OpenMP for now)
def run_test(func_num: int, test_type: str, exe_name: str, args: List[int], num_threads: int) -> Result:
    # Make a copy of the current environment so we can modify it locally for this specific test (OMP_NUM_THREADS variable)
    process_env = env.copy()
    
    # Convert all provided arguments to strings, since command-line arguments must be strings
    stringified_args = [str(arg) for arg in args]

    subdir = exe_name.rsplit('_', 1)[0]  # npr. "feynman_omp_1d" -> "feynman_omp"
    exe_path = join(BUILD_DIR, subdir, exe_name)

    # If the test type is OpenMP
    if test_type == 'sequential':
        # Prepare the base arguments for executing the program (path to executable + function number)
        process_args = [exe_path]
        # Build the log file name based on the command arguments
        log_filename = ' '.join(process_args + stringified_args + [str(num_threads)])
    elif test_type == 'omp':
        # Set the number of OpenMP threads in the environment for this process
        process_env['OMP_NUM_THREADS'] = str(num_threads)
        # Prepare the base arguments for executing the program (path to executable + function number)
        process_args = [exe_path, str(func_num)]
        # Build the log file name based on the command arguments
        log_filename = ' '.join(process_args + stringified_args + [str(num_threads)])
    elif test_type == 'pthreads':       # IMPORTANT: PTHREADS reads number of threads from OMP_NUM_THREADS variable
        # Set the number of OpenMP threads in the environment for this process
        process_env['OMP_NUM_THREADS'] = str(num_threads)
        # Prepare the base arguments for executing the program (path to executable + function number)
        process_args = [exe_path]
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
    # is redirected so it can be read from Python code through process.stdout

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


def run_tests(test_name: str, test_data: Dict[str, Any], func_index: int = -1):
    # Print the name of the test being run
    print('Running', test_name, 'tests')

    # Number of functions to test (default to 1 if 'funcs' not specified)
    num_funcs = test_data['funcs'] if 'funcs' in test_data else 1
    
    # Arguments for the tests (default to empty list if 'args' not specified)
    test_args = test_data['args_N'] if 'args_N' in test_data else [[]]

    # Type of test (e.g., OpenMP)
    test_type = test_data['type']
    
    # List of thread counts to test
    threads = test_data['threads']

    # Iterate over the functions to test (if more than one function to test)
    for func_num in range(num_funcs):
        if func_index != -1 and func_index != func_num:
            continue

        # give OS time to do his stuff
        # time.sleep(5)

        # get name of sequential test
        seq_name = test_name.split("_")
        seq_name[1] = "sequential"
        seq_name = "_".join(seq_name)

        seq_res = TESTS[seq_name]["results"]
        x_axis = TESTS[seq_name]["x_axis"]            # [0, 1, 2 ,3]

        # Set up the plot figure size - sirina 15 i visina 6 inca
        plt.figure(figsize=(15, 7))
        shown_labels = set()

        # Iterate over the argument sets for the function
        for args in test_args:
            
            # give OS time to do his stuff  
            # time.sleep(5)
            
            # Iterate over the number of threads for the test
            for num_threads in threads:
                if num_threads == threads[0]:
                    continue

                print('Running test with function', func_num, 'arguments', args, 'and', num_threads, 'threads')
                
                # time.sleep(2)
                
                # Run the test with the current number of threads
                results = run_test(func_num, test_type, test_name, args, num_threads)
                
                # If results are empty, print error and exit
                if len(results) == 0:
                    print('An error occurred while getting results for ', func_num, args, num_threads, file=stderr)
                    exit(1)

                # If results do not match sequential results, print error and exit
                if not get_same(seq_res[f"{args}"], results):
                    print('Results mismatch for function ', func_num, args, num_threads, seq_res[f"{args}"], results, file=stderr)
                    print('Test FAILED')
                    exit(2)

                # Calculate speedups based on sequential and parallel results
                speedups = get_y_SPEEDUP(results, seq_res[f"{args}"])      # sta treba da pise na baru

                # Calculate bar width for plotting
                bar_width_disp = GROUP_WIDTH / (len(threads) - 1)
                bar_width = BAR_WIDTH_INDEX * GROUP_WIDTH / (len(threads) - 1)

                # Adjust x-axis positions for bars based on thread count
                x_my = x_axis[test_args.index(args)] - (GROUP_WIDTH / 2) + (threads.index(num_threads) - 1) * bar_width_disp + (bar_width / 2) + (1 - BAR_WIDTH_INDEX) * bar_width_disp /2

                # Create a bar plot for the current thread count
                label = f"threads={num_threads}"
                color = thread_colors[num_threads]

                if label not in shown_labels:
                    bar = plt.bar(x_my, speedups, label=label, width=bar_width, color=color)
                    shown_labels.add(label)
                else:
                    bar = plt.bar(x_my, speedups, width=bar_width, color=color)

    
                # Label the bars with speedup values
                plt.bar_label(bar, [round(speedup, 1) for speedup in speedups])
            
        # Set the title, labels, and ticks for the plot
        plt.title(f'Results for dimension {test_name.split("_")[2][0]} function index {func_num}')
        plt.xlabel('$N$')
        plt.ylabel('Speedup')
        plt.grid(axis='y', linestyle='--', alpha=0.2)
        plt.xticks(x_axis, test_args)
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),
           ncol=len(threads), frameon=False)


        short_name = test_name[:-3]
        svg_dir = join(BUILD_DIR, short_name)

        # Save the plot as an SVG file
        svg_filename = f'results-{test_name}-func-{func_num}.svg'
        plt.savefig(join(svg_dir, svg_filename))
        plt.close()
    
    # After all tests are run, print that the test has passed
    print('Test PASSED')


def run_sequential_tests(test_name_in: str):       # feynman_omp_1d
    # run sequential implementation
    # Print the name of the test being run
    if test_name_in is not None:
        print('Running sequential', test_name_in, 'tests')
        tmp = test_name_in.split("_")
        tmp[1] = "sequential"
        test_name_in = "_".join(tmp)            # feynman_sequential_1d
    else:
        print('Running sequential tests')

    for test_name, test_data in TESTS.items():
        if test_data["type"] != "sequential":
            continue
        if test_name_in is not None and test_name != test_name_in:
            continue

        for args in test_data["args_N"]:

            results = test_data["results"]

            results[f"{args}"] = run_test(0, test_data["type"], test_name, args, 1)
            # TODO: remove x_labels
            test_data["x_axis"] = np.array([])
            test_data["x_axis"] = np.arange(len(test_data["args_N"])) * WIDTH  # Generate x-axis values based on number of labels         #[0, 1, 2 ,3]



def main():
    if len(sys.argv) > 1:
        test_name = sys.argv[1]
        if test_name not in TESTS:
            print('Invalid test name.')
            exit(3)

        run_sequential_tests(test_name)

        if len(sys.argv) > 2:
            for i in range(2, len(sys.argv)):
                # TODO: provere 
                run_tests(test_name, TESTS[test_name], int(sys.argv[i]))
        else:
            run_tests(test_name, TESTS[test_name], -1)

    else:
        # runs all tests
        run_sequential_tests(None)
        for test_name, test_data in TESTS.items():
            if test_data["type"] == "sequential":
                continue
            run_tests(test_name, test_data)



# To run the script, use 'python run.py' to run all tests or 'python run.py test_name' to run a specific test (e.g., 'feynman_omp_1d'). Results are saved as .svg charts and log files in the 'gen' directory.
if __name__ == "__main__":
    run_make()
    main()


