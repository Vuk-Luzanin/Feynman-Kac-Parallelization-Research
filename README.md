# feynman-kac-parallel-research

Repository focused on research and implementation of parallel algorithms for the Feynman-Kac problem. Covers CPU and GPU parallelization using programming models like OpenMP, Pthreads, CUDA, and others. The goal is to explore performance optimizations and scalability across different hardware.

## Structure

- `run.py` is a test runner that automatically builds the project, runs tests, plots speedup graphs, and logs the output.  
- `Makefile` builds all tasks.  
- `src` contains source files named as `name.c`.  
- Tests are executed with:

```bash
python3 run.py <test_name>
```

Where <test_name> is optional (prime, feynman, moldyn). If omitted, all tests will be executed.

Compiled binaries, logs, and generated SVG charts are placed in the gen/ directory.

## Report

Attached is the PDF document [Report.pdf](./Reports/Report.pdf) (in Reports/ directory), in which the method of parallelization and discussion are detailed in Serbian, along with examples of result graphs.
