# Feynman Kac Parallel Research

Repository focused on research and implementation of parallel algorithms for the Feynman-Kac problem. Covers CPU and GPU parallelization using programming models like OpenMP, Pthreads, CUDA, and others. The goal is to explore performance optimizations and scalability across different hardware.

## Structure

- `run.py` is a test runner that automatically builds the project, runs tests, plots speedup graphs, and logs the output.  
- `Makefile` builds all tasks.  
- `src` contains source files named as `name.c`.  
- Tests are executed with:


```bash
python3 run.py <test_name>
```

Where <test_name> is optional (feynman_omp_1d, feynman_omp_2d, feynman_omp_3d, feynman_pthreads_1d, feynman_pthreads_2d, feynman_pthreads_3d). If omitted, all tests will be executed.


In the `result/` folder, you can find:

- Compiled binaries
- Log files
- Acceleration plots relative to the sequential implementation

These files are **organized into subfolders** based on the **technology used**.  
The folder is **automatically generated**.


## Report

Attached are the PDF documents [ReportOpenMP.pdf](./Reports/ReportOpenMP.pdf) and [ReportPthreads.pdf](./Reports/ReportPthreads.pdf) (in Reports/ directory), in which the method of parallelization and discussion are detailed in Serbian, along with examples of result graphs.
