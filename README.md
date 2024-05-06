
# DNA Storage System: Closest Reads/N-reads

## Overview

## User Information

### Dependencies

1. CUDA version 12.2
2. Python version 3.0+ (tested with version 3.10)
3. Compiler for C++ (tested with version 9.4)

### Instalations

In order to install the libraries used in python type the following code in the terminal:
```
pip install -r requirements.txt
```

### Organization

The repository is organized as such:
- `edit_distance`: this directory contains the 4 algorithms:
  - `gpu_opt.cu`: this file contains the optimized closest read algorithm written in parallel using `CUDA`.
  - `n_reads_gpu_opt.cu`: this file contains the optimized closest N-reads algorithm written in parallel using `CUDA`.
  - `gpu_unopt.cu`: this file contains the naive(unoptimized) closest read algorithm written in parallel using `CUDA`.
  - `cpu_unopt.cpp`: this file contains the naive(unoptimized) closest read algorithm written with `C++`
- `python_files`: this directory contains files written in python that are used for preparing the data, analysing results and creating figures.
