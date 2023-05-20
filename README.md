# Parallel Sort using OMP and MPI

This repository contains implementations of parallel sorting algorithms, benchmark tests for them and results of tests.

## Requirements

- OS : Linux
- C compiler, preferably gcc
- Any OMP implementation
- Any MPI implementation

## Details

### Theory

All implementations are based on [**MergeSort**](https://en.wikipedia.org/wiki/Merge_sort).

Repository provides three possible implementations: **OMP**, **MPI** and **MIX** (**OMP** + **MPI**).

**OMP** being the easiest one uses the fact, that every subarray can be sorted independently.

**MPI** sorts array on every process independently using [**MergeSort**](https://en.wikipedia.org/wiki/Merge_sort). After that arrays are merged through [sorting network](https://en.wikipedia.org/wiki/Sorting_network). Here two approaches were considered:
1. [Batcher odd–even mergesort](https://en.wikipedia.org/wiki/Batcher_odd–even_mergesort)
2. [Bitonic sorter](https://en.wikipedia.org/wiki/Bitonic_sorter) (alternative one)

While the [bitonic sorter](https://en.wikipedia.org/wiki/Bitonic_sorter) is easier to build on powers of 2 and harder on other number of processes, the [odd–even mergesort](https://en.wikipedia.org/wiki/Batcher_odd–even_mergesort) is built universally on any number of processes.

### Tests

Tests were run on a university cluster using one node with 40 cores.

## Contents

- ```prlsort_omp.c```, ```prlsort_mpi.c```, ```prlsort_mix.c``` - implementations of algorithms
- ```test_omp.c```, ```test_mpi.c```, ```test_mix.c``` - benchmark tests, each is independent with own ```main``` function. Build and run instructions are on top of each file
- ```qs_prlsort``` - **sbatch** file launching tests on a cluster
- ```test_results``` directory contains results from running tests on a cluster
- ```utility``` and ```headers``` directories contain ```utility``` library and ```.h``` files for them
  - ```ordering.c``` implements [Batcher odd–even mergesort](https://en.wikipedia.org/wiki/Batcher_odd–even_mergesort) and [bitonic sorter](https://en.wikipedia.org/wiki/Bitonic_sorter)
  - ```sort_algs.c``` implements [**HeapSort**](https://en.wikipedia.org/wiki/Heapsort), [**MergeSort**](https://en.wikipedia.org/wiki/Merge_sort) and [**QuickSort**](https://medium.com/@nehasangeetajha/3-way-quick-sort-18d2dcc5b06b) (three-way partition)
  - ```test_gen.c``` provides functions to create tests: sorted, reversed, randomized and slightly randomized
  - ```test_sizes.c``` contains function ```GetTestSizes``` which returns sorted array of sizes for tests to generate
- ```makelibutility.sh``` - ```bash``` script which builds ```utility``` libraries if there were changes