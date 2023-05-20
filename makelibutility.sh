#!/bin/bash
gcc -c utility/ordering.c -o utility/ordering.o
gcc -c utility/sort_algs.c -o utility/sort_algs.o
gcc -c utility/test_gen.c -o utility/test_gen.o
gcc -c utility/test_sizes.c -o utility/test_sizes.o
ar rcs utility/libutility.a utility/ordering.o utility/sort_algs.o utility/test_gen.o utility/test_sizes.o
rm utility/ordering.o
rm utility/sort_algs.o
rm utility/test_gen.o
rm utility/test_sizes.o