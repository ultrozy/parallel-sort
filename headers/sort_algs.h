#ifndef SORT_ALGS_H_
#define SORT_ALGS_H_

#include "typedefs.h"

void HeapSort(ValueType *const array, const IndexType n);
ValueType *MergeSort(ValueType * array, ValueType *, const IndexType n);
void QuickSort(ValueType *const array, const IndexType left, const IndexType right);

double MySqrt(double x);

#endif