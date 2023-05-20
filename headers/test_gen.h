#ifndef TEST_GEN_H_
#define TEST_GEN_H_

#include "typedefs.h"

int *GetSortedInts(const IndexType n);
int *GetRandomInts(const IndexType n, const IndexType min, const IndexType max, unsigned int seed);
void SetRandomInts(int* array, const IndexType n, const IndexType min, const IndexType max, unsigned int seed);
ValueType *GetReverseTest(ValueType *const array, const IndexType n);
ValueType *GetRandomizeTest(ValueType *const array, const IndexType n, unsigned int seed);
ValueType *GetSlightlyRandomizedTest(ValueType *array, const IndexType n, IndexType window, unsigned int seed);

#endif