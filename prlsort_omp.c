#include <omp.h>

#define SWAP(type, temp, x, y) type(temp) = (x); (x) = (y); (y) = (temp)

typedef int ValueType;
typedef int IndexType;

/*
Takes as arguments:
  (array) - pointer to elements to sort
  (auxil) - auxiliary space to perform merge
  (size)  - size of (array)
(array) and (auxil) must have capacity of at least (size).
Returns pointer to sorted elements. Could be either (array) or (auxil).
*/
ValueType *MergeSortOMP(ValueType *array, ValueType *auxil, const IndexType size) {
  // SORTED LENGTH LOOP. sequence of sort_len will come up AFTER the iteration
  for (IndexType sort_len = 2; sort_len < (size << 1); sort_len <<= 1) {
    // MERGE LOOP
#pragma omp parallel for
    for (IndexType start1 = 0; start1 < size; start1 += sort_len) {
      IndexType start2 = start1 + (sort_len >> 1);
      IndexType end1 = start2, end2 = start1 + sort_len, i1 = start1, i2 = start2, j = start1;
      if (end1 > size) {
        end1 = size;
      }
      if (end2 > size) {
        end2 = size;
      }
      // MERGE
      if (i2 < end2) {
        while (1) {
          if (array[i2] < array[i1]) {
            auxil[j++] = array[i2++];
            if (i2 == end2) {
              break;
            }
          } else {
            auxil[j++] = array[i1++];
            if (i1 == end1) {
              break;
            }
          }
        }
      }
      // FINISH MERGE
      while (i1 < end1) {
        auxil[j++] = array[i1++];
      }
      while (i2 < end2) {
        auxil[j++] = array[i2++];
      }
    }
    // To avoid memcpy between (array) and (auxil) SWAP pointers
    SWAP(ValueType *, temp, array, auxil);
  }
  return array;
}
