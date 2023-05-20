#include <stdlib.h>

#include "../headers/typedefs.h"

IndexType* GetTestSizes(IndexType start, IndexType mid_size, IndexType max_size, double multiplier, IndexType* pn) {
  if (start > mid_size) {
    *pn = 0;
    return NULL;
  }
  const IndexType num_add = mid_size - start;
  if (multiplier <= 1.0 || max_size <= mid_size) {
    IndexType* sizes = (IndexType*)malloc((num_add + 1) * sizeof(IndexType));
    for (IndexType i = start; i <= mid_size; ++i) {
      sizes[i - start] = i;
    }
    *pn = num_add + 1;
    return sizes;
  }
  IndexType num_mult = 1;
  double prev = mid_size, next = (IndexType)(mid_size * multiplier + 0.5);
  if (next == prev) {
    ++next;
  }
  while (next < max_size) {
    ++num_mult;
    prev = next;
    next = (IndexType)(next * multiplier + 0.5);
    if (next == prev) {
      ++next;
    }
  }
  if (next - max_size < 2 * (max_size - prev)) {
    ++num_mult;
  }

  // Write
  IndexType* sizes = (IndexType*)malloc((num_add + num_mult) * sizeof(IndexType));
  for (IndexType i = 0; i < num_add; ++i) {
    sizes[i] = i + start;
  }
  prev = mid_size;
  next = (IndexType)(mid_size * multiplier + 0.5);
  if (next == prev) {
    ++next;
  }
  for (IndexType i = 0; i < num_mult; ++i) {
    sizes[num_add + i] = prev;
    prev = next;
    next = (IndexType)(next * multiplier + 0.5);
    if (next == prev) {
      ++next;
    }
  }
  sizes[num_add + num_mult - 1] = max_size;
  *pn = num_add + num_mult;
  return sizes;
}