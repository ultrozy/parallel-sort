#include "../headers/ordering.h"

#include <stdint.h>
#include <stdlib.h>

struct MyVector GetVector(int32_t capacity) {
  if (capacity <= 0) {
    struct MyVector vec = {0, 0, NULL};
    return vec;
  }
  struct MyVector vec = {0, capacity, (short*)malloc(capacity * sizeof(short))};
  return vec;
}
void PushBack(struct MyVector* pvec, short value) {
  if (pvec->data == NULL) {
    ++pvec->capacity;
    pvec->data = (short*)malloc(sizeof(short));
  } else if (pvec->capacity == pvec->size) {
    pvec->capacity = (pvec->capacity * 3 + 1) >> 1;
    pvec->data = (short*)realloc(pvec->data, pvec->capacity * sizeof(short));
  }
  pvec->data[pvec->size] = value;
  ++pvec->size;
}

void OEmergenet(SortNet* all_mergenets, short n) {
  all_mergenets[n] = (SortNet)malloc(n * sizeof(struct MyVector));
  SortNet net = all_mergenets[n];
  for (short i = 0; i < n; ++i) {
    net[i] = GetVector(0);
  }
  short a_size = ((n + 1) >> 1), b_size = n - a_size;
  short a_odd = (a_size >> 1), b_odd = (b_size >> 1), odds = a_odd + b_odd;
  short a_even = a_size - a_odd, b_even = b_size - b_odd, evens = a_even + b_even;

  // MERGE EVEN
  //// A
  for (short i = 0; i < a_size; i += 2) {
    short converted = i >> 1;
    for (int32_t j = 0; j < all_mergenets[evens][converted].size; ++j) {
      short elem = all_mergenets[evens][converted].data[j];
      PushBack(&net[i], elem < a_even ? elem << 1 : a_size + ((elem - a_even) << 1));
    }
  }
  //// B
  for (short i = a_size; i < n; i += 2) {
    short converted = ((i - a_size) >> 1) + a_even;
    for (int32_t j = 0; j < all_mergenets[evens][converted].size; ++j) {
      short elem = all_mergenets[evens][converted].data[j];
      PushBack(&net[i], elem < a_even ? elem << 1 : a_size + ((elem - a_even) << 1));
    }
  }
  // MERGE ODD
  //// A
  for (short i = 1; i < a_size; i += 2) {
    short converted = i >> 1;
    for (int32_t j = 0; j < all_mergenets[odds][converted].size; ++j) {
      short elem = all_mergenets[odds][converted].data[j];
      PushBack(&net[i], elem < a_odd ? (elem << 1) + 1 : a_size + 1 + ((elem - a_odd) << 1));
    }
  }
  //// B
  for (short i = a_size + 1; i < n; i += 2) {
    short converted = ((i - a_size) >> 1) + a_odd;
    for (int32_t j = 0; j < all_mergenets[odds][converted].size; ++j) {
      short elem = all_mergenets[odds][converted].data[j];
      PushBack(&net[i], elem < a_odd ? (elem << 1) + 1 : a_size + 1 + ((elem - a_odd) << 1));
    }
  }

  // EMERGE
  for (short i = 1; i < n - 1; i += 2) {
    PushBack(&net[i], i + 1);
    PushBack(&net[i + 1], i);
  }
}
void OEsortnet(SortNet* all_sortnets, SortNet* all_mergenets, short n) {
  all_sortnets[n] = (SortNet)malloc(n * sizeof(struct MyVector));
  short a_size = ((n + 1) >> 1), b_size = n - a_size;
  SortNet net = all_sortnets[n], merge = all_mergenets[n];
  SortNet asort = all_sortnets[a_size], bsort = all_sortnets[b_size];

  // A
  for (short i = 0; i < a_size; ++i) {
    net[i] = GetVector(asort[i].size + merge[i].size);
    for (int32_t j = 0; j < asort[i].size; ++j) {
      PushBack(&net[i], asort[i].data[j]);
    }
    for (int32_t j = 0; j < merge[i].size; ++j) {
      PushBack(&net[i], merge[i].data[j]);
    }
  }

  // B
  for (short i = a_size; i < n; ++i) {
    net[i] = GetVector(bsort[i - a_size].size + merge[i].size);
    for (int32_t j = 0; j < bsort[i - a_size].size; ++j) {
      PushBack(&net[i], bsort[i - a_size].data[j] + a_size);
    }
    for (int32_t j = 0; j < merge[i].size; ++j) {
      PushBack(&net[i], merge[i].data[j]);
    }
  }
}
SortNet OddEvenOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    SortNet sortnet = (SortNet)malloc(sizeof(struct MyVector));
    sortnet[0] = GetVector(0);
    return sortnet;
  }

  // Sort(k)
  SortNet* all_sortnets = (SortNet*)malloc((n + 1) * sizeof(SortNet));
  all_sortnets[0] = NULL;
  all_sortnets[1] = (SortNet)malloc(sizeof(struct MyVector));
  all_sortnets[1][0] = GetVector(0);

  // [k] -> Merge(ceil(k/2), floor(k/2))
  SortNet* all_mergenets = (SortNet*)malloc((n + 1) * sizeof(SortNet));
  all_mergenets[0] = NULL;
  all_mergenets[1] = (SortNet)malloc(sizeof(struct MyVector));
  all_mergenets[1][0] = GetVector(0);
  all_mergenets[2] = (SortNet)malloc(2 * sizeof(struct MyVector));
  all_mergenets[2][0] = GetVector(1);
  all_mergenets[2][1] = GetVector(1);
  PushBack(&all_mergenets[2][0], 1);
  PushBack(&all_mergenets[2][1], 0);

  for (short k = 3; k <= n; ++k) {
    OEmergenet(all_mergenets, k);
  }
  for (short k = 2; k <= n; ++k) {
    OEsortnet(all_sortnets, all_mergenets, k);
  }

  SortNet result = all_sortnets[n];
  for (short k = 1; k < n; ++k) {
    for (short j = 0; j < k; ++j) {
      free(all_sortnets[k][j].data);
      free(all_mergenets[k][j].data);
    }
    free(all_sortnets[k]);
    free(all_mergenets[k]);
  }
  for (short j = 0; j < n; ++j) {
    free(all_mergenets[n][j].data);
  }
  free(all_mergenets[n]);
  free(all_mergenets);
  free(all_sortnets);
  return result;
}

SortNet BitonicOrdering(short n) {
  if (n < 1) {
    return NULL;
  } else if (n == 1) {
    SortNet sortnet = (SortNet)malloc(sizeof(struct MyVector));
    sortnet[0] = GetVector(0);
    return sortnet;
  }
  short p = 0, nn = 1;
  while (nn < n) {
    p += 1;
    nn <<= 1;
  }
  SortNet sortnet = (SortNet)malloc(n * sizeof(struct MyVector));
  for (short i = 0; i < n; ++i) {
    sortnet[i] = GetVector((p * (p - 1)) >> 1);
  }
  for (short halflen = 1, len = 2; halflen < nn; halflen <<= 1, len <<= 1) {
    // First Block
    for (short start = 0; start < n; start += len) {
      for (short i = 0; i < halflen; ++i) {
        short j = len - i - 1;
        if (i + start < n && j + start < n) {
          PushBack(&sortnet[start + i], start + j);
          PushBack(&sortnet[start + j], start + i);
        }
      }
    }
    // Other Blocks
    for (short dist = halflen >> 1, delta = halflen; dist > 0; dist >>= 1, delta >>= 1) {
      for (short start = 0; start < n; start += delta) {
        for (short i = start; i < start + dist; ++i) {
          if (i + dist < n) {
            PushBack(&sortnet[i], i + dist);
            PushBack(&sortnet[i + dist], i);
          }
        }
      }
    }
  }
  return sortnet;
}

void FreeOrdering(struct MyVector* ordering, short n) {
  for (short i = 0; i < n; ++i) {
    free(ordering[i].data);
  }
  free(ordering);
}