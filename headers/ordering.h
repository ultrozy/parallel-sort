#ifndef ORDERING_H_
#define ORDERING_H_

#include <stdint.h>

struct MyVector {
  int32_t size;
  int32_t capacity;
  short* data;
};
typedef struct MyVector* SortNet;

struct MyVector GetVector(int32_t capacity);
void PushBack(struct MyVector* pvec, short value);

SortNet OddEvenOrdering(short n);
SortNet BitonicOrdering(short n);
void FreeOrdering(struct MyVector* ordering, short n);

#endif