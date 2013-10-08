#ifndef types_hh
#define types_hh

#include <TH1F.h>
#include <TH2F.h>

#include <halffloat.hh>

// @TODO use halffloat for grid data 

// typedef half* smallreal_array_t;
// typedef half  smallreal_t;

typedef float* real_array_t; 
typedef float  real_t;

typedef float* float_array_t; 
typedef float  float_t;

typedef int* integer_array_t;
typedef int  integer_t;

typedef TH1F hist1D;
typedef TH2F hist2D;

#endif 
