#ifndef grid_hh
#define grid_hh

#include "types.hh"

class igrid_base
{
public:
  igrid_base() {}
  virtual ~igrid_base() {}

  //
  // get the grid point coordinates
  //
  virtual const real_array_t* get_data( ) const = 0;

  //
  // get the grid point weights 
  //
  virtual const real_t* get_weights( ) const = 0;
  
  //
  // get number of grid points
  //
  virtual integer_t get_npoints( ) const = 0;

  //
  // get resource allocation
  //
  virtual integer_t get_resource_size( ) const = 0;

  //
  // get number of dimensions 
  //
  virtual integer_t get_ndim( ) const = 0;
  virtual integer_t get_ndim_effective( ) const = 0;
};
  
#endif
