#ifndef agf_hh
#define agf_hh

#include "types.hh"

// #include <boost/function.hpp>

// @FIXME set these as configurable parameters
// #define WEIGHTS_TOL 0.01
// #define WEIGHTS_MAXITER 5000
// #define MINFILTERW 1.0e-3
// #define MAXFILTERW 1.0e+3

namespace agf
{
  class configuration
  {
  public:
    static real_t  _WEIGHTS_TOL;
    static integer_t _WEIGHTS_MAXITER;
    static real_t  MINFILTERW;
    static real_t  MAXFILTERW;
    static integer_t KMIN;
  };

  struct diagnostics
  {
    integer_t nd;	// number of iterations
    real_t f;		// ratio of min. weight to max.
    real_t W;		// total weight
    real_t V;           // final filter variance
  };

  // namespace metrics // @TODO allow metric plugins ... just force Euclidean for now 
  // {
  //   enum emetric_t {
  //     EUCLID = 0,    
  //   };
  //
  //   template<typename FLOAT_T>
  //   boost::function<void ( FLOAT_T* x, FLOAT_T* y, integer_t ndim )> get_metric( emetric_t t ) { /* @TODO */ }
  //
  //   template <typename FLOAT_T>
  //   void euclidean( FLOAT_T* x,
  // 		       FLOAT_T* y, 
  // 		       integer_t ndim );
  // }

  template <typename FLOAT_T>
  FLOAT_T metric( FLOAT_T* x,
		  FLOAT_T* y, 
		  integer_t ndim )
  {
    FLOAT_T d=0;
    for( integer_t idim=0; idim<ndim; ++idim )
    {
      d += std::pow( x[idim] - y[idim], 2 );
    }
    return d;
  }
  

  template <typename FLOAT_T>
  void deltaw( FLOAT_T var,      // filter width
	       void* param,      // parameter array
	       FLOAT_T* wdiff,   // function value
	       FLOAT_T* dWdv );  // function derivative 

   template <typename FLOAT_T>
   integer_t optimize( FLOAT_T* dsqr,    // distances squared
		       integer_t npts,   // number of distances
		       FLOAT_T Wc,       // objective total weight
		       FLOAT_T varlo,    // initial filter width
		       FLOAT_T varhi,  
		       FLOAT_T* weight,  // returned weights
		       FLOAT_T& var_f,   // returned final filter variance
		       FLOAT_T* clwts ); // point weights    
  
  namespace pdf
  {
    template <typename FLOAT_T>
    class grid_helper
    {
    public:
      grid_helper( FLOAT_T **grid, integer_t idim, integer_t npts ) :
	_grid( grid ),
	_idim( idim ),
	_npts( npts ) { }
      
      FLOAT_T operator[]( integer_t ipt ) const { return _grid[ipt][_idim]; }

      FLOAT_T get_max( ) const 
      {
	FLOAT_T m = std::numeric_limits<FLOAT_T>::min();
	for( integer_t ipt=0; ipt<_npts; ++ipt )
	{
	  if( (*this)[ipt] > m ) m = (*this)[ipt];
	}
	return m;
      }

      FLOAT_T get_min( ) const 
      {
	FLOAT_T m = std::numeric_limits<FLOAT_T>::max();
	for( integer_t ipt=0; ipt<_npts; ++ipt )
	{
	  if( (*this)[ipt] < m ) m = (*this)[ipt];
	}
	return m;
      }
      
      FLOAT_T** _grid;
      integer_t _idim;
      integer_t _npts;
    };

    template <typename FLOAT_T>
    FLOAT_T adaptive( FLOAT_T** grid,	// array of grid points
		      FLOAT_T* clwts,	// weights for each point
		      integer_t ndim,   // number of dimensions in the grid
		      integer_t n,      // number of grid points
		      FLOAT_T* vec,	// vector (test point) at which to estimate the grid density
		      FLOAT_T varlo,    // lower filter width
		      FLOAT_T varhi,    // upper filter width
		      integer_t k,      // number of grid points to use to determine filter width
		      FLOAT_T Wc,       // weighted sum (equivalent to k in kNN)
		      diagnostics* diag_param ); 

    template <typename FLOAT_T>
    FLOAT_T knn( FLOAT_T** grid,		// array of grid points
		 FLOAT_T* clwts,		// weights for each point
		 integer_t ndim,		// number of dimensions in the grid
		 integer_t n,			// number of grid points
		 FLOAT_T* vec,			// vector (test point) at which to estimate the grid density
		 integer_t k,			// number of nearest neighbours
		 const FLOAT_T* dimensions=NULL );	// grid dimensions (for volume normalization)

  }
}

#endif
