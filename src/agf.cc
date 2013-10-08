
//
// author: Doug Schouten <doug dot schouten at triumf dot ca>
//
// adapted from libAGF (by Peter Mills)
//

#include <algorithm>
#include <string>
#include <cmath>
#include <iterator>

#include <gsl/gsl_sf_gamma.h>

#include "agf.hh"
#include "logger.hh"
#include "numerical.hh"

////////////////////////////////////////////////////////////////////////////////

void dummy_gsl_handler( const char *reason, 
			const char * file,
			int line, 
			int gsl_errno ) 
{
  logger::log() << msg::ERROR << "GSL error in [" << file << "] @L" << line << " : " << reason;
}

////////////////////////////////////////////////////////////////////////////////

namespace agf
{
  real_t      configuration::_WEIGHTS_TOL	= 0.01;
  integer_t   configuration::_WEIGHTS_MAXITER	= 5000;
  real_t      configuration::MINFILTERW		= 1.0e-3;
  real_t      configuration::MAXFILTERW		= 1.0e+3;
  integer_t   configuration::KMIN               = 10;
  
  template <typename FLOAT_T>
  void deltaw( FLOAT_T var, void* param, FLOAT_T* wdiff, FLOAT_T* dWdv )
  {
    void **parcache = (void **)param;
    
    FLOAT_T W0;
    FLOAT_T W;
    integer_t k;
    FLOAT_T* dsqr;
    FLOAT_T* wts;
    FLOAT_T* clwts;
    
    W0    = *((FLOAT_T* ) parcache[0]);
    k     = *((integer_t*) parcache[1]);
    dsqr  = (FLOAT_T* ) parcache[2];
    wts   = (FLOAT_T* ) parcache[3];
    clwts = (FLOAT_T* ) parcache[4];
    
    W     = 0;
    *dWdv = 0;
    for( integer_t i=0; i<k; ++i )
    {
      wts[i] = std::exp(-dsqr[i]/var/2);
      wts[i] = wts[i] * clwts[i];
      (*dWdv) += dsqr[i] * wts[i] / var / var / 2;
      W += wts[i];
    }
    (*wdiff) = W-W0;
  } 
  
  ////////////////////////////////////////////////////////////////////////////////
  
  template <typename FLOAT_T>
  integer_t optimize( FLOAT_T* dsqr,   // distances squared
		      integer_t npts,  // number of distances
		      FLOAT_T Wc,      // objective total weight
		      FLOAT_T varlo,   // initial filter width
		      FLOAT_T varhi,  
		      FLOAT_T* weight, // returned weights
		      FLOAT_T& var_f,  // returned final filter variance
		      FLOAT_T* clwts ) // point weights    
  {
    long niter; // number of iterations to find optimal filter width
    long nmov;  // number of times bracket has been moved
    
    FLOAT_T wlo, whi, dWdvlo, dWdvhi;
    
    gsl_error_handler_t* old_handler;
    
    void* param[5];
    param[0] = &Wc;
    param[1] = &npts;
    param[2] = dsqr;
    param[3] = weight;
    param[4] = clwts;
    
    for( nmov=0; nmov<configuration::_WEIGHTS_MAXITER; ++nmov )
    {
      deltaw<FLOAT_T>( varlo, param, &wlo, &dWdvlo );
      if (wlo <= Wc) break;
      if (varlo <= configuration::MINFILTERW ) break;
      varlo = varlo/2;
      logger::log() << msg::WARN << "agf::optimize: lower bracket decreased to " << varlo;
    }
    for ( ; nmov<configuration::_WEIGHTS_MAXITER; ++nmov )
    {
      deltaw<FLOAT_T>( varhi, param, &whi, &dWdvhi );
      if (whi >= Wc) break;
      if (varhi >= configuration::MAXFILTERW ) break;
      varhi = varhi*2;
      logger::log() << msg::WARN << "agf::optimize: upper bracket increased to " << varhi;
    }
    
    // fprintf( stderr, "whi = %g, Wc = %g, wlo = %g", whi, Wc, wlo );
    
    if( wlo > Wc || whi < Wc )
    {
      logger::log() << msg::ERROR << "agf::optimize: failed to bracket root";
      logger::log() << msg::ERROR << "\t wlo(varlo=" << varlo << ")="<< wlo << "; whi(varhi=" << varhi << ")=" << whi;
      return -1;
    }
    
    old_handler = gsl_set_error_handler( &dummy_gsl_handler );
    
    var_f = supernewton( &deltaw<FLOAT_T>, (void *) param, varlo, varhi, (FLOAT_T)0.,
			 (FLOAT_T)configuration::_WEIGHTS_TOL, configuration::_WEIGHTS_MAXITER, niter, wlo, dWdvlo, whi, dWdvhi );
    
    // set error handler back to previous one:
    gsl_set_error_handler(old_handler);
    
    return niter+nmov;
  }
  
  namespace pdf 
  {
    
    ////////////////////////////////////////////////////////////////////////////////
    //
    // estimates the pdf at a test point based on a grid
    //
    // grid  : set of n vectors with m elements (i.e., the grid points)
    // clwts : the weight for each grid point 
    // vec   : the vector to classify
    // k     : the number of nearest neighbours to use in the classification
    // minw  : the minimum un-normalized sum of the weights with which to tune the filter width
    //
    // returns the esimated pdf using the adaptive Gaussian filter technique
    //
    template <typename FLOAT_T>
    FLOAT_T adaptive( FLOAT_T** grid, 
		      FLOAT_T* clwts, 
		      integer_t ndim, 
		      integer_t npts, 
		      FLOAT_T* vec,
		      FLOAT_T varlo, 
		      FLOAT_T varhi, 
		      integer_t k, 
		      FLOAT_T Wc, 
		      diagnostics* diag_param )
    {
      FLOAT_T var_f;		// final value of the filter width (as variance)
      FLOAT_T totw;		// total weight
      FLOAT_T* knearest;	// distances of k nearest neighbours
      FLOAT_T* clwtnearest;     // cluster/point weights of k nearest neighbours
      FLOAT_T* weight;		// the current value for the weights
      FLOAT_T scale, norm;	// normalisation coeff.
      FLOAT_T pdf;		// final calculated value of pdf
      
      // first we calculate all the distances:
      std::vector< std::pair<FLOAT_T, FLOAT_T> > dsqr;
      
      for( integer_t i=0; i<npts; ++i ) 
      {
	dsqr.push_back( std::pair<FLOAT_T,FLOAT_T>( metric(vec, grid[i], ndim), clwts[i] ) );
      }
      
      knearest    = new FLOAT_T[k]; 
      clwtnearest = new FLOAT_T[k];
      
      if( k < npts )
      {
	std::partial_sort( dsqr.begin(), dsqr.begin() + k, dsqr.end() ); 
      }

      for( integer_t i=0; i<k; ++i )
      {
	knearest[i] = dsqr[i].first;
	clwtnearest[i] = dsqr[i].second;
      }
      
      // calculate the weights using the central "engine":
      
      weight = new FLOAT_T[k];
      diag_param->nd = optimize(knearest, k, Wc, varlo, varhi, weight, var_f, clwtnearest);
      totw = 0;
      for( integer_t i=0; i<k; ++i ) 
      {
	totw+=weight[i];
      }
      
      //fprintf( stderr, "totw=%g, Wc=%g", totw, Wc );
      
      // use the final filter width to normalize the pdf:
      
      scale = std::sqrt(var_f*M_PI*2);
      norm = 1;
      for( integer_t i=0; i<ndim; ++i )
      {
	norm *= scale;
      }
      
      //norm = std::pow( var_f*M_PI*2, m/2. );
      //printf( "var_f=%g, totw=%g, norm=%g\n", var_f, totw, norm );
      
      pdf = totw/norm/npts;
      
      // set the diagnostic parameters:
      diag_param->f = weight[k-1]/weight[0];
      diag_param->W = totw;
      diag_param->V = var_f;
      
      delete [] knearest;
      delete [] clwtnearest;
      delete [] weight;
      
      return pdf;  
    }
    
    template real_t adaptive<real_t>( real_t** grid,
				      real_t *clwts, 
				      integer_t ndim, 
				      integer_t n, 
				      real_t *vec, 
				      real_t varlo,
				      real_t varhi, 
				      integer_t k, 
				      real_t Wc, 
				      diagnostics* diag_param );
    
    ////////////////////////////////////////////////////////////////////////////////
    //
    // estimates the pdf at a test point based on a grid
    //
    // grid  : set of n vectors with m elements (i.e., the grid points)
    // clwts : the weight for each grid point 
    // vec   : the vector to classify
    // k     : the number of nearest neighbours to use in the classification
    // minw  : the minimum un-normalized sum of the weights with which to tune the filter width
    //
    // returns the esimated pdf using the k-nearest neighbours
    //
    template <typename FLOAT_T>
    FLOAT_T knn( FLOAT_T** grid,
		 FLOAT_T* clwts,
		 integer_t ndim, 
		 integer_t npts,
		 FLOAT_T* vec,
		 integer_t k,
		 const FLOAT_T* dimensions )
    {
      FLOAT_T* knearest;		// distances of k nearest neighbours
      FLOAT_T* clwtnearest;             // cluster/point weights of k nearest neighbours
      FLOAT_T pdf;			// final calculated value of pdf

      long double totw;			// total weight
      long double volume;               // volume of n-ball containing the k-nearest points
      long double norm;                 // volume of entire grid

      std::vector< std::pair<FLOAT_T, FLOAT_T> > dsqr; // weighted distances to all grid points

      for( integer_t i=0; i<npts; ++i ) 
      {
	dsqr.push_back( std::pair<FLOAT_T,FLOAT_T>( metric(vec, grid[i], ndim), clwts[i] ) );
      }
      
      knearest    = new FLOAT_T[k]; 
      clwtnearest = new FLOAT_T[k];

      std::vector<FLOAT_T> ldim( ndim, 1 ); // length scale for each dimension
      if( dimensions == NULL )
      {
	for( integer_t i=0; i<ndim; ++i )
	{
	  ldim[i] = std::fabs( grid_helper<FLOAT_T>( grid, i, npts ).get_max() - 
			       grid_helper<FLOAT_T>( grid, i, npts ).get_min() );
	}
      }
      else
      {
	std::copy( dimensions, dimensions + ndim, ldim.begin() );
      }
      
      std::partial_sort( dsqr.begin(), dsqr.begin() + k, dsqr.end() ); 
      
      for( integer_t i=0; i<k; ++i )
      {
	knearest[i] = dsqr[i].first;
	clwtnearest[i] = dsqr[i].second;
      }
      
      // calculate the weights of the k-nearest neighbours
      totw = 0.;
      for( integer_t i=0; i<k; ++i ) 
      {
	totw+=clwtnearest[i];
      }
            
      // use the final filter width to normalize the pdf:
      
      volume = pow( M_PI, ndim/2.0 ) / gsl_sf_gamma( 1 + ndim/2.0 );
      norm = volume;
      for( integer_t i=0; i<ndim; ++i )
      {
	volume *= std::sqrt(knearest[k-1]); // volume of the n-ball containing nearest neighbours
	norm *= ldim[i]; // total volume of grid
      }
            
      pdf = (FLOAT_T)(totw/volume/norm);
      
      delete [] knearest;
      delete [] clwtnearest;
      
      return pdf;  
    }
    
    template real_t knn<real_t>( real_t** grid, 
				 real_t *clwts, 
				 integer_t ndim, 
				 integer_t n, 
				 real_t *vec, 
				 integer_t k, 
				 const real_t* dimensions );
  }
  
}
