#ifndef foam_hh
#define foam_hh

#include <string>
#include <vector>

#ifndef __CINT__
#include "logger.hh"
#endif
#include "types.hh"

namespace TMVA
{
  class PDEFoamEvent;
  class PDEFoamEventDensity;
}

class foam
{
private:
  foam() { }
  foam( const foam& ) { }
  foam& operator=( const foam& ) { return (*this); }
  
public:  
  typedef float float_t;
  
  enum kernel_type { NONE = 0,
		     GAUSS = 1,
		     AGF = 2 };
  
  //
  // create a new foam from an array of data
  //
  foam( const std::string& name, 
	const real_array_t* data, 
	const real_t* weights,
	unsigned ndim,
	unsigned npoints,
	double tail_frac=0.01, 
	double vol_frac=1.0/30., 
	unsigned ncells=1000, 
	unsigned nmin=100,
	unsigned nsampl=2000,
	unsigned nbin=5,
	kernel_type t=NONE );
  
  //
  // load the foam back from a ROOT file
  //
  foam( const std::string& filename, const std::string& name );
  
  ~foam();
  
  std::string get_name() const { return _name; }
  
  //
  // calculate the probability density for a given test point (x_{0}, ... x_{14}),
  // i.e., for grids with dimension <= 15
  // 
  real_t pdf( real_t  x0   , real_t  x1=0., real_t  x2=0., real_t  x3=0., real_t  x4=0., 
	      real_t  x5=0., real_t  x6=0., real_t  x7=0., real_t  x8=0., real_t  x9=0.,
	      real_t x10=0., real_t x11=0., real_t x12=0., real_t x13=0., real_t x14=0 ) const;
  
  //
  // calculate the probability density for a given test point vec=(x_{0}, ... x_{ndim-1})
  //
  real_t pdf( const std::vector<float_t>& vec ) const;
  
  //
  // set the kernel density estimator
  // 
  void set_kernel_type( kernel_type t=NONE ) { /* @TODO */ };
  
  //
  // store the foam in a ROOT format file
  //
  void save( const std::string& filename ) { /* @TODO */ };
  
  //
  // access the PDEFoam objects directly
  //
  const TMVA::PDEFoamEvent* get_source( ) const ;
  const TMVA::PDEFoamEventDensity* get_density( ) const ;
  
private:  
  std::vector<double> _xmin; // min value along each dimension
  std::vector<double> _xmax; // max value along each dimension
  
  //
  // calculate the dimensions of the box that contains fraction 
  // of the total phase space, store ranges for each dimension
  //
  std::vector<double> get_box( const real_array_t* grid,
			       const real_t* weights, 
			       unsigned ndim, 
			       unsigned npoints, 
			       double tail_frac,
			       double vol_frac );
  
  TMVA::PDEFoamEvent* _pdefoam; // pointer to the PDEFoam object
  TMVA::PDEFoamEventDensity* _pdedensity; // pointer to the PDEFoam density object
  
  unsigned _ndim; // dimensionality of the foam
  
  std::string _name; // name of the foam
  
  kernel_type _kernel; // type of kernel to use when sampling the foam
};

#endif
