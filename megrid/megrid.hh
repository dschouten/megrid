#ifndef megrid_hh
#define megrid_hh

#include <string>
#include <vector>

#ifndef __CINT__
#include "logger.hh"
#endif
#include "grid.hh"
#include "types.hh"
#include "transformations.hh"
#include "common.hh"

#define WC_DEFAULT (real_t)20.

/*
 *
 * @TODO
 *
 * X include weights for each grid point in AGF calculation
 *       X (lower priority) port the AGF calculations into own library ... do things better and avoid complex dependencies
 * X implement re-sampling with fast-cluster to compactify the grid
 * X export PDEFoam object
 * - test uniform normalization using empirical quantile function
 * - (lower priority) use boost::ndimarray instead of dynamic arrays (forego memory issues ...)
 * - (lower priority) allow to re-normalize the grid on the fly
 * - (lower priority) factorize the grid reading/storing into separate classes (ROOT, txt, HDF5 ...)
 *
 *
 *
 *
 */ 

class foam;

class megrid : public igrid_base
{
public:
  
  enum pdf_type { KNN = 0,
		  AGF = 1 };
  
  enum cluster_type { KMEANS       = 0,
		      HIERARCHICAL = 1 };
  
  megrid( std::string name, integer_t ndim ) ;
  
  megrid( std::string name,
	  std::string file_pattern,
	  std::string tree_name, 
	  std::string branch_names,
	  std::string weight_expression="", 
	  integer_t max_entries = -1 ) ;
  
  virtual ~megrid( );
  
  //
  // calculate the probability density for a given test point (x_{0}, ... x_{9}),
  // i.e., for grids with dimension <= 10
  // 
  real_t pdf( real_t  x0   , real_t  x1=0., real_t  x2=0., real_t  x3=0., real_t  x4=0., 
	      real_t  x5=0., real_t  x6=0., real_t  x7=0., real_t  x8=0., real_t  x9=0., 
	      real_t x10=0., real_t x11=0., real_t x12=0., real_t x13=0., real_t x14=0. ) const;
  
  //
  // calculate the probability density for a given test point vec=(x_{0}, ... x_{ndim-1})
  //
  real_t pdf( const std::vector<real_t>& vec ) const;
  
  //
  // get & set W_{c} value, equivalent of K in k-nn scheme 
  //
  real_t get_wc( ) const { return _wc; }
  void set_wc( real_t wc ) { _wc = wc > 0 ? (real_t)wc : WC_DEFAULT; }
  
  //
  // get & set max_{k}, max number of grid points to sample in PDF calculation
  //
  integer_t get_maxk( ) const { return _maxk; }
  void set_maxk( integer_t k ) { _maxk = k; }

  //
  // get & set PDF option
  //
  pdf_type get_option_pdf( ) const { return _option_pdf; }
  void set_option_pdf( pdf_type t ) { _option_pdf = t; }
  
  //
  // get the grid point coordinates
  //
  virtual const real_array_t* get_data( ) const { return _data; }
  
  //
  // get the grid point weights 
  //
  virtual const real_t* get_weights( ) const { return _wgts; }
  
  //
  // get number of grid points
  //
  virtual integer_t get_npoints( ) const { return _npoints; }
  
  //
  // get resource allocation
  //
  virtual integer_t get_resource_size( ) const { return _nreserved * sizeof(real_t) * _ndim; }
  
  //
  // get number of dimensions 
  //
  virtual integer_t get_ndim( ) const { return _ndim; }
  virtual integer_t get_ndim_effective( ) const { return _transformations.size() == 0 ? _ndim : _transformations.front()->get_in_dim(); }
  
  //
  // get & set name of the grid
  //
  std::string get_name( ) const { return _name; }
  void set_name( const std::string& name ) { _name = name; }
  
  //
  // set the lock flag
  //
  void set_locked( ) { _is_locked = true; }
  
  //
  // add point to the grid
  //
  bool add_point( const std::vector<real_t>& p, real_t w=1.0 );
  megrid& operator+=( const std::vector<real_t>& unweighted_pt );
  megrid& operator+=( const std::pair<std::vector<real_t>, real_t>& weighted_pt );
  
  //
  // apply a transformation to the grid coordinates
  // eg., a scaling or rotation
  //  
  bool apply_transformation( transform::itransformation_base* t, bool retain=true, bool force_configure=false );

  //
  // get list of transformations applied to this grid
  // 
  std::vector<const transform::itransformation_base*> get_transformations( ) const;
  
  //
  // resample the grid, iteratively merging closest grid points (brute-force method)
  //
  // niter is the number of merge operations (grid size is reduced by factor of <= 2^{niter})
  // maxd is the maximum distance between grid points that are allowed to be merged (<= 0 implies no maximum)
  //
  bool resample( unsigned niter, double maxd=-1 );
  
  // 
  // resample the grid using various clustering algorithms
  //
  bool cluster( cluster_type t, unsigned int nclusters, char method, char metric='e' );
  
  // 
  // create a 'foam' integration grid from this grid
  //
  foam* create_foam( real_t tail_cut=0.01, 
		     real_t vol_frac=1.0/30., 
		     unsigned ncells=1000, 
		     unsigned nmin=100,
		     unsigned nsampl=2000,
		     unsigned nbin=5 ) const;

  //
  // integrate the local p-value using transfer functions for observables
  // (i.e., grid coordinates)
  //
  // real_t integrate( const std::vector<real_t>&, const std::vector<itransfer_function_base*>& ) const;

  //
  // draw a projection of the grid
  //
  hist1D* project1D( integer_t ix, integer_t nbins, real_t a, real_t b, bool dyn=false ) const;
  hist2D* project2D( integer_t ix, integer_t iy, integer_t nbins, 
		     real_t ax, real_t bx, real_t ay, real_t by, bool dyn=false ) const;

  //
  // print a subset of the grid data to stdout
  //
  void scan( integer_t istart, integer_t nrows=50, const char* delim=" | " ) const;
  
  //
  // save grid into output file
  //
  // file_name is name of file
  // tree_name is the name of output tree
  // branch_names are the branches to associate to each dimension (separated by ':' delimiters)
  //
  bool store( const std::string& file_name, 
	      const std::string& tree_name, 
	      const std::string& branch_names ); 

  //
  // add another grid to this one
  //
  bool add_grid( const megrid& g );
  megrid& operator+=( const megrid& g );
  
  //
  // set the verbosity of the logging
  //
  static void set_loglevel( const std::string& lvl );

private:
  megrid() { }
  megrid( const megrid& ) { }
  megrid& operator=( const megrid& ) { return (*this); }
  
  //
  // load a pre-computed grid from a ROOT file
  //
  bool load( std::string file_pattern,
	     std::string tree_name,
	     std::string branch_names, 
	     std::string selection="",
	     integer_t max_entries = -1 );
  
  //
  // calculate the location of the point in the transformed grid
  // 
  std::vector<real_t> get_transformed_point( const std::vector<real_t>& x ) const;
  void get_transformed_point( const real_t* x, real_t* xp ) const;
  
  std::vector<transform::itransformation_base*> _transformations;
  
  std::string _name; // name of the grid
  real_t _wc; // W_{c} value used for calculating PDF's
  real_t _mindelta; // defines k in calculating PDF's such that w_{k} - w_{k+1} / (w_{k} + w_{k+1}) < mindelta, if pdf_type==MINDELTA
  integer_t _maxk; // defines k used in optimizing sigma and calculating PDF's if pdf_type==AGF
  
  pdf_type _option_pdf; // option for algorithm used in calculating PDF's
  
  bool _is_locked; // flag indicates grid is locked (no further operations can be performed)
  
  real_t* _wgts; // weights for each grid point
  real_array_t* _data; // grid data
  
  std::vector<real_t> _grid_dimensions; // length scale for each grid axis
  real_t _sum_wgts;

  integer_t _ndim; // dimensionality of the grid
  integer_t _npoints; // number of points in the grid
  integer_t _nreserved; // length of allocated data & weight arrays
  
  //
  // clear all internal data
  //
  void clear( );
  
  //
  // release memory for a given array
  //
  template<typename T>
  void release( T*& arr ) const;
  template<typename T>
  void release( T*& arr, integer_t len ) const;
  
  //
  // reserve space in memory
  //
  bool reserve( unsigned m );
  
  //
  // copy an array structure
  //
  real_array_t* copy( const real_array_t* arr, integer_t len, integer_t dim ) const;
  
  //
  // k-means clustering
  //
  bool kmeans_cluster( unsigned int nclusters, char method, char metric='e', unsigned npass=1 );
  
  //
  // hierarchical clustering
  //
  bool hierarchical_cluster( unsigned int nclusters, char method, char metric='e' );

  //
  // calculate grid dimensions
  //
  void set_metadata( );
};

#endif
