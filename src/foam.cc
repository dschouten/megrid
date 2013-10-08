
#include <cmath>
#include <limits>
#include <map>
#include <algorithm>
#include <stdexcept>

#include "foam.hh"
#include "kernels.hh"

#include "TMVA/PDEFoamEvent.h"
#include "TMVA/PDEFoamEventDensity.h"
#include "TMVA/Event.h"
#include "TFile.h"

typedef float float_t;

//////////////////////////////////////////////////////////////////////

foam::foam( const std::string& name,
	    const real_array_t* data, 
	    const real_t* weights,
	    unsigned ndim,
	    unsigned npoints,
	    double tail_frac,
	    double vol_frac,
	    unsigned ncells,
	    unsigned nmin,
	    unsigned nsampl,
	    unsigned nbin,
	    kernel_type t ) : _ndim( ndim ), _name( name ), _kernel( t )
{
  _pdefoam = new TMVA::PDEFoamEvent( name.c_str() );
  _pdedensity = new TMVA::PDEFoamEventDensity( get_box( data, weights, ndim, npoints, tail_frac, vol_frac ) );
  
  _pdefoam->SetDensity( _pdedensity );
  
  switch( logger::msg_level )
  {
    case log4cpp::Priority::DEBUG:
      _pdefoam->Log().SetMinType( TMVA::kDEBUG );
      break;
    case log4cpp::Priority::INFO:
      _pdefoam->Log().SetMinType( TMVA::kINFO );
      break;
    case log4cpp::Priority::WARN:
      _pdefoam->Log().SetMinType( TMVA::kWARNING );
      break;
    case log4cpp::Priority::ERROR:
      _pdefoam->Log().SetMinType( TMVA::kERROR );
      break;
    default:      
      _pdefoam->Log().SetMinType( TMVA::kINFO );
      break;
  }
  
  _pdefoam->SetDim(ndim);
  _pdefoam->SetnCells(ncells);    
  _pdefoam->SetnSampl(nsampl);    
  _pdefoam->SetnBin(nbin);      
  _pdefoam->SetNmin(nmin);
  _pdefoam->SetMaxDepth(0);
  
  _pdefoam->Initialize();
  
  for( unsigned idim=0; idim<ndim; idim++) 
  { 
    _pdefoam->SetXmin(idim, _xmin[idim]);
    _pdefoam->SetXmax(idim, _xmax[idim]);
  }
  
  std::vector<float_t> fl_buff( ndim, 0. );
  std::vector<TMVA::Event> ev_buff( 0 );      
  
  for( unsigned ipt=0; ipt<npoints; ++ipt )
  {    
    std::copy( data[ipt], data[ipt] + ndim, fl_buff.begin() );
    ev_buff.push_back( TMVA::Event( fl_buff, 0, weights[ipt] ) );
    _pdefoam->FillBinarySearchTree( &( ev_buff.back() ) );
  }
  
  logger::log() << msg::INFO << "creating PDEFoam from [" << ev_buff.size() 
		<< "] grid points in [" << fl_buff.size() << "] dimensions";
  
  _pdefoam->Create();
  
  for( unsigned ipt=0; ipt < npoints; ++ipt )
  {
    _pdefoam->FillFoamCells( &( ev_buff[ipt] ), weights[ipt] );
  }
  
  _pdefoam->Finalize();
  _pdefoam->DeleteBinarySearchTree();
}
//////////////////////////////////////////////////////////////////////

foam::foam( const std::string& filename, const std::string& name ) : _ndim( 0 ), _name( name ), _kernel( NONE )
{
  TFile* f = TFile::Open( filename.c_str(), "read" );
  if( f==NULL || f->IsZombie() )
  {
    logger::log() << msg::ERROR << "could not open file [" << filename << "]";
  }
  else
  {
    _pdefoam = dynamic_cast<TMVA::PDEFoamEvent*>( f->Get( name.c_str() ) );
    if( _pdefoam==0x0 )
    {
      logger::log() << msg::ERROR << "could not retrieve [" << name << "] from file [" << filename << "]";
    }
    else
    {
      _ndim = _pdefoam->GetTotDim();
      for( unsigned idim=0; idim<_ndim; ++idim )
      {
	_xmax.push_back( _pdefoam->GetXmax(idim) );
	_xmin.push_back( _pdefoam->GetXmin(idim) );
      }
    }
    _pdedensity = 0x0; // NB the density is not available anymore ... 
  }
}

//////////////////////////////////////////////////////////////////////

foam::~foam() 
{
  if( _pdefoam != 0x0 ) delete _pdefoam;
  if( _pdedensity != 0x0 ) delete _pdedensity;
}

//////////////////////////////////////////////////////////////////////

std::vector<double> foam::get_box( const real_array_t* grid,
				   const real_t* weights, 
				   unsigned ndim, 
				   unsigned npoints, 
				   double tail_frac, 
				   double vol_frac )
{
  _xmin.clear();
  _xmax.clear();
  
  real_t* xmin = new real_t[ndim];
  real_t* xmax = new real_t[ndim];
  
  for( unsigned idim=0; idim<ndim; ++idim) 
  {
    xmax[idim] = std::numeric_limits<double>::min();
    xmin[idim] = std::numeric_limits<double>::max();
  }
  
  unsigned ntail = (unsigned)(npoints * tail_frac / ndim);
  unsigned nbins = (unsigned)(1.0 / std::max( 0.1e-6, tail_frac ));
  
  for( unsigned ipt=0; ipt<npoints; ++ipt )
  {
    for( unsigned idim=0; idim<ndim; ++idim) 
    {
      if( grid[ipt][idim] < xmin[idim] )
	xmin[idim] = grid[ipt][idim];
      if( grid[ipt][idim] > xmax[idim] )
	xmax[idim] = grid[ipt][idim];
    }
  }
  
  TH1F* range_h = new TH1F[ndim]; 
  for( unsigned idim=0; idim<ndim; idim++ )
  {
    range_h[idim] = TH1F( Form("range%i", idim), "range", nbins, xmin[idim], xmax[idim] );
  }
  
  for( unsigned ipt=0; ipt<npoints; ++ipt )
  {
    for( unsigned idim=0; idim<ndim; ++idim) 
    {
      range_h[idim].Fill( grid[ipt][idim], weights[ipt] );
    }
  }
  
  for( unsigned idim=0; idim<ndim; idim++ ) 
  { 
    for( int i=1; i<(int)(nbins+1); ++i ) 
    {
      if (range_h[idim].Integral(0, i) > ntail)
      {
	xmin[idim]=range_h[idim].GetBinLowEdge(i);
	break;
      }
    }
    for( int i=nbins; i>0; --i )
    { 
      if (range_h[idim].Integral(i, (nbins+1)) > ntail)
      {
	xmax[idim]=range_h[idim].GetBinLowEdge(i+1);
	break;
      }
    }
    logger::log() << msg::INFO << "using range [" << xmin[idim] << "," << xmax[idim] << "] for dimension " << idim;
  }
  
  for( unsigned idim=0; idim<ndim; idim++ ) 
  { 
    _xmin.push_back(xmin[idim]);
    _xmax.push_back(xmax[idim]);
  }
  
  delete[] xmin;
  delete[] xmax;  
  delete[] range_h;
  
  std::vector<double> box;
  for( unsigned idim = 0; idim<_ndim; ++idim) 
  {
    box.push_back(std::fabs(_xmax[idim] - _xmin[idim]) * vol_frac);
  }
  
  return box;
}

//////////////////////////////////////////////////////////////////////

real_t foam::pdf( const std::vector<float_t>& vec ) const
{
  TMVA::PDEFoamKernelBase* kernel = 0x0; // new TMVA::PDEFoamKernelGauss( 0.05 ); 
  
  for( unsigned int idim=0; idim<_ndim; ++idim )
  {
    if( vec[idim] < _xmin[idim] || vec[idim] > _xmax[idim] )
    {
      logger::log() << msg::WARN << "PDF extrapolated from outside of stored phase space" ;
    }
  }
  if( _pdefoam != 0x0 )
  {
    return _pdefoam->GetCellValue( vec, TMVA::kValueDensity, kernel );
  }
  else
  {
    logger::log() << msg::ERROR << "foam not defined";
    return -1;
  }
}

//////////////////////////////////////////////////////////////////////

real_t foam::pdf( real_t  x0, real_t  x1, real_t  x2, real_t  x3, real_t  x4, 
		  real_t  x5, real_t  x6, real_t  x7, real_t  x8, real_t  x9,
		  real_t x10, real_t x11, real_t x12, real_t x13, real_t x14 ) const
{
  if( _ndim > 15 ) 
  {
    logger::log() << msg::ERROR << "too many dimensions, used pdf( vector<double> ) method instead" ;
    return -1;
  } 
  
  std::vector<float_t> vec;
  if( _ndim >= 1 ) vec.push_back( x0 );
  if( _ndim >= 2 ) vec.push_back( x1 );
  if( _ndim >= 3 ) vec.push_back( x2 );
  if( _ndim >= 4 ) vec.push_back( x3 );
  if( _ndim >= 5 ) vec.push_back( x4 );
  if( _ndim >= 6 ) vec.push_back( x5 );
  if( _ndim >= 7 ) vec.push_back( x6 );
  if( _ndim >= 8 ) vec.push_back( x7 );
  if( _ndim >= 9 ) vec.push_back( x8 );
  if( _ndim >= 10) vec.push_back( x9 );
  if( _ndim >= 11) vec.push_back( x10);
  if( _ndim >= 12) vec.push_back( x11);
  if( _ndim >= 13) vec.push_back( x12);
  if( _ndim >= 14) vec.push_back( x13);
  if( _ndim >= 15) vec.push_back( x14);
  return pdf( vec );
}

//////////////////////////////////////////////////////////////////////

const TMVA::PDEFoamEvent* foam::get_source( ) const 
{ 
  return _pdefoam; 
}

//////////////////////////////////////////////////////////////////////

const TMVA::PDEFoamEventDensity* foam::get_density( ) const
{
  return _pdedensity; 
}
