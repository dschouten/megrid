#include "transformations.hh"
#include "grid.hh"
#include "logger.hh"

#include <algorithm>
#include <vector>
#include <string>
#include <cmath>
#include <limits>

#include <boost/tokenizer.hpp>

using std::vector;
using std::string;

using namespace transform;

//////////////////////////////////////////////////////////////////////////////// 

////////////////////////////////////////////////////////////////////////////////

//
// specify a vector rotation by identifying the coordinates
// of the input vectors
//
// possible eg. is 'X0:Y0:X1:X2:Y1:Y2'
//
// the vector defined by (X0,Y0) defines the rotation, which is
// then applied to the other 2D vectors (X1,Y1), (X2,Y2), &c.
//
// z-coordinates are, of course, left untouched
//

bool rotate_phi::configure(const igrid_base* g, bool force)
{
  if(_configured && !force)
  {
    return g->get_ndim() == (integer_t)get_in_dim();
  }

  _configured = false;

  _i_pairs.clear();
  _paired_indexes.clear();

  std::transform(_coord_def.begin(), _coord_def.end(), _coord_def.begin(), ::tolower);

  vector<string> coord_ids;  
  boost::char_separator<char> sep(":");
  boost::tokenizer< boost::char_separator<char> > tokens(_coord_def, sep);
  std::copy(tokens.begin(), tokens.end(), std::back_inserter(coord_ids));
  
  bool flag = true;
  
  // this array of pairs stores the i,j that are the coordinate indexes for
  // the x,y for each vector

  vector< std::pair<int,int> > 
    point_coord_indices(coord_ids.size(), std::pair<int,int>(-1,-1));
  
  for(unsigned icoord=0; icoord<coord_ids.size(); ++icoord)
  {
    if(coord_ids[icoord].size() != 2 || (coord_ids[icoord][0] != 'x' && 
					   coord_ids[icoord][0] != 'y' && 
					   coord_ids[icoord][0] != 'z')) 
    {
      logger::log() << msg::ERROR << "coordinate [" << coord_ids[icoord][0] << "] is not understood";
      flag = false; // only cartesian coordinate systems allowed
      break; 
    }
    int ip = atoi(&(coord_ids[icoord][1]));
    if(ip < 0 || ip > (int) coord_ids.size())
    {
      logger::log() << msg::ERROR << "index [" << ip << "] is out of bounds";
      flag = false;
      break;
    }
    if(coord_ids[icoord][0] == 'x') // cartesian coordinate this token is for
      point_coord_indices[ip].first = icoord;
    else if(coord_ids[icoord][0] == 'y')
      point_coord_indices[ip].second = icoord;
  }
  
  for(unsigned ip=0; ip<point_coord_indices.size(); ++ip)
  {
    if(point_coord_indices[ip].first >= 0 && point_coord_indices[ip].second >= 0) // valid object coordinates
    {
      _i_pairs.push_back(std::pair<unsigned,unsigned>(point_coord_indices[ip].first,
							point_coord_indices[ip].second));
      _paired_indexes.insert(point_coord_indices[ip].first );
      _paired_indexes.insert(point_coord_indices[ip].second);
      logger::log() << msg::DEBUG << "adding vector " << point_coord_indices[ip].first << "," << point_coord_indices[ip].second;
    }
  }
  
  // _i_pairs now contains the grid dimension number for the x,y pairs of each vector;
  // the rotation proceeds by setting the 0'th vector along the y-axis and rotating all the
  // others accordingly
  
  if(flag)
  {
    set_in_dim(g->get_ndim());
    set_out_dim(g->get_ndim() - 1);
    
    _configured = true;
  }  

  return (_configured && g->get_ndim() == (integer_t)get_in_dim());
}

////////////////////////////////////////////////////////////////////////////////

real_t* rotate_phi::operator()(const real_t* x) const
{
  if(!_configured) 
  {
    return NULL;
  }
  real_t* xp = new real_t[get_out_dim()];
  (*this)(x, xp);
  return xp;
}

////////////////////////////////////////////////////////////////////////////////

void rotate_phi::operator()(const real_t* x, real_t* xp) const
{  
  if(!_configured || _i_pairs.size() == 0)
  {
    std::copy(x, x+get_out_dim(), xp);
  }
  else
  {
    double X0 = x[ _i_pairs[0].first  ];
    double Y0 = x[ _i_pairs[0].second ]; // (x,y) for 'reference' vector
    
    double ctheta = X0 / std::sqrt(pow(X0,2) + pow(Y0,2));
    double stheta = Y0 / std::sqrt(pow(X0,2) + pow(Y0,2));
        
    for(unsigned icoord=0; icoord<get_in_dim(); ++icoord)
    {
      xp[ icoord > _i_pairs[0].first ? icoord-1 : icoord ] = x[icoord];
    }
    
    for(unsigned ip=0; ip<_i_pairs.size(); ++ip)
    {
      X0 = x[ _i_pairs[ip].first  ];
      Y0 = x[ _i_pairs[ip].second ];
      // NB: shift the indexes if to the right of the redundant coordinate X0
      xp[ _i_pairs[ip].first   > _i_pairs[0].first ? _i_pairs[ip].first-1  : _i_pairs[ip].first   ] = X0*ctheta - Y0*stheta; // x'
      xp[ _i_pairs[ip].second  > _i_pairs[0].first ? _i_pairs[ip].second-1 : _i_pairs[ip].second  ] = X0*stheta + Y0*ctheta; // y'
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool rotate_phi::operator==(const itransformation_base* t) const
{
  const rotate_phi* other = dynamic_cast<const rotate_phi*>(t);
  if(other != NULL)
  {
    return (other->get_vector_indices() == get_vector_indices());
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////// 

////////////////////////////////////////////////////////////////////////////////

bool norm_gauss::configure(const igrid_base* g, bool force) 
{
  if(_configured && !force)
  {
    return g->get_ndim() == (integer_t)get_in_dim();
  }

  _configured = false;
  real_t sumw=0;
  for(integer_t ipt=0; ipt<g->get_npoints(); ++ipt)
  {
    sumw += g->get_weights()[ipt];
    for(unsigned idim=0; idim<get_in_dim(); ++idim)
    {
      _means[idim] += g->get_data()[ipt][idim];
      _std_deviations[idim] += std::pow(g->get_data()[ipt][idim], 2);
    }
  }
  if(sumw > std::numeric_limits<real_t>::epsilon())
  {
    for(unsigned idim=0; idim<get_in_dim(); ++idim)
    {
      _means[idim] /= sumw;
      _std_deviations[idim] = std::sqrt(_std_deviations[idim] / sumw);
    }
    _configured = true;
  }
  return (_configured && g->get_ndim() == (integer_t)get_in_dim());
}

////////////////////////////////////////////////////////////////////////////////

real_t* norm_gauss::operator()(const real_t* x) const
{
  if(!_configured) 
  {
    return NULL;
  }
  real_t* xp = new real_t[get_out_dim()];
  (*this)(x, xp);
  return xp;
}

////////////////////////////////////////////////////////////////////////////////

void norm_gauss::operator()(const real_t* x, real_t* xp) const
{
  if(_configured)
  {
    for(unsigned idim=0; idim<get_in_dim(); ++idim)
    {
      xp[idim] = (x[idim] - _means[idim]) / _std_deviations[idim];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool norm_gauss::operator==(const itransformation_base* t) const
{
  const norm_gauss* other = dynamic_cast<const norm_gauss*>(t);
  if(other != NULL)
  {
    return (other->get_std_deviations() == _std_deviations &&
	     other->get_means() == _means);
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

real_t* norm_range::operator()(const real_t* x) const
{
  if(!_configured) 
  {
    return NULL;
  }
  real_t* xp = new real_t[get_out_dim()];
  (*this)(x, xp);
  return xp;
}

////////////////////////////////////////////////////////////////////////////////

void norm_range::operator()(const real_t* x, real_t* xp) const
{
  if(_configured)
  {
    integer_t l, m, r, d;
    
    for(unsigned idim=0; idim<get_in_dim(); ++idim)
    {
      l = 0;
      m = _sorted_data.size() / 2;
      r = _sorted_data.size();
      
      // binary search to find index of coordinate
      
      // @TODO interpolate between bracketing points ... 
      
      while(true)
      {
	if(x[idim] > _sorted_data[idim][m])
	{
	  d = std::abs(r - m) / 2;
	  l = m;
	  m += d;
	  if(d <= 1)
	  {
	    break;
	  }
	}
	if(x[idim] < _sorted_data[idim][m])
	{
	  d = std::abs(m - l) / 2;
	  r = m;
	  m -= d;
	  if(d <= 1)
	  {
	    break;
	  }
	}
	if(x[idim] == _sorted_data[idim][m])
	  break;
      }      
      xp[idim] = ((real_t)m) / _sorted_data.size();
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

bool norm_range::configure(const igrid_base* g, bool force)
{
  if(_configured && !force)
  {
    return g->get_ndim() == (integer_t)get_in_dim();
  }
  _configured = false;
  _sorted_data.clear();
  for(integer_t ipt=0; ipt<g->get_npoints(); ++ipt)
  {
    if(ipt % _sampling_frequency == 0)
    {
      for(unsigned idim=0; idim<get_in_dim(); ++idim)
      {
	_sorted_data[idim].push_back(g->get_data()[ipt][idim]);
      }
    }
  } 
  for(unsigned idim=0; idim<get_in_dim(); ++idim)
  {
    std::sort(_sorted_data[idim].begin(), _sorted_data[idim].end());
  }
  if(_sorted_data[0].size() > 1)
  {
    _configured = true;
  }
  return (_configured && g->get_ndim() == (integer_t)get_in_dim());
}

////////////////////////////////////////////////////////////////////////////////

bool norm_range::operator==(const itransformation_base* t) const
{
  const norm_range* other = dynamic_cast<const norm_range*>(t);
  if(other != NULL)
  {
    return (other->get_cdf() == _sorted_data);
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

real_t* scale_fixed::operator()(const real_t* x) const
{
  real_t* xp = new real_t[get_out_dim()];
  (*this)(x, xp);
  return xp;
}

////////////////////////////////////////////////////////////////////////////////

void scale_fixed::operator()(const real_t* x, real_t* xp) const
{
  for(unsigned int idim=0; idim<get_in_dim(); ++idim)
  {
    xp[idim] = x[idim] * _scales[idim];
  }
}

////////////////////////////////////////////////////////////////////////////////

bool scale_fixed::configure(const igrid_base* g, bool force)
{
  return (g->get_ndim() == (integer_t)get_in_dim());
}

////////////////////////////////////////////////////////////////////////////////

bool scale_fixed::operator==(const itransformation_base* t) const
{
  const scale_fixed* other = dynamic_cast<const scale_fixed*>(t);
  if(other != NULL)
  {
    return (other->get_scale_factors() == _scales);
  }
  return false;
}
