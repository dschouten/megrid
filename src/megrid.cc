#include "megrid.hh"
#include "foam.hh"

#include <list>
#include <limits>
#include <memory>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <stdexcept>

#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeFormula.h>
#include <TROOT.h>
#include <TDirectory.h>
#include <TEventList.h>

#include "agf.hh"
#include "cluster.hh" 

// #include "agf_lib.h"

#define INIT_RESERVE 1000

using std::string;
using std::vector;
using std::pair;
using std::list;

//////////////////////////////////////////////////////////////////////

megrid::megrid(string name, integer_t ndim) : 
  igrid_base(),
  _name (name),
  _wc (WC_DEFAULT),
  _mindelta(-1),
  _maxk(0),
  _option_pdf(AGF),
  _is_locked (false),
  _wgts (NULL),
  _data (NULL),
  _grid_dimensions(ndim, 1),
  _sum_wgts(0.),
  _ndim(ndim),
  _npoints(0),
  _nreserved(0)
{               
  
  
}

//////////////////////////////////////////////////////////////////////

megrid::megrid(string name,
	       string file_pattern, 
	       string tree_name, 
	       string branch_names, 
	       string weight_expression,
	       integer_t max_entries) :
  igrid_base(),
  _name (name),
  _wc (WC_DEFAULT),
  _mindelta(-1),
  _maxk(0),
  _option_pdf(AGF),
  _is_locked (false),
  _wgts (NULL),
  _data (NULL),
  _grid_dimensions(0),
  _sum_wgts(0.),
  _ndim(0),
  _npoints(0),
  _nreserved(0)
{ 
  if(!load(file_pattern, tree_name, branch_names, weight_expression, max_entries)) clear(); 
}

//////////////////////////////////////////////////////////////////////

megrid::~megrid() 
{
  clear();
}

//////////////////////////////////////////////////////////////////////

bool megrid::load(string file_pattern,
		  string tree_name,
		  string branch_names, 
		  string weight_expression,
		  integer_t max_entries)
{
  clear();
  
  // boost::shared_ptr<TFile> fptr(new TFile(file_name.c_str(), "read"));
  // if(fptr->IsZombie())
  // {
  //   logger::log() << msg::ERROR << "could not open file [" << file_name << "] for reading" ;
  //   return false;
  // }
  
  // TTree* tptr = dynamic_cast<TTree*>(fptr->Get(tree_name.c_str()));
  // if(!tptr)
  // {
  //   logger::log() << msg::ERROR << "could not retrieve tree [" << tree_name << "] from input file" ;
  //   return false;
  // }
  
  boost::shared_ptr<TChain> tptr(new TChain(tree_name.c_str()));
  tptr->Add(file_pattern.c_str());  
  tptr->SetBranchStatus("*", true);
  
  logger::log() << msg::INFO << "begin reading in grid data" ;  
  
  if(weight_expression.size() == 0)
  {
    weight_expression = "1.00";
  }
  
  max_entries = min(max_entries, tptr->GetEntries());
  if(max_entries <= 0)
  {
    max_entries = tptr->GetEntries();
  }
  reserve(max_entries + INIT_RESERVE); 
  
  // _npoints = min(tptr->Draw("", weight_expression.c_str(), "goff"), max_entries);
  
  logger::log() << msg::INFO << "will read in maximum [" << max_entries << "] points" ;
  
  vector<string> dim_exprs;  
  boost::char_separator<char> sep("|");
  boost::tokenizer< boost::char_separator<char> > tokens(branch_names, sep);
  std::copy(tokens.begin(), tokens.end(), std::back_inserter(dim_exprs));     
  _ndim = dim_exprs.size();
  vector< boost::shared_ptr<TTreeFormula> > dim_formulas;  
  for(unsigned ibr = 0; ibr < (unsigned) _ndim; ++ibr)
  {
    logger::log() << msg::INFO << "set dimension " << ibr << " to [" << dim_exprs[ibr] << "]" ;
    dim_formulas.push_back(boost::shared_ptr<TTreeFormula>(new TTreeFormula(Form("br_%i", ibr), dim_exprs[ibr].c_str(), tptr.get())));
    dim_formulas.back()->SetQuickLoad(true);
  }
  boost::shared_ptr<TTreeFormula> wgt_formula(new TTreeFormula("weight", weight_expression.c_str(), tptr.get()));
  wgt_formula->SetQuickLoad(true);
  
  integer_t ipt = 0;
  for(integer_t ientry = 0; ientry < max_entries; ++ientry)
  {    
    int nb = tptr->GetEntry((int)ientry);
    if(nb < 0)
    {
      logger::log() << msg::ERROR << "error reading entry #" << ientry << " in ntuple" ;
      return false;
    }
    
    if(max_entries < 10 || ientry % (max_entries/10) == 0)
    {
      logger::log() << msg::INFO << "processed [" << ientry << "] entries in input ntuple table" ;
      logger::log() << msg::INFO << "stored [" << ipt << "] grid points" ;
    }
    
    if(wgt_formula->EvalInstance(0) == 0)
    {
      continue;
    }
    
    _wgts[ipt] = wgt_formula->EvalInstance(0);
    
    _data[ipt] = new real_t[_ndim];
    for(integer_t idim = 0; idim < _ndim; ++idim)
    {
      _data[ipt][idim] = dim_formulas[idim]->EvalInstance(0);
    }
    ipt += 1;
  }
  
  _npoints = ipt;
  
  set_metadata();
  
  logger::log() << msg::INFO << "stored [" << ipt << "] grid points in total" ;
  logger::log() << msg::INFO << "successfully loaded " << _ndim << "-dimensional grid from [" << file_pattern << "]";
  
  return true;
}

//////////////////////////////////////////////////////////////////////

real_t megrid::pdf(real_t  x0, real_t  x1, real_t  x2, real_t  x3, real_t  x4, 
		   real_t  x5, real_t  x6, real_t  x7, real_t  x8, real_t  x9, 
		   real_t x10, real_t x11, real_t x12, real_t x13, real_t x14) const
{
  if(_ndim > 15) 
  {
    logger::log() << msg::ERROR << "too many dimensions, used pdf(vector<double>) method instead" ;
    return -1;
  } 
  
  integer_t _ndim_in = get_ndim_effective();
  
  vector<real_t> vec; 
  if(_ndim_in >= 1) vec.push_back(x0);
  if(_ndim_in >= 2) vec.push_back(x1);
  if(_ndim_in >= 3) vec.push_back(x2);
  if(_ndim_in >= 4) vec.push_back(x3);
  if(_ndim_in >= 5) vec.push_back(x4);
  if(_ndim_in >= 6) vec.push_back(x5);
  if(_ndim_in >= 7) vec.push_back(x6);
  if(_ndim_in >= 8) vec.push_back(x7);
  if(_ndim_in >= 9) vec.push_back(x8);
  if(_ndim_in >= 10) vec.push_back(x9);
  if(_ndim_in >= 11) vec.push_back(x10);
  if(_ndim_in >= 12) vec.push_back(x11);
  if(_ndim_in >= 13) vec.push_back(x12);
  if(_ndim_in >= 14) vec.push_back(x13);
  if(_ndim_in >= 15) vec.push_back(x14);
  return pdf(vec);
}

//////////////////////////////////////////////////////////////////////

real_t megrid::pdf(const vector<real_t>& vec) const
{
  if(_npoints == 0 || _data == NULL || _wgts == NULL)
  {
    logger::log() << msg::ERROR << "grid not loaded" ;
    return -1;
  }
  
  vector<real_t> xarr;
  std::copy(vec.begin(), vec.end(), back_inserter(xarr));
  
  if((integer_t) xarr.size() != get_ndim_effective())
  {
    logger::log() << msg::ERROR << "dimension mismatch in test point" ;
    return -1;
  }    
  
  xarr = get_transformed_point(xarr);
  
  if((integer_t) xarr.size() != get_ndim())
  {
    logger::log() << msg::ERROR << "dimension mismatch in test point" ;
    return -1;
  }
  
  real_t* ptr_vec = &(xarr[0]);
  
  integer_t kmax = (_maxk == 0 ? _npoints : _maxk);
  
  if(_option_pdf == AGF)
  {    
    real_t d=0;
    real_t w=0;
    for(integer_t ipt=0; ipt<_npoints; ++ipt)
    {
      d += agf::metric(_data[ipt], ptr_vec, _ndim); 
      w += _wgts[ipt];
    }
    if(w > 0)
      d /= w; // total variance
    else
      return -1;
    
    double vlo = d / pow(_npoints, 2.0/_ndim) / 32.0; // NB the division factor is empirical
    double vhi = d; 
    
    logger::log() << msg::DEBUG << "using variance bounds [" << vlo << "," << vhi << "]";
    
    agf::diagnostics diag_params;
    
    if(_mindelta < 0)
    {
      double pdf = agf::pdf::adaptive<real_t>(_data, _wgts, _ndim, _npoints, ptr_vec, vlo, vhi, kmax, _wc, &diag_params);
      logger::log() << msg::DEBUG << "filter width: " << diag_params.V << ", total W: " << diag_params.W;
      
      // real_t varbnds[2] = { vlo, vhi };
      // agf_diag_param dpar;
      // real_t pdf = agf_calc_pdf<real_t>(_data, _ndim, _npoints, ptr_vec, varbnds, kmax, _wc, &dpar); 
      
      return pdf;
    }
    else
    {
      std::vector<real_t> testpdf(kmax, -1);
      for(integer_t k=0; k<kmax; ++k)
      {
      	testpdf[k] = agf::pdf::adaptive<real_t>(_data, _wgts, _ndim, _npoints, ptr_vec, vlo, vhi, k + agf::configuration::KMIN, _wc, &diag_params);
      	if(k > 0)
      	{
      	  if(std::fabs(testpdf[k] - testpdf[k-1]) / (testpdf[k] + testpdf[k-1]) < _mindelta)
      	    return testpdf[k];
      	}
      }
      logger::log() << msg::INFO << "failed to find k within bounds [" << agf::configuration::KMIN << "," << kmax << "]";
      return testpdf.back();
    }
  }
  else
  {    
    if(_mindelta < 0)
    {
      return agf::pdf::knn<real_t>(_data, _wgts, _ndim, _npoints, ptr_vec, kmax, &(_grid_dimensions[0]));
    }
    else
    {
      std::vector<real_t> testpdf(kmax, -1);
      for(integer_t k=0; k<kmax; ++k)
      {
	testpdf[k] = agf::pdf::knn<real_t>(_data, _wgts, _ndim, _npoints, ptr_vec, k + agf::configuration::KMIN, &(_grid_dimensions[0]));
	if(k > 0)
	{
	  if(std::fabs(testpdf[k] - testpdf[k-1]) / (testpdf[k] + testpdf[k-1]) < _mindelta)
	    return testpdf[k];
	}
      }
      logger::log() << msg::INFO << "failed to find k within bounds [" << agf::configuration::KMIN << "," << kmax << "]";
      return testpdf.back();
    }
  }
}

//////////////////////////////////////////////////////////////////////

void megrid::clear() 
{
  release<real_array_t>(_data, _nreserved);
  release<real_t>(_wgts);
  
  _npoints   = 0;
  _nreserved = 0;
  
  _is_locked = false;
}

//////////////////////////////////////////////////////////////////////

template<typename T>
void megrid::release(T*& arr, integer_t len) const
{
  // logger::log() << msg::DEBUG << "clearing array @ " << arr << " of size [" << len << "]";
  if(arr != 0x0)
  {
    for(integer_t irow = 0; irow < len; ++irow) 
    {
      // std::cout << irow << std::endl;
      delete[] arr[irow]; 
    }
    delete[] arr;
    // logger::log() << msg::DEBUG << "cleared grid with space for [" << len << "] points";
  }
  arr = 0x0;
}

//////////////////////////////////////////////////////////////////////

template<typename T>
void megrid::release(T*& arr) const
{
  if(arr != 0x0)
  {
    delete[] arr;
  }
  arr = 0x0;
}

//////////////////////////////////////////////////////////////////////

real_array_t* megrid::copy(const real_array_t* arr, integer_t len, integer_t dim) const
{
  if(arr != 0x0)
  {
    real_array_t* clone = new real_array_t[len]; 
    for(integer_t irow = 0; irow < len; ++irow) 
    {
      clone[irow] = new real_t[dim];
      std::copy(arr[irow], arr[irow] + dim, clone[irow]);
    }
    return clone;
  }
  return 0x0;
}

//////////////////////////////////////////////////////////////////////

bool megrid::resample(unsigned int nsample, double dmax)
{
  // if(_is_locked)
  // {
  //   logger::log() << msg::ERROR << "grid is locked, cannot resample";
  //   return false;
  // }
  
  if(_npoints < (integer_t)nsample)
  {
    logger::log() << msg::ERROR << "grid size too small, cannot resample";
    return false;
  }
  
  real_t dmin, dbuf=0;
  
  list< pair< vector<real_t>, real_t> > output_cache;
  list< pair< vector<real_t>, real_t> > local_cache;
  list< pair< vector<real_t>, real_t> >::iterator iitr;
  list< pair< vector<real_t>, real_t> >::iterator jitr;
  list< pair< vector<real_t>, real_t> >::iterator jitr_min;
  
  vector<real_t> apoint(_ndim, 0);
  vector<real_t> bpoint(_ndim, 0);
  
  for(integer_t ipt=0; ipt<_npoints; ++ipt) 
  {
    std::copy(_data[ipt], _data[ipt] + _ndim, apoint.begin());
    local_cache.push_back(pair< vector<real_t>, real_t>(apoint, _wgts[ipt]));
  }
  
  logger::log() << msg::DEBUG << "created local cache of size " << local_cache.size();
  
  for(unsigned isample=0; isample < nsample; ++isample)
  {
    output_cache.clear();
    
    iitr = local_cache.begin();
    jitr = local_cache.begin();
    
    integer_t ipt = 0;
    while(iitr != local_cache.end())
    {
      dmin = -1;
      for(jitr = local_cache.begin(); jitr != local_cache.end(); ++jitr)
      {
	if(jitr == iitr) 
	  continue;
        dbuf = agf::metric<real_t>(&((*iitr).first[0]), &((*jitr).first[0]), _ndim);
	if(dbuf < dmin || dmin < 0)
	{
	  dmin = dbuf;
	  jitr_min = jitr;
	}
      } // nested loop over grid points
      if((dmax <= 0 || dmin < dmax) && dmin > 0)
      {
	apoint = (*iitr).first;
	bpoint = (*jitr_min).first;
	
	double weight = (*iitr).second + (*jitr_min).second;
	
	jitr = local_cache.erase(jitr_min);
	iitr = local_cache.erase(iitr);
	
	for(unsigned int ix = 0; ix < (unsigned)_ndim; ++ix)
	{
	  apoint[ix] = (apoint[ix] + bpoint[ix]) / 2.0;
	}
	
	output_cache.push_back(pair<vector<real_t>, real_t>(apoint, weight));
      }
      else
      {
	++iitr;
      }
      
      if((++ipt) % (_npoints / 20) == 0)
	logger::log() << msg::DEBUG << "processed " << (real_t(ipt)/_npoints)*200 << "% of grid points";
      
    } // loop over grid points
    
    std::copy(output_cache.begin(), output_cache.end(), std::back_inserter(local_cache));
    
    logger::log() << msg::INFO << "resampling iteration " << isample << ", new grid population is " << local_cache.size();
    
  } // loop over re-sampling iterations
  
  _npoints = local_cache.size();
  
  release<real_array_t>(_data, _nreserved);
  release<real_t>(_wgts);
  
  _nreserved = _npoints + INIT_RESERVE;
  
  _data = new real_array_t[_nreserved];
  _wgts = new real_t[_nreserved];
  
  integer_t ipt = 0;
  for(iitr = local_cache.begin(); iitr != local_cache.end(); ++iitr)
  {
    _wgts[ipt] = iitr->second;
    _data[ipt] = new real_t[_ndim];
    std::copy(iitr->first.begin(), iitr->first.end(), (_data[ipt++]));
  } // reconstruct the grid
  
  set_metadata();
  
  return true;
}

//////////////////////////////////////////////////////////////////////

bool megrid::apply_transformation(transform::itransformation_base* t, bool retain, bool force_configure)
{
  if(t->configure(this, force_configure) &&
     t->get_out_dim() != 0                 &&
     (integer_t)t->get_in_dim() == get_ndim())
  {
    real_array_t* local_data = new real_array_t[_nreserved];
    for(integer_t ipt=0; ipt<_nreserved; ++ipt)
    {
      local_data[ipt] = new real_t[t->get_out_dim()];
      if(ipt<_npoints)
      {
	(*t)(_data[ipt], local_data[ipt]);
	if(logger::msg_level >= msg::DEBUG)
	{
	  if(_npoints <= 10 || ipt % (_npoints / 10) == 0)
	  {
	    logger::log() << msg::DEBUG << "transformed [" << ipt << "] of [" << _npoints << "]";
	    std::stringstream sstr;
	    for(integer_t idim=0; idim<get_ndim(); ++idim)
	    {
	      sstr << _data[ipt][idim] << "->" << local_data[ipt][idim] << " ";
	    }
	    logger::log() << msg::DEBUG << sstr.str();
	  }
	}
      }
      delete[] _data[ipt];
    }
    delete[] _data; // release<real_array_t>(_data, _nreserved);
    _data = local_data;    
    
    logger::log() << msg::INFO << "applied transformation to grid with dimension [" 
		  << get_ndim() << "], final grid has dimension [" << t->get_out_dim() << "]";
    
    _ndim = t->get_out_dim();
    
    set_metadata();
    
    if(retain)
    {
      _transformations.push_back(t);
    }
    _is_locked = true;
    return true;
  }
  else
  {
    logger::log() << msg::WARN << "failed to configure transformation";
    return false;
  }
}

//////////////////////////////////////////////////////////////////////

std::vector<const transform::itransformation_base*> megrid::get_transformations() const
{
  std::vector<const transform::itransformation_base*> buff;
  for(unsigned int it=0; it<_transformations.size(); ++it)
  {
    buff.push_back(_transformations[it]);
  }
  return buff;
}

//////////////////////////////////////////////////////////////////////

bool megrid::store(const string& file_name,
		   const string& tree_name,
		   const string& branch_names /* bool cached */)
{
  if(_npoints == 0)
  {
    logger::log() << msg::ERROR << "grid not defined, not saving to output file";
    return false;
  }
  
  TFile* fptr = new TFile(file_name.c_str(), "recreate");
  if(fptr->IsZombie())
  {
    logger::log() << msg::ERROR << "could not open file [" << file_name << "] for reading" ;
    delete fptr;
    return false;
  }
  
  TTree* tptr = new TTree(tree_name.c_str(), "");
  
  vector<string> br_names(_ndim, "");
  vector<real_t> br_refs(_ndim, 0.0);
  
  real_t weight = 0.0;
  
  boost::char_separator<char> sep("|");
  boost::tokenizer< boost::char_separator<char> > tokens(branch_names, sep);
  std::copy(tokens.begin(), tokens.end(), br_names.begin());
  
  for(integer_t ibr = 0; ibr < _ndim; ++ibr)
  {
    tptr->Branch(br_names[ibr].c_str(), &(br_refs[ibr]), (br_names[ibr] + "/F").c_str());
  }
  tptr->Branch("weight", &weight, "weight/F");
  
  real_array_t xarr[1];
  xarr[0] = new real_t[_ndim];
  
  for(integer_t ipt = 0; ipt < _npoints; ++ipt)
  {
    weight = _wgts[ipt];
    
    //
    // @FIXME get un-transformed points ... (?)
    //
    
    std::copy(_data[ipt], _data[ipt] + _ndim, br_refs.begin());
    tptr->Fill();
  }
  
  delete[] xarr[0];
  
  tptr->Write();
  fptr->Close();
  
  gROOT->cd();
  
  return true;
}

//////////////////////////////////////////////////////////////////////

vector<real_t> megrid::get_transformed_point(const vector<real_t>& x) const
{
  
  // 
  // loop over transformations (in order) and apply them
  //
  
  vector<real_t> xa(x), xb(x);
  for(unsigned itransform=0; itransform<_transformations.size(); ++itransform)
  {
    xa = xb;
    xb = vector<real_t>(_transformations[itransform]->get_out_dim());
    (*(_transformations[itransform]))(&(xa[0]), &(xb[0]));
    
    if(logger::msg_level == msg::DEBUG)
    {
      std::stringstream sstra;
      std::stringstream sstrb;
      for(integer_t idim=0; idim<(integer_t)xa.size(); ++idim)
	sstra << xa[idim] << " ";
      for(integer_t idim=0; idim<(integer_t)xb.size(); ++idim)
	sstrb << xb[idim] << " ";
      logger::log() << msg::DEBUG << sstra.str();
      logger::log() << msg::DEBUG << sstrb.str();    
    }
  }
  
  return xb;
} 

//////////////////////////////////////////////////////////////////////

void megrid::get_transformed_point(const real_t* x, real_t* xp) const
{
  for(unsigned itransform=0; itransform<_transformations.size(); ++itransform)
  {
    (*(_transformations[itransform]))(x, xp);
  }
}

//////////////////////////////////////////////////////////////////////

foam* megrid::create_foam(real_t tail_cut,
			  real_t vol_frac,
			  unsigned ncells,
			  unsigned nmin,
			  unsigned nsampl,
			  unsigned nbin) const
{
  if(_npoints != 0)
    return new foam(_name, _data, _wgts, _ndim, _npoints, tail_cut, vol_frac, ncells, nmin, nsampl, nbin);
  return NULL;
}

//////////////////////////////////////////////////////////////////////

bool megrid::reserve(unsigned m)
{
  if(m == 0) return false;
  
  integer_t N = _nreserved + m;
  
  real_t* local_wgts = new real_t[N];
  real_array_t* local_data = new real_array_t[N];
  
  if(_wgts != 0x0)
  {
    std::copy(_wgts, _wgts + _npoints, local_wgts);
  }  
  for(integer_t ipt=0; ipt<N; ++ipt)
  {
    local_data[ipt] = new real_t[_ndim];
    if(_data != 0x0 && ipt < _npoints)
    {
      std::copy(_data[ipt], _data[ipt] + _ndim, local_data[ipt]);
    }
  }
  
  release<real_array_t>(_data, _nreserved);
  release<real_t>(_wgts);
  
  _nreserved = N;
  
  _data = local_data;
  _wgts = local_wgts;
  
  return true;
}

//////////////////////////////////////////////////////////////////////

bool megrid::add_point(const vector<real_t>& p, real_t w)
{
  if(_is_locked)
  {
    logger::log() << msg::ERROR << "grid is locked, cannot add point";
    return false;
  }
  
  if((integer_t)p.size() != get_ndim())
    return false; 
  
  if(_nreserved <= _npoints + 1)
  {
    reserve(max(INIT_RESERVE, _npoints/100)); 
  }
  
  std::copy(p.begin(), p.end(), _data[_npoints]); 
  _wgts[_npoints] = (real_t)w; 
  
  _npoints += 1;  
  
  set_metadata();
  
  return true;
}

//////////////////////////////////////////////////////////////////////

megrid& megrid::operator+=(const std::vector<real_t>& unweighted_pt)
{
  add_point(unweighted_pt);
  return (*this);
}

//////////////////////////////////////////////////////////////////////

megrid& megrid::operator+=(const std::pair<std::vector<real_t>, real_t>& weighted_pt)
{
  add_point(weighted_pt.first, weighted_pt.second);
  return (*this);
}

//////////////////////////////////////////////////////////////////////

bool megrid::add_grid(const megrid& g)
{
  if(g.get_ndim() != _ndim)
  {
    logger::log() << msg::ERROR << "cannot add grids with different dimensions";
    return false;
  }
  
  //
  // check that all transformations are equivalent
  //  
  std::vector<const transform::itransformation_base*> tbuff = g.get_transformations();
  
  if(_transformations.size() != tbuff.size())
  {
    logger::log() << msg::ERROR << "cannot add grids with different coordinate transformations applied";
    return false;
  }
  for(unsigned int it=0; it<_transformations.size(); ++it)
  {
    if(!((_transformations[it])->operator==(tbuff[it])))
    {
      logger::log() << msg::ERROR << "cannot add grids with different coordinate transformations applied";
      return false;
    }
  }  
  
  //
  // add grid data and weights
  //
  if(_nreserved < _npoints + g.get_npoints())
  {
    integer_t nadd = (_npoints + g.get_npoints()) - _nreserved;
    reserve(nadd + INIT_RESERVE);
  }
  for(integer_t ipt=0; ipt<g.get_npoints(); ++ipt)
  {
    std::copy(g.get_data()[ipt], g.get_data()[ipt] + _ndim, _data[ipt+_npoints]);
    _wgts[ipt+_npoints] = g.get_weights()[ipt];
  }
  _npoints += g.get_npoints();
  
  //
  // re-configure all transformations
  //
  for(unsigned int it=0; it<_transformations.size(); ++it)
  {
    _transformations[it]->configure(this);
  }
  
  set_metadata();
  
  return true;
}

//////////////////////////////////////////////////////////////////////

megrid& megrid::operator+=(const megrid& g)
{
  add_grid(g);
  return (*this);
}

//////////////////////////////////////////////////////////////////////

bool megrid::cluster(cluster_type t, unsigned int nclusters, char method, char metric)
{
  switch(t)
  {
    case KMEANS:
      return kmeans_cluster(nclusters, method, metric);
      break;
    case HIERARCHICAL:
      return hierarchical_cluster(nclusters, method, metric);
      break;
    default:
      return false;
  }
}

//////////////////////////////////////////////////////////////////////

bool megrid::kmeans_cluster(unsigned int nclusters, char method, char metric, unsigned int npass)
{  
  if(_npoints < (integer_t)nclusters)
  {
    logger::log() << msg::ERROR << "grid size too small, cannot cluster";
    return false;
  }
  
  logger::log() << msg::INFO << "begin clustering using k-means";  
  
  std::vector<integer_t> mask_buff(_ndim, 0);
  std::vector<integer_array_t> masks(_npoints, NULL);
  for(integer_t ipt=0; ipt < _npoints; ++ipt)
  {
    masks[ipt] = &(mask_buff[0]);
  }
  std::vector<integer_t> clids(_npoints, 0);
  real_t clerr;
  integer_t clfound;
  kcluster(nclusters, _npoints, _ndim, _data, &(masks[0]), _wgts, false, npass, method, metric, &(clids[0]), &clerr, &clfound);
  
  if(clfound <= 0)
  {
    logger::log() << msg::ERROR << "k-means clustering failed with status code [" << clfound << "]";
    return false;
  }
  
  real_array_t* local_data = new real_array_t[nclusters];
  real_t* local_wgts = new real_t[nclusters];
  for(integer_t icl=0; icl<(integer_t)nclusters; ++icl)
  {
    local_wgts[icl] = 0.;
    local_data[icl] = new real_t[_ndim];
  }
  
  std::cout << "clusters: ";
  for(integer_t ipt=0; ipt<_npoints; ++ipt) 
  {
    integer_t icl = clids[ipt]; // index of the cluster that this grid point belongs to
    for(integer_t idim=0; idim<_ndim; ++idim)
    {
      local_data[icl][idim] = (((local_wgts[icl] * local_data[icl][idim]) + (_wgts[ipt] * _data[ipt][idim])) / 
			       (_wgts[ipt] + local_wgts[icl])); // weighted centroid of merged cluster
    }
    std::cout << icl << " ";
    local_wgts[icl] += _wgts[ipt];
  }
  std::cout << std::endl;
  
  logger::log() << msg::INFO << "finished clustering using k-means";  
  logger::log() << msg::DEBUG << "overwrite existing grid";
  
  release<real_array_t>(_data, _nreserved);
  release<real_t>(_wgts);
  
  _npoints = nclusters;
  _nreserved = nclusters;
  
  _data = local_data;
  _wgts = local_wgts;
  
  reserve(INIT_RESERVE); // add padding s.t. adding new points is fast ... 
  
  set_metadata();
  
  return true;
}

//////////////////////////////////////////////////////////////////////

bool megrid::hierarchical_cluster(unsigned int nclusters, char method, char metric)
{  
  if(_npoints < (integer_t)nclusters)
  {
    logger::log() << msg::ERROR << "grid size too small, cannot cluster";
    return false;
  }
  
  logger::log() << msg::INFO << "begin clustering using tree-clustering";  
  
  std::vector<integer_t> mask_buff(_ndim, 0);
  std::vector<integer_array_t> masks(_npoints, NULL);
  for(integer_t ipt=0; ipt < _npoints; ++ipt)
  {
    masks[ipt] = &(mask_buff.at(0));
  }
  
  Node* root_node = treecluster(_npoints, _ndim, _data, &(masks[0]), _wgts, false, metric, method, NULL);
  
  if(root_node == NULL)
  {
    logger::log() << msg::ERROR << "tree-clustering failed (memory allocation error)";  
    return false;
  }
  
  std::vector<integer_t> clids(_npoints, 0);
  
  cuttree(_npoints, root_node, nclusters, &(clids[0]));
  
  delete[] root_node;
  
  real_array_t* local_data = new real_array_t[nclusters];
  real_t* local_wgts = new real_t[nclusters];
  for(integer_t icl=0; icl<(integer_t)nclusters; ++icl)
  {
    local_wgts[icl] = 0.;
    local_data[icl] = new real_t[_ndim];
  }
  
  std::cout << "clusters: ";
  for(integer_t ipt=0; ipt<_npoints; ++ipt) // 
  {
    integer_t icl = clids[ipt];
    for(integer_t idim=0; idim<_ndim; ++idim)
    {
      local_data[icl][idim] = (((local_wgts[icl] * local_data[icl][idim]) + (_wgts[ipt] * _data[ipt][idim])) / 
			       (_wgts[ipt] + local_wgts[icl])); // weighted centroid of merged cluster
    }
    std::cout << icl << " ";
    local_wgts[icl] += _wgts[ipt];
  }
  std::cout << std::endl;
  
  logger::log() << msg::INFO << "finished tree-clustering";  
  logger::log() << msg::DEBUG << "overwrite existing grid";
  
  release<real_array_t>(_data, _nreserved);
  release<real_t>(_wgts);
  
  _npoints = nclusters;
  _nreserved = nclusters;
  
  _data = local_data;
  _wgts = local_wgts;
  
  reserve(INIT_RESERVE); // add padding s.t. adding new points is fast ... 
  
  set_metadata();
  
  return true;
}

//////////////////////////////////////////////////////////////////////

hist1D* megrid::project1D(integer_t ix, integer_t nbins, real_t a, real_t b, bool dyn) const 
{
  if(ix >= _ndim || ix < 0)
    return NULL;
  hist1D* h = new hist1D(Form("h_proj_%d", ix), Form("d_{%d}", ix), nbins, a, b);
  if(dyn)
    h->SetBit(TH1::kCanRebin);
  for(integer_t ipt=0; ipt<_npoints; ++ipt) // 
  {
    h->Fill(_data[ipt][ix], _wgts[ipt]);
  }
  return h;
}

//////////////////////////////////////////////////////////////////////

hist2D* megrid::project2D(integer_t ix,  integer_t iy, integer_t nbins, 
			  real_t ax, real_t bx, real_t ay, real_t by, bool dyn) const 
{
  if(ix >= _ndim || ix < 0 || 
     iy >= _ndim || iy < 0)
    return NULL;
  hist2D* h = new hist2D(Form("h_proj_%d", ix), Form("d_{%d}", ix), nbins, ax, bx, nbins, ay, by);
  if(dyn)
    h->SetBit(TH1::kCanRebin);
  for(integer_t ipt=0; ipt<_npoints; ++ipt) // 
  {
    h->Fill(_data[ipt][ix], _data[ipt][iy], _wgts[ipt]);
  }
  return h;
}

//////////////////////////////////////////////////////////////////////

void megrid::scan(integer_t istart, integer_t nrows, const char* delim) const
{
  if(istart >= _npoints) return;
  integer_t iend = min(istart + nrows, _npoints);
  
  std::streamsize p = std::cout.precision();
  std::streamsize w = std::cout.width();
  
  std::cout << std::setprecision(3);
  for(integer_t idim=0; idim<_ndim; ++idim)
  {
    std::cout << std::setw(12) << idim << delim;
  }
  std::cout << std::setw(12) << "weight" << delim;
  std::cout << std::endl;
  for(integer_t ipt=istart; ipt<iend; ++ipt)
  {
    for(integer_t idim=0; idim<_ndim; ++idim)
    {
      std::cout << std::setw(12) << std::scientific << _data[ipt][idim] << delim;
    }
    std::cout << std::setw(12) << std::scientific << _wgts[ipt] << delim;
    std::cout << std::endl;
  }
  std::cout << std::setw(w) << std::setprecision(p) << std::fixed;
}

//////////////////////////////////////////////////////////////////////

void megrid::set_loglevel(const std::string& lvl)
{ 
  if(lvl == "DEBUG")
    logger::set_level(log4cpp::Priority::DEBUG);
  if(lvl == "INFO")
    logger::set_level(log4cpp::Priority::INFO);
  if(lvl == "WARNING")
    logger::set_level(log4cpp::Priority::WARN);
  if(lvl == "ERROR")
    logger::set_level(log4cpp::Priority::ERROR);
}

//////////////////////////////////////////////////////////////////////
void megrid::set_metadata() 
{
  _grid_dimensions = std::vector<real_t>(_ndim, 1);
  for(integer_t idim=0; idim<_ndim; ++idim)
  {
    agf::pdf::grid_helper<real_t> h(_data, idim, _npoints);
    _grid_dimensions[idim] = std::fabs(h.get_max() - h.get_min());
  }
  _sum_wgts = std::accumulate(_wgts, _wgts + _npoints, 0.);
}
