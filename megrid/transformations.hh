#ifndef transformations_hh
#define transformations_hh

#include "types.hh"

#include <vector>
#include <string>
#include <set>

class igrid_base;

namespace transform 
{
  
  //////////////////////////////////////////////////////////////////////////////////
  
  class itransformation_base
  {
  public:
    itransformation_base(unsigned in_dim=0, unsigned out_dim=0) : 
      _configured(false),
      _in_dim(in_dim),
      _out_dim(out_dim) { }
    virtual ~itransformation_base() { }
    
    virtual real_t* operator()(const real_t* x) const = 0;
    virtual void operator()(const real_t* x, real_t* xp) const = 0;
    
    virtual bool configure(const igrid_base*, bool force=false) = 0;

    virtual bool operator==(const itransformation_base* t) const = 0;
    
    virtual unsigned get_in_dim() const { return _in_dim; }
    virtual unsigned get_out_dim() const { return _out_dim; }
    
  protected:
    void set_in_dim(unsigned n) { _in_dim = n; }
    void set_out_dim(unsigned n) { _out_dim = n; }
    
    bool _configured;
    
  private:
    unsigned _in_dim;
    unsigned _out_dim;
  };
  
  //////////////////////////////////////////////////////////////////////////////////
  
  class rotate_phi : public itransformation_base
  {
  public:
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
    rotate_phi(std::string coord_def) : 
      itransformation_base(),
      _coord_def(coord_def) { }
    virtual ~rotate_phi() { }
    
    virtual real_t* operator()(const real_t* x) const;
    virtual void operator()(const real_t* x, real_t* xp) const;
    
    virtual bool configure(const igrid_base*, bool force=false);

    virtual bool operator==(const itransformation_base* t) const;

    const std::vector< std::pair<unsigned,unsigned> >& get_vector_indices() const { return _i_pairs; }

  private:
    std::string _coord_def;
    std::vector< std::pair<unsigned,unsigned> > _i_pairs; 
    std::set<unsigned> _paired_indexes;
  };
  
  //////////////////////////////////////////////////////////////////////////////////
  
  class norm_gauss : public itransformation_base
  {
  public:
    norm_gauss(unsigned idim) : 
      itransformation_base(idim, idim),
      _std_deviations(idim, 0), 
      _means(idim, 0) { }
    virtual ~norm_gauss() { }
    
    virtual real_t* operator()(const real_t* x) const;
    virtual void operator()(const real_t* x, real_t* xp) const;
    
    virtual bool configure(const igrid_base*, bool force=false);

    virtual bool operator==(const itransformation_base* t) const;

    const std::vector<real_t>& get_std_deviations() const { return _std_deviations; }
    const std::vector<real_t>& get_means() const { return _means; }
    
  private:
    std::vector<real_t> _std_deviations;
    std::vector<real_t> _means;
  };
  
  //////////////////////////////////////////////////////////////////////////////////
  
  class norm_range : public itransformation_base
  {
  public:
    norm_range(unsigned idim, unsigned sampling_frequency = 1) : 
      itransformation_base(idim, idim),
      _sampling_frequency(sampling_frequency) { }
    virtual ~norm_range() { }
    
    virtual real_t* operator()(const real_t* x) const;
    virtual void operator()(const real_t* x, real_t* xp) const;
    
    virtual bool configure(const igrid_base*, bool force=false);

    virtual bool operator==(const itransformation_base* t) const;

    const std::vector< std::vector<real_t> >& get_cdf() const { return _sorted_data; }
    unsigned get_sampling_frequency() const { return _sampling_frequency; }
    
  private:
    unsigned _sampling_frequency;
    std::vector< std::vector<real_t> > _sorted_data;
  };

  //////////////////////////////////////////////////////////////////////////////////
  
  class scale_fixed : public itransformation_base
  {
  public:
    scale_fixed(const real_t* scales, integer_t ndim) :
      itransformation_base(ndim, ndim),
      _scales(ndim,1.0) { std::copy(scales, scales+ndim, _scales.begin()); _configured=true; }
    scale_fixed(const std::vector<real_t>& scales) : 
      itransformation_base(scales.size(), scales.size()),
      _scales(scales) { }
    virtual ~scale_fixed() { }
    
    virtual real_t* operator()(const real_t* x) const;
    virtual void operator()(const real_t* x, real_t* xp) const;
    
    virtual bool configure(const igrid_base*, bool force=false);

    virtual bool operator==(const itransformation_base* t) const;

    const std::vector<real_t> get_scale_factors() const { return _scales; }
    
  private:
    std::vector<real_t> _scales;
  };

}

#endif
