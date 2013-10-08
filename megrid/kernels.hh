
#include "TMVA/PDEFoamKernelGauss.h"
#include "TMVA/PDEFoamKernelLinN.h" 

/* 
   
 *
 * @TODO ...
 * 

class PDEFoamKernelAGF : public TMVA::PDEFoamKernelBase
{
public:
  PDEFoamKernelAGF(Float_t wc);    
  PDEFoamKernelAGF(const PDEFoamKernelAGF&); 
  virtual ~PDEFoamKernelAGF() { };  
  
  // kernel estimator
  virtual Float_t Estimate(PDEFoam*, std::vector<Float_t>&, ECellValue);
  
  ClassDef(PDEFoamKernelAGF, 1); // AGF style PDEFoam kernel estimator
protected:
  Float_t fWC; // weight goal
  
  // Square function (fastest implementation)
  template<typename T> T Sqr(T x) const { return x * x; }
  
  // calculate gaussian weight
  Float_t WeightGaus(PDEFoam*, PDEFoamCell*, std::vector<Float_t>&);
  
  // estimate the cell value by its neighbors
  Float_t GetAverageNeighborsValue(PDEFoam*, std::vector<Float_t>&, ECellValue);
}; 

*/
