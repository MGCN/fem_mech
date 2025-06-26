#include "GrowthTrace.h"
#include "FEProblem.h"
#include "libmesh/quadrature_gauss.h"
#include <limits>
#include "GrowthMaterial.h"

registerMooseObject("MECHApp", GrowthTrace);

template <>


InputParameters
validParams<GrowthTrace>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

GrowthTrace::GrowthTrace(const InputParameters &parameters)
    :  AuxKernel(parameters),
      _Se(getMaterialProperty<RealTensorValue>("SecondPiolaElastic")),
      _Ce(getMaterialProperty<RealTensorValue>("CauchyGreenElastic"))

      
{
  
}


Real
GrowthTrace::computeValue()
{
    Real result=0;  

    result = _Se[_qp](0,0) * _Ce[_qp](0,0) + _Se[_qp](1,1) * _Ce[_qp](1,1)  + _Se[_qp](2,2) * _Ce[_qp](2,2) ;
    
    return result;

}



