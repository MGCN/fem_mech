#pragma once
#include "Material.h"
	
class MonodomainMaterialEllipsoid;
	
template <>
InputParameters validParams<MonodomainMaterialEllipsoid>();

class MonodomainMaterialEllipsoid : public Material
{
public:
    static InputParameters validParams();
    
    MonodomainMaterialEllipsoid(const InputParameters &parameters);
	
protected:

  virtual void computeQpProperties();

protected:

  const MaterialProperty<RealTensorValue> &_fXf, &_sXs, &_nXn;

    RealVectorValue const sigma_i;
    
    RealVectorValue const sigma_e;
    
    Real const C_m;
    Real const Chi;
   
    //RealVectorValue const _mono;

  MaterialProperty<RealTensorValue> &_Kmono;
    
  MaterialProperty<Real> &time_coefficient;
};
