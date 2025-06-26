#ifndef EIKONAL_H

#define EIKONAL_H
#include "Kernel.h"
	
//#include "ElasticMaterial.h"

class Eikonal;
	
template <>

InputParameters validParams<Eikonal>();

class Eikonal : public Kernel
	
{

public:
	
  Eikonal(const InputParameters &parameters);

//  virtual ~MonodomainDiffusion();
	
protected:

  virtual Real computeQpResidual();
	
  virtual Real computeQpJacobian();
    
    Real _tau;
    
    Real _c0;
	
    //Real _theta;

    RealTensorValue _sigma;

    const MaterialProperty<RealTensorValue> &_conductivity;
	
//  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
	
  // unsigned const _numComp;

//  RealVectorValue zeros;

//  RealTensorValue H;

//private:
	
//  Real _surface_to_volume;
	
//
	
};
	
#endif
