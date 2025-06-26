#ifndef MONODOMAINCONDUCTIVITY_H
#define MONODOMAINCONDUCTIVITY_H

#include "Material.h"

class MonodomainConductivity;

template <>
InputParameters validParams<MonodomainConductivity>();

/**
 * Material for supplying a conductivity tensor.
 * \f$G=\sum_{i=f,s,n}\sigma_i(\vec{e}_i\otimes\vec{e}_i)\f$
 * Conductivities in fibre, sheet and sheet-normal directions are read from the
 * input file.
 * For sensible values, see e.g. \ref Potse2006 "Potse, 2006, Table I".
 */
class MonodomainConductivity : public Material
{
public:
  MonodomainConductivity(const InputParameters &parameters);

protected:
  virtual void computeQpProperties();

protected:



  const MaterialProperty<RealTensorValue> &_fXf, &_sXs, &_nXn;


  std::vector<Real> _intraConductivities;
  std::vector<Real> _extraConductivities;
  std::vector<Real> _monoConductivities;

  MaterialProperty<RealTensorValue> &_Gmono;
  MaterialProperty<RealTensorValue> &_Gintra;

};

#endif
