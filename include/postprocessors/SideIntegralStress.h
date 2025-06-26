#ifndef SIDEINTEGRALSTRESS_H
#define SIDEINTEGRALSTRESS_H

#include "SideIntegralPostprocessor.h"
#include "MooseVariableInterface.h"
#include "ElasticMaterial.h"
// Forward Declarations
class SideIntegralStress;


template <>
InputParameters validParams<SideIntegralStress>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class SideIntegralStress : public SideIntegralPostprocessor
{
public:
  SideIntegralStress(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  const MaterialProperty<RealTensorValue> &_stress;
};

#endif

