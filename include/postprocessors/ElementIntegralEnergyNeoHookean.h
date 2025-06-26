/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef ElementIntegralEnergyNeoHookean_H
#define ElementIntegralEnergyNeoHookean_H

#include "ElementIntegralPostprocessor.h"
#include "ElasticMaterial.h"
#include "IncompressibleNeoHookean.h"


// Forward Declarations
class ElementIntegralEnergyNeoHookean;

template <>
InputParameters validParams<ElementIntegralEnergyNeoHookean>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class ElementIntegralEnergyNeoHookean : public ElementIntegralPostprocessor
{
public:
  ElementIntegralEnergyNeoHookean(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  Real _mu;
  Real _kappa;


  MaterialProperty<Real> const & _J;
 
  MaterialProperty<RealTensorValue> const &_C;
 
  MaterialProperty<Real> const & _J23;
    
  MaterialProperty<RealTensorValue> const & _F;
    
  const VariableValue & _pressure;


};

#endif
